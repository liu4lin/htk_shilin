/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/* developed at:                                               */
/*                                                             */
/*      Speech Vision and Robotics group                       */
/*      Cambridge University Engineering Department            */
/*      http://svr-www.eng.cam.ac.uk/                          */
/*                                                             */
/*      Entropic Cambridge Research Laboratory                 */
/*      (now part of Microsoft)                                */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright: Microsoft Corporation                    */
/*          1995-2000 Redmond, Washington USA                  */
/*                    http://www.microsoft.com                 */
/*                                                             */
/*              2001  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*  File: HCompV.c: HMM global mean/variance initialisation    */
/* ----------------------------------------------------------- */

char *hcompv_version = "!HVER!HCompV:   3.4.1 [CUED 12/03/09]";
char *hcompv_vc_id = "$Id: HCompV.c,v 1.1.1.1 2006/10/11 09:54:59 jal58 Exp $";


/* 
   This program calculates a single overall variance vector from a
   sequence of training files.  This vector is then copied into
   all of the components of all the states of the hidden Markov
   model specified on the command line.  Optionally, HCompV
   will also update the mean vector. HCompV.c can be used as
   the initial step in Fixed Variance and Grand Variance training schemes,
   and to initialise HMMs for "flat-start" training.  
   The structure of the HMM to be initialised must be
   defined beforehand (eg using a prototype HMM def).  

   It can also be used to generate variance floor vectors.
*/ 

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HAudio.h"
#include "HWave.h"
#include "HVQ.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HUtil.h"
#include "pthread.h"
/* -------------------------- Trace Flags & Vars ------------------------ */

#define T_TOP     0001           /* basic progress reporting */
#define T_COVS    0002           /* show covariance matrices */
#define T_LOAD    0004           /* trace data loading */
#define T_SEGS    0010           /* list label segments */
#define T_CMV     0012           /* doing CMV */
#define BIG_FLOAT 1E20      /* Limit for float arguments */
static int  trace    = 0;           /* trace level */
static ConfParam *cParm[MAXGLOBS];   /* configuration parameters */
static int nParm = 0;               /* total num params */

/* ---------------------- Global Variables ----------------------- */

static char *segLab = NULL;         /* segment label if any */
static LabId  segId  = NULL;        /* and its id */
static LabId  hmmId  = NULL;        /* id of model */
static char *labDir = NULL;         /* label file directory */
static char *labExt = "lab";        /* label file extension */
static float minVar  = 0.0;         /* minimum variance */
static FileFormat dff=UNDEFF;       /* data file format */
static FileFormat lff=UNDEFF;       /* label file format */
static char *hmmfn=NULL;            /* HMM definition file name */
static char *outfn=NULL;            /* output HMM file name (name only) */
static char *outDir=NULL;           /* HMM output directory */
static long totalCount=0;           /* total number of vector samples*/
static Boolean meanUpdate = FALSE;  /* update means  */
static Boolean saveBinary = FALSE;  /* save output in binary  */
static float vFloorScale = 0.0;     /* if >0.0 then vFloor scaling */

/* Major Data Structures */
static MLink macroLink;             /* Link to specific HMM macro */
static HLink hmmLink;               /* Link to the physical HMM */
static HMMSet hset;                 /* HMM to be initialised with */

static MemHeap iStack;

static int xpwin=1;
static Boolean onlyStatic=FALSE;
static MemHeap hmmStack;   /*For Storage of all dynamic structures created...*/
HMMSet synrhset;
char gmmclusterListFn[MAXFNAMELEN];
Boolean gotgmmcluster=FALSE;
HLink* gmmcluster=NULL;

/* Storage for mean and covariance accumulators */
typedef struct {
   Vector       meanSum;            /* acc for mean vector value */
   Covariance   squareSum;          /* acc for sum of squares */
   Covariance   fixed;              /* fixed (co)variance values */
} CovAcc;
static CovAcc accs[SMAX];           /* one CovAcc for each stream */
static Boolean fullcNeeded[SMAX];   /* true for each stream that needs full
                                       covariance calculated */
static Observation obs;             /* storage for observations  */

typedef struct _AccLoad{
	Vector tempV;
	Vector meanSum;
	TriMat squareSum;
	char** fileList;/*subset of filelist*/
	int numFiles;
	long numFrames;
}*AccLoad;

/* ------------ structures for cmn ------------ */

typedef struct{
   char SpkrName[MAXSTRLEN];             /* speaker name */
   int NumFrame;                         /* number of frames for speaker */
   Vector meanSum;                       /* mean accumulate structure for speaker */
   Vector squareSum;                     /* variance structure for speaker */
}SpkrAcc;                  

typedef struct SpkrAccListItem{
   SpkrAcc *sa;                          /* speaker accumulate */
   struct SpkrAccListItem *nextSpkr;     /* next pointer */
}SpkrAccListItem;

static SpkrAccListItem *salist = NULL;   /* global speaker accumulate list */
static int vSize = 0;                    /* target observation vector size */
static char spPattern[MAXSTRLEN];        /* speaker mask */
static char pathPattern[MAXSTRLEN];      /* path mask */
static char oflags[MAXSTRLEN] = "m";     /* export flags for CMV */  
static char cmDir[MAXSTRLEN];            /* directory to export CMV */
static char TargetPKStr[MAXSTRLEN];      /* target parm kind string */
static Boolean DoCMV = FALSE;            /* switch from old HCompV to CMV */
static Boolean DoTrim=FALSE;
static int trimNum=0;
static Boolean forceFc=FALSE;
static Boolean forceLog=FALSE;
static Boolean tiedCov=FALSE;

static char** fileList=NULL;
static char** fileList2=NULL;
static Boolean stereo=FALSE;
static int pca=0;

/* ------------- Process Command Line and Check Data ------------ */
Group* GMMClusterFile(char** fileList, int numFiles, int numClust, int vSize, int numThreads, Boolean fc, Boolean tiedCov);
Group* GMMClusterFileX(char** fileList, int numFiles, int numClust, int vSize, int numThreads, Boolean fc, Boolean tiedCov);
void AccParallel(char** fileList, int numFiles, int vSize, int numThreads);
/* SetConfParms: set conf parms relevant to HCompV  */
void CalcCovs(void);
void SetConfParms(void)
{
   Boolean b,c;
   int i;
   double d;
   
   nParm = GetConfig("HCOMPV", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfBool(cParm,nParm,"UPDATEMEANS",&b)) meanUpdate = b;
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&c)) saveBinary = c;
      if (GetConfFlt(cParm,nParm,"MINVARFLOOR",&d)) minVar = d;
   }
}

void ReportUsage(void)
{
   printf("\nUSAGE: HCompV [options] [hmmFile] trainFiles...\n" );
   printf(" Option                                       Default\n\n");
   printf(" -c dir  Set output directiry for CMV         none\n");
   printf(" -f f    Output vFloor as f * global var      none\n");
   printf(" -k s    spkr pattern for CMV                 none\n");
   printf(" -l s    Set segment label to s               none\n");
   printf(" -m      Update means                         off\n");
   printf(" -o fn   Store new hmm def in fn (name only)  outDir/srcfn\n");
   printf(" -p s    path pattern for CMV                 none\n");
   printf(" -q nmv  output type flags for CMV            m\n");
   printf(" -v f    Set minimum variance to f            0.0\n");
   printf(" -trim int Trim head&tail features to calc noise model 0\n");
   printf(" -fc     force generate full covariance       FALSE\n");
   printf(" -os     only calculate static coeff			 FALSE\n");
   printf(" -nt int number of threads to apply			 1\n");
   printf(" -nf int number of utterances(required)       none\n");
   printf(" -gmm int number of gmm clusters(required)     none\n");
   printf(" -E int  context expansion, window size       1\n");
   printf(" -W f    MMF as initialisation of gmm cluster none\n");
   printf(" -e f    MMF label list                       none\n");
   PrintStdOpts("BCFGHILMX");
   printf("\n\n");
}


/* ------------------------ Initialisation ----------------------- */

/* CheckVarianceKind: set fullcNeeded[s] for each non-diag stream s */
void CheckVarianceKind(void)
{
   int i,s,m;
   StateElem *se;
   StreamElem *ste;
   MixtureElem *me;
   
   for (s=1;s<=hset.swidth[0];s++)
      fullcNeeded[s]=FALSE;
   for (i=2,se=hmmLink->svec+2; i < hmmLink->numStates; i++,se++)
      for (s=1,ste=se->info->pdf+1; s <= hset.swidth[0]; s++,ste++)
         for (m=1,me = ste->spdf.cpdf+1; m<=ste->nMix; m++, me++)
            if (me->mpdf->ckind == FULLC) 
               fullcNeeded[s] = TRUE;
            else if(forceFc){
            	hset.ckind=FULLC;
            	fullcNeeded[s]=TRUE;
            	me->mpdf->ckind=FULLC;
            	me->mpdf->cov.inv=CreateSTriMat(hset.hmem,hset.vecSize);
            }
}

/* Initialise: load HMMs and create accumulators */
void Initialise(void)
{
   int s,V;
   Boolean eSep;
   char base[MAXSTRLEN];
   char path[MAXSTRLEN];
   char ext[MAXSTRLEN];

   /* Load HMM defs */     
   if(MakeOneHMM(&hset,BaseOf(hmmfn,base))<SUCCESS)
      HError(2028,"Initialise: MakeOneHMM failed");
   if(LoadHMMSet(&hset,PathOf(hmmfn,path),ExtnOf(hmmfn,ext))<SUCCESS)
      HError(2028,"Initialise: LoadHMMSet failed");
   SetParmHMMSet(&hset);
   if (hset.hsKind==DISCRETEHS || hset.hsKind==TIEDHS)
      HError(2030,"Initialise: HCompV only uses continuous models");

   /* Create a heap to store the input data */
   CreateHeap(&iStack,"InBuf", MSTAK, 1, 0.5, 100000, LONG_MAX);
   
   /* Get a pointer to the physical HMM */
   hmmId = GetLabId(base,FALSE);
   macroLink = FindMacroName(&hset,'h',hmmId);
   if (macroLink==NULL)
      HError(2020,"Initialise: cannot find hmm %s in hset",hmmfn);
   hmmLink = (HLink)macroLink->structure;

   /* Find out for which streams full covariance is needed */
   CheckVarianceKind( );

   /* Create accumulators for the mean and variance */
   for (s=1;s<=hset.swidth[0]; s++){
      V = hset.swidth[s];
      accs[s].meanSum=CreateVector(&gstack,V);
      ZeroVector(accs[s].meanSum);
      if (fullcNeeded[s]) {
         accs[s].squareSum.inv=CreateSTriMat(&gstack,V);
         accs[s].fixed.inv=CreateSTriMat(&gstack,V);
         ZeroTriMat(accs[s].squareSum.inv);
      }
      else {
         accs[s].squareSum.var=CreateSVector(&gstack,V);
         accs[s].fixed.var=CreateSVector(&gstack,V);
         ZeroVector(accs[s].squareSum.var);
      }
   }

   /* Create an object to hold the input parameters */
   SetStreamWidths(hset.pkind,hset.vecSize,hset.swidth,&eSep);
   obs=MakeObservation(&gstack,hset.swidth,hset.pkind,FALSE,eSep);
   if(segLab != NULL) {
      segId = GetLabId(segLab,TRUE);
   }

   if (trace&T_TOP) {
      printf("Calculating Fixed Variance\n");
      printf("  HMM Prototype: %s\n",hmmfn);
      printf("  Segment Label: %s\n",(segLab==NULL)?"None":segLab);
      printf("  Num Streams  : %d\n",hset.swidth[0]);
      printf("  UpdatingMeans: %s\n",(meanUpdate)?"Yes":"No");
      printf("  Target Direct: %s\n",(outDir==NULL)?"Current":outDir);     
   }
}

/* ----------------------[Co]Variance Estimation ---------------------- */

/* CalcCovs: calculate covariance of speech data */
void CalcCovs(void)
{
   int x,y,s,V;
   float meanx,meany,varxy,n;
   Matrix fullMat;
   
   if (totalCount<2)
      HError(2021,"CalcCovs: Only %d speech frames accumulated",totalCount);
   if (trace&T_TOP)
      printf("%ld speech frames accumulated\n", totalCount);
   n = (float)totalCount;     /* to prevent rounding to integer below */
   for (s=1; s<=hset.swidth[0]; s++){  /* For each stream   */
      V = hset.swidth[s];
      for (x=1; x<=V; x++)            /* For each coefficient ... */{
         accs[s].meanSum[x] /= n;         /* ... calculate mean */
      }
      for (x=1;x<=V;x++) {
         meanx = accs[s].meanSum[x];      /* ... and [co]variance */
         if (fullcNeeded[s]) {
            for (y=1; y<=x; y++) {
               meany = accs[s].meanSum[y];
               varxy = accs[s].squareSum.inv[x][y]/n - meanx*meany;
               accs[s].squareSum.inv[x][y] =
                  (x != y || varxy > minVar) ? varxy : minVar;    
            }
         }
         else {
            varxy = accs[s].squareSum.var[x]/n - meanx*meanx;
            accs[s].fixed.var[x] = (varxy > minVar) ? varxy :minVar;
         }
      }
      if (fullcNeeded[s]) { /* invert covariance matrix */
    	 if(pca){
    		 CopyTriMat(accs[s].squareSum.inv,accs[s].fixed.inv);
    	 }
    	 else{
    		 fullMat=CreateMatrix(&gstack,V,V);
    		 ZeroMatrix(fullMat);
    		 CovInvert(&gstack, accs[s].squareSum.inv,fullMat);
    		 Mat2Tri(fullMat,accs[s].fixed.inv);
    		 FreeMatrix(&gstack,fullMat);
    	 }
      }
      if (trace&T_COVS) {
         printf("Stream %d\n",s);
         if (meanUpdate)
            ShowVector(" Mean Vector ", accs[s].meanSum,12);
         if (fullcNeeded[s]) {
            ShowTriMat(" Covariance Matrix ",accs[s].squareSum.inv,12,12);
         } else
            ShowVector(" Variance Vector ", accs[s].fixed.var,12);
      }
   }
}

/* TriDiag2Vector: Copy diagonal from m into v */
void TriDiag2Vector(TriMat m, Vector v)
{
   int i,size;

   if (TriMatSize(m) != (size=VectorSize(v)))
      HError(2090,"TriDiag2Vector: Covariance sizes differ %d vs %d",
             TriMatSize(m),VectorSize(v));
   for (i=1; i<=size; i++)
      v[i] = m[i][i];
}

/* SetCovs: set covariance values in hmm */
void SetCovs(void)
{
   int i,s,m;
   StateElem *se;
   StreamElem *ste;
   MixtureElem *me;
   MixPDF *mp;

   CalcCovs();
   if (trace&T_TOP) {
      printf("Updating HMM ");
      if (meanUpdate) printf("Means and ");
      printf("Covariances\n");
   }
   for (i=2,se=hmmLink->svec+2; i < hmmLink->numStates; i++,se++)
      for (s=1,ste=se->info->pdf+1; s <= hset.swidth[0]; s++,ste++)
         for (m=1,me = ste->spdf.cpdf+1; m<=ste->nMix; m++, me++) {
            mp = me->mpdf;
            if (meanUpdate && !IsSeenV(mp->mean)){      /* meanSum now holds mean */
               CopyVector(accs[s].meanSum,mp->mean); 
               TouchV(mp->mean);
            }
            if (!IsSeenV(mp->cov.var)){
               if (mp->ckind==FULLC)
                  CopyMatrix(accs[s].fixed.inv,mp->cov.inv);
               else if (fullcNeeded[s])  /* dont need full cov, but its all we have */                
                  TriDiag2Vector(accs[s].fixed.inv,mp->cov.var);
               else
                  CopyVector(accs[s].fixed.var,mp->cov.var);
               TouchV(mp->cov.var);
            }
         }
   ClearSeenFlags(&hset,CLR_ALL);
}

/* PutVFloor: output variance floor vectors */
void PutVFloor(void)
{
   int i,s;
   char outfn[MAXSTRLEN],vName[32],num[10];
   FILE *f;
   Vector v;
   
   MakeFN("vFloors",outDir,NULL,outfn);
   if ((f = fopen(outfn,"w")) == NULL)
      HError(2011,"PutVFloor: cannot create %s",outfn);
   for (s=1; s <= hset.swidth[0]; s++) {
      v = CreateVector(&gstack,hset.swidth[s]);
      sprintf(num,"%d",s); 
      strcpy(vName,"varFloor"); strcat(vName,num);
      fprintf(f,"~v %s\n",vName);
      if (fullcNeeded[s])              
         TriDiag2Vector(accs[s].squareSum.inv,v);
      else
         CopyVector(accs[s].fixed.var,v);
      for (i=1; i<=hset.swidth[s]; i++)
         v[i] *= vFloorScale;
      fprintf(f,"<Variance> %d\n",hset.swidth[s]);
      WriteVector(f,v,FALSE);
      FreeVector(&gstack,v);
   }
   fclose(f);
   if (trace&T_TOP)
      printf("Var floor macros output to file %s\n",outfn);
}

/* ---------------- Load Data and Accumulate Stats --------------- */

/* AccVar:  update global accumulators with given observation */
void AccVar(Observation obs)
{
   int x,y,s,V;
   float valx,valy;
   Vector v;

   totalCount++;
   for (s=1; s<=hset.swidth[0]; s++){
      v = obs.fv[s]; V = hset.swidth[s];
      for (x=1;x<=V;x++) { 
         valx=forceLog?(v[x]>0?log(v[x]):LSMALLP):v[x];
         accs[s].meanSum[x] += valx;     /* accumulate mean */
         if (fullcNeeded[s]) {          /* accumulate covar */ 
            accs[s].squareSum.inv[x][x] += valx*valx;
            for (y=1;y<x;y++){
               valy=forceLog?(v[y]>0?log(v[y]):LSMALLP):v[y];
               accs[s].squareSum.inv[x][y] += valx*valy;
            }
         } else                         /* accumulate var */
            accs[s].squareSum.var[x] += valx*valx;
      }
   }
}

/* CheckData: check data file consistent with HMM definition */
void CheckData(char *fn, BufferInfo info) 
{
   if (info.tgtVecSize!=hset.vecSize)
      HError(2050,"CheckData: Vector size in %s[%d] is incompatible with hmm %s[%d]",
             fn,info.tgtVecSize,hmmfn,hset.vecSize);
   if (info.tgtPK != hset.pkind)
      HError(2050,"CheckData: Parameterisation in %s is incompatible with hmm %s",
             fn,hmmfn);
}

void LoadSpaf(char *spafn, int vSize){
	int idx,i,j,dnn_row,dnn_col;
	char ch;
	char dnn_buf[32];
	float value;
	FILE* dnn_f;
	if((dnn_f=fopen(spafn,"r"))==NULL){
		HError(1,"Cannot open file: %s",spafn);
	}
	fscanf(dnn_f,"%d %d\n",&dnn_row,&dnn_col);
	assert(dnn_col==vSize);
	for(i=1;i<=dnn_row;i++){
		for(idx=1;idx<=vSize;idx++){/*set all small posterior to be floor rather than zero to avoid LZERO*/
			obs.fv[1][idx]=SMALLP;
		}
		j=0;
		while((ch=getc(dnn_f))!=EOF){
			if(ch==' '||ch=='\n'){
				dnn_buf[j++]='\0';
				if(j>1){
					sscanf(dnn_buf,"%d:%e",&idx,&value);
					obs.fv[1][idx]=value;
				}
				j=0;
				if(ch=='\n')
					break;
			}
			else{
				dnn_buf[j++]=ch;
			}
		}
		AccVar(obs);
	}
	fclose(dnn_f);
}

/* LoadFile: load whole file or segments and accumulate variance */
void LoadFile(char *fn, Point plist, int *pntCount)
{
   ParmBuf pbuf;
   BufferInfo info;
   char labfn[80];
   Transcription *trans;
   long segStIdx,segEnIdx;  
   int i,j,ncas,nObs;
   LLink p;
   int vecSize=hset.swidth[1];
   if (segId == NULL)  {   /* load whole parameter file */
      if((pbuf=OpenBuffer(&iStack, fn, 0, dff, FALSE_dup, FALSE_dup))==NULL)
         HError(2050,"LoadFile: Config parameters invalid");
      GetBufferInfo(pbuf,&info);
      CheckData(fn,info);
      nObs = ObsInBuffer(pbuf);
      if(DoTrim){
          for (i=0; i<nObs; i++){
        	  if(i<trimNum||i>=nObs-trimNum){
                  ReadAsTable(pbuf,i,&obs);
                  AccVar(obs);
                  if(plist!=NULL){
                	  Vector vec=CreateVector(&gstack,vecSize);
                	  CopyVector(obs.fv[1],vec);
                	  Point pnt=CreatePoint(&gstack,vec,NULL,NULL);
                	  Point temp=plist->next;
                	  plist->next=pnt;
                	  pnt->next=temp;
                	  (*pntCount)++;
                  }
        	  }
          }
      }
      else{
          for (i=0; i<nObs; i++){
             ReadAsTable(pbuf,i,&obs);
             AccVar(obs);
             if(plist!=NULL){
           	  Vector vec=CreateVector(&gstack,vecSize);
           	  CopyVector(obs.fv[1],vec);
           	  Point pnt=CreatePoint(&gstack,vec,NULL,NULL);
           	  Point temp=plist->next;
           	  plist->next=pnt;
           	  pnt->next=temp;
           	  (*pntCount)++;
             }
          }
      }
      if (trace&T_LOAD) {
         printf(" %d observations loaded from %s\n",nObs,fn);
         fflush(stdout);
      }        
      CloseBuffer(pbuf);
   }
   else {                  /* load segment of parameter file */
      MakeFN(fn,labDir,labExt,labfn);
      trans = LOpen(&iStack,labfn,lff);
      ncas = NumCases(trans->head,segId);
      if ( ncas > 0) {
         if((pbuf=OpenBuffer(&iStack, fn, 0, dff, FALSE_dup, FALSE_dup))==NULL)
            HError(2050,"LoadFile: Config parameters invalid");
         GetBufferInfo(pbuf,&info);
         CheckData(fn,info);
         for (i=1,nObs=0; i<=ncas; i++) {
            p = GetCase(trans->head,segId,i);
            segStIdx= (long) (p->start/info.tgtSampRate);
            segEnIdx  = (long) (p->end/info.tgtSampRate);
            if (trace&T_SEGS)
               printf(" loading seg %s [%ld->%ld]\n",
                      segId->name,segStIdx,segEnIdx);
            if (segEnIdx >= ObsInBuffer(pbuf))
               segEnIdx = ObsInBuffer(pbuf)-1;
            if (segEnIdx >= segStIdx) {
               for (j=segStIdx;j<=segEnIdx;j++) {
                  ReadAsTable(pbuf,j,&obs);
                  AccVar(obs); ++nObs;
               }
            }
         }        
         CloseBuffer(pbuf);
         if (trace&T_LOAD)
            printf(" %d observations loaded from %s\n",nObs,fn);
      }  
   }
   ResetHeap(&iStack);
}

/* ------------------------- Save Model ----------------------- */

/* SaveModel: save HMMSet containing one model */
void SaveModel(char *outfn)
{
   if (outfn != NULL)
      macroLink->id = GetLabId(outfn,TRUE);
   if(SaveHMMSet(&hset,outDir,NULL,NULL,saveBinary)<SUCCESS)
      HError(2011,"SaveModel: SaveHMMSet failed");
}



/* ------------- Cepstral Mean Substraction & Variance Normalisation ------------ */

/* initialise and return an instance of SpkrAcc */  
SpkrAcc *InitSpkrAcc(void)
{
   SpkrAcc *sa;

   sa = New(&gstack,sizeof(SpkrAcc));
   sa->meanSum = CreateVector(&gstack,vSize);
   ZeroVector(sa->meanSum);
   sa->squareSum = CreateVector(&gstack,vSize);
   ZeroVector(sa->squareSum);
   sa->NumFrame = 0;
   return sa;
}

/* reset an instance of SpkrAcc type */ 
void ClrSpkrAcc(SpkrAcc *sa)
{
   ZeroVector(sa->meanSum);
   ZeroVector(sa->squareSum);
   sa->NumFrame = 0;
}
  
/* Accumulate stats from an utterance file */
SpkrAcc *AccGenUtt(char *SpkrPattern, char *UttFileName, SpkrAcc *sa)
{
   char SpkrName[MAXSTRLEN];
   ParmBuf pbuf;
   BufferInfo info;
   short swidth[SMAX];
   Boolean eSep;
   Vector tempV;
   int i;

   if (MaskMatch(SpkrPattern,SpkrName,UttFileName)==TRUE){
      /* open buffer and construct observation */
      pbuf = OpenBuffer(&iStack,UttFileName,0,dff,FALSE_dup,FALSE_dup);
      GetBufferInfo(pbuf,&info);
      if ((info.tgtPK & HASZEROM) && strchr(oflags,'m')) {
         HError(-2021,"HCompV: AccGenUtt: qualifier _Z not appropriate when calculating means!\n");
      }
      /* treat as single stream system though a bit weird */
      ZeroStreamWidths(1,swidth);
      SetStreamWidths(info.tgtPK,info.tgtVecSize,swidth,&eSep);
      obs = MakeObservation(&gstack,swidth,info.tgtPK,FALSE,eSep);
      if (info.tgtVecSize != vSize){
         vSize = info.tgtVecSize;
         /* if needed init a SpkrAcc */
         sa = InitSpkrAcc();
         fprintf(stdout,"Target observation vector size set to %d ......\n",info.tgtVecSize);
         fflush(stdout);
      }
      ParmKind2Str(info.tgtPK,TargetPKStr);
      /* accumulate stats for current utterance file */
      StartBuffer(pbuf);
      while (BufferStatus(pbuf) != PB_CLEARED) 
         {
            /* copy current observation and set vector ptr to first stream */
            ReadAsBuffer(pbuf,&obs);
            tempV = obs.fv[1];
            for (i=1;i<=vSize;i++){
               sa->meanSum[i] += tempV[i];
               sa->squareSum[i] += tempV[i]*tempV[i];
            }
            sa->NumFrame += 1;
         }
      CloseBuffer(pbuf);
      strcpy(sa->SpkrName,SpkrName);
      if (trace&T_CMV){
         fprintf(stdout,"Utterance %s accumulate generated for speaker %s\n",UttFileName,sa->SpkrName);
         fflush(stdout);
      }
      ResetHeap(&iStack);
      return sa;
   }
   else {
      HError(2039,"HCompV: AccGenUtt: speaker pattern matching failure on file: %s\n",UttFileName);
      return NULL;
   }
}


/* Append speaker accumulate structure list */
SpkrAccListItem *AppendSpkrAccList(SpkrAccListItem *sal, SpkrAcc *sa)
{
   SpkrAccListItem *temp;
   int i;

   temp = New(&gstack,sizeof(SpkrAccListItem));
   temp->sa = InitSpkrAcc();
   for (i=1;i<=vSize;i++){
      temp->sa->meanSum[i] = sa->meanSum[i];
      temp->sa->squareSum[i] = sa->squareSum[i];
   }
   temp->sa->NumFrame = sa->NumFrame;
   strcpy(temp->sa->SpkrName,sa->SpkrName);
   temp->nextSpkr = sal;
   sal = temp;
  
   if (trace&T_CMV){
      fprintf(stdout,"Creating entry for speaker %s ......\n",sa->SpkrName);
      fflush(stdout);
   }

   return sal;
}

/* 
   search the speaker accumulate list and if there is an entry for the
   current speaker given in sa, then update the entry in the list; otherwise
   insert sa into the list as a new speaker entry
*/
SpkrAccListItem *UpdateSpkrAccList(SpkrAccListItem *sal, SpkrAcc *sa)
{ 
   SpkrAccListItem *p;
   int i;
   Boolean tag = FALSE;

   p = sal;

   while (p != NULL){
      if (strcmp(p->sa->SpkrName,sa->SpkrName)==0){
         for (i=1;i<=vSize;i++){
            p->sa->meanSum[i] += sa->meanSum[i];
            p->sa->squareSum[i] += sa->squareSum[i];
         }
         p->sa->NumFrame += sa->NumFrame;
         tag = TRUE;
         break;
      }
      p = p->nextSpkr;
   }
   if (tag == FALSE){
      sal = AppendSpkrAccList(sal,sa);
   }

   return sal;
}


/* Compute Mean and Var */
void UpdateMeanVar(SpkrAccListItem *sal)
{
   int i;
   SpkrAccListItem *p;
  
   p = sal;
   while (p != NULL){
      for (i=1;i<=vSize;i++){
         p->sa->meanSum[i] /= ((float)(p->sa->NumFrame));
      }
      for (i=1;i<=vSize;i++){
         p->sa->squareSum[i] /= ((float)(p->sa->NumFrame));
         p->sa->squareSum[i] -= (p->sa->meanSum[i])*(p->sa->meanSum[i]);
      }
      p = p->nextSpkr;
   }  
}

/* Report output type NumFrame/cms/cvn */
void ReportOutput()
{
   if (strcmp(oflags,"m")==0){
      fprintf(stdout,"\n  HCompV: Comptuting side based cepstral mean ......\n\n");
   }
   else if (strcmp(oflags,"v")==0){
      fprintf(stdout,"\n  HCompV: Computing side based cepstral variance normaliser ......\n\n");
   }
   else if (strcmp(oflags,"mv")==0){
      fprintf(stdout,"\n  HCompV: Comptuting side based cepstral mean & variance normaliser ......\n\n");
   }
   else if (strcmp(oflags,"nv")==0){
      fprintf(stdout,"\n  HCompV: Comptuting side based frame number & variance normaliser ......\n\n");
   }
   else if (strcmp(oflags,"nmv")==0){
      fprintf(stdout,"\n  HCompV: Comptuting side based frame number & cepstral mean & variance normaliser ......\n\n");
   }
   else {
      fprintf(stdout,"\n  HCompV: ReportOutput: Unrecognisable output flag setting: %s\n",oflags);
      fflush(stdout);
      HError(2019,"HCompV: ReportOutput: Unrecognisable output flag setting: %s\n",oflags);
   }
   fflush(stdout);
}


/* Export number of frames, mean, var to a given directory */
void ExportNMV(SpkrAccListItem *sal, char *OutDirName, char *tgtPKStr) 
{
   FILE *oFile;
   Boolean isPipe;
   char oFileName[MAXSTRLEN];
   char pathBuffer1[MAXSTRLEN];
   char pathBuffer2[MAXSTRLEN];
   SpkrAccListItem *p;
   int i;

   p = sal;
   while(p != NULL){
      /* create output file name for current spkr index */    
      if ( pathPattern[0] != '\0'){
         if ( MaskMatch(pathPattern,pathBuffer1,p->sa->SpkrName) != TRUE ){
            HError(2039,"HCompV: ExportNMV: path pattern matching failure on speaker: %s\n",p->sa->SpkrName);
         }
         MakeFN(pathBuffer1,OutDirName,NULL,pathBuffer2); 
         MakeFN(p->sa->SpkrName,pathBuffer2,NULL,oFileName);
      }
      else
         MakeFN(p->sa->SpkrName,OutDirName,NULL,oFileName);

      /* open and write */
      oFile = FOpen(oFileName,NoOFilter,&isPipe);
      if (oFile == NULL){
         HError(2011,"HCompV: ExportNMV: output file creation error %s",oFileName);
      }
    
      /* write header */
      fprintf(oFile,"<CEPSNORM> <%s>",tgtPKStr);
    
      /* write number frames */
      if (strchr(oflags,'n')){
         fprintf(oFile,"\n<NFRAMES> %d",p->sa->NumFrame);
      }
      /* write mean */
      if (strchr(oflags,'m')){
         fprintf(oFile,"\n<MEAN> %d\n",vSize);
         for (i=1;i<=vSize;i++){
            fprintf(oFile," %e",(p->sa->meanSum[i]));
         }
      }
      /* write variance */
      if (strchr(oflags,'v')){   
         fprintf(oFile,"\n<VARIANCE> %d\n",vSize);
         for (i=1;i<=vSize;i++){
            fprintf(oFile," %e",(p->sa->squareSum[i]));
         }
      }
      fprintf(oFile,"\n");
      FClose(oFile,isPipe);
      p = p->nextSpkr;
   }
}
void AssignGMM(HLink hmm, Group* grps, int M)
{
   int n,m,s,S,i,j;
   StateElem *ste;
   MixtureElem *me;
   StreamElem *se;
   MixPDF *mp;
   /*if(forceFc)
	   CalcCovs();*/
   ste = hmm->svec+2; S = hmm->owner->swidth[0];
   for (n=2; n<hmm->numStates; n++,ste++) {
      se = ste->info->pdf+1;
      for (s=1; s<=S; s++,se++){
    	 se->nMix=M;
    	 se->spdf.cpdf=CreateCME(&hset,M);
         me = se->spdf.cpdf+1;
         me->rwgt=NULL;
         for (m=1; m<=se->nMix; m++,me++){
        	 me->weight=me->sweight=grps[m-1]->weight;
        	 mp=me->mpdf= (MixPDF *)New(hset.hmem,sizeof(MixPDF));
        	 mp->nUse = 0; mp->hook = NULL;
        	 mp->mIdx = 0; mp->stream = 0; mp->vFloor = NULL; mp->info = NULL;
        	 mp->mean=CreateSVector(hset.hmem,hset.vecSize);
        	 for(i=1;i<=hset.vecSize;i++){
        		 mp->mean[i]=grps[m-1]->mean[i];
        	 }
        	 /*CopyVector(grps[m-1]->mean,mp->mean);*/
        	 /*only support tied full covariance case*/
        	 if(forceFc){
                 mp->ckind=FULLC;
                 mp->cov.inv=CreateSTriMat(hset.hmem,hset.vecSize);
                 /*CopyMatrix(accs[s].fixed.inv,mp->cov.inv);*/
                 for(i=1;i<=hset.vecSize;i++){
                	 for(j=1;j<=i;j++){
                		 mp->cov.inv[i][j]=grps[m-1]->inv[i][j];
                	 }
                 }
                 /*CopyMatrix(grps[m-1]->inv,mp->cov.inv);*/
                 FixFullGConst(mp,-CovDet(&gstack, mp->cov.inv));
        	 }
        	 else{
            	 mp->ckind=DIAGC;
            	 mp->cov.var=CreateSVector(hset.hmem,hset.vecSize);
            	 for(i=1;i<=hset.vecSize;i++){
            		 mp->cov.var[i]=grps[m-1]->var[i];
            	 }
            	 /*CopyVector(grps[m-1]->var,mp->cov.var);*/
            	 FixDiagGConst(mp);
        	 }
         }
      }
   }
}

void AddPCA2InputXForm(){
	int i,j,vSize=hset.swidth[1];
	char newFn[MAXSTRLEN];
	Matrix fV=CreateMatrix(&gstack,vSize,vSize);
	Vector fd=CreateVector(&gstack,vSize);
	SMatrix pcaXForm=CreateSMatrix(hset.hmem,pca,vSize);
	Vector mu=CreateVector(&gstack,vSize);
	TriMat cov=CreateTriMat(&gstack,vSize);
	SVector pcaBias=CreateSVector(&gstack,vSize);
	float varxy,n;

	if (totalCount<2)
		HError(2021,"CalcCovs: Only %d speech frames accumulated",totalCount);
	if	(trace&T_TOP)
		printf("%ld speech frames accumulated\n", totalCount);
	n = (float)totalCount;     /* to prevent rounding to integer below */
	for (i=1; i<=vSize; i++)            /* For each coefficient ... */
		mu[i]=accs[1].meanSum[i] / n;         /* ... calculate mean */
	for (i=1;i<=vSize;i++) {
		for (j=1; j<=i; j++) {
			varxy = accs[1].squareSum.inv[i][j]/n - mu[i]*mu[j];
			if(!finite(varxy))
				HError(7191, "Infinite varxy!");
			if (isnan(varxy))
				HError(1, "varxy isnan..");
			cov[i][j] = varxy;
		}
	}
	EigenDecomposition(cov, fV, fd);
	for(i=1;i<=pca;i++){
		printf("Eigenvector index=%d, eigenvalue=%e\n",vSize-i+1,fd[vSize-i+1]);
		for(j=1;j<=vSize;j++){
			pcaXForm[i][j]=fV[j][vSize-i+1]/sqrt(fd[vSize-i+1]);/*iCov*PCA*(x-mu), Cov is the covariance of transformed features*/
		}
		pcaBias[i]=0;
		for(j=1;j<=vSize;j++){/*bias=-iCov*PCA*mu*/
			pcaBias[i]-=pcaXForm[i][j]*mu[j];
		}
	}
	MakeFN("pca",outDir,NULL,newFn);
	hset.xf = (InputXForm *)New(hset.hmem,sizeof(InputXForm));
	hset.xf->fname = CopyString(hset.hmem,newFn);
	hset.xf->xformName = "pca";
	hset.xf->mmfIdMask = "*";
	hset.xf->pkind = hset.pkind;
	hset.xf->preQual = FALSE;
	hset.xf->xform = (LinXForm *)New(hset.hmem,sizeof(LinXForm));
	hset.xf->xform->bias=pcaBias;
	hset.xf->xform->xform= (SMatrix*)New(hset.hmem, sizeof(SMatrix)*2);
	hset.xf->xform->xform[1]=pcaXForm;
	hset.xf->xform->blockSize = CreateIntVec(hset.hmem,1);
	hset.xf->xform->blockSize[1]=1;
	hset.xf->xform->vFloor=NULL;
	hset.xf->xform->nUse=0;
	hset.xf->xform->vecSize=vSize;
	hset.xf->nUse=0;
	SaveInputXForm(&hset,hset.xf,newFn,FALSE);
}


/* main func */
int main(int argc, char *argv[])
{
   char *datafn, *s;
   void Initialise(void);
   void LoadFile(char *fn, Point plist, int *pntCount);
   void SetCovs(void);
   void PutVFloor(void);
   void SaveModel(char *outfn);
   SpkrAcc *sa = NULL;
   Point plist=NULL;
   Point* points;
   Point iter;
   int pntCount=0;
   int gmm=1;
   int i;
   int numFiles=0,numThreads=0;

   if(InitShell(argc,argv,hcompv_version,hcompv_vc_id)<SUCCESS)
      HError(2000,"HCompV: InitShell failed");
   InitMem();   InitLabel();
   InitMath();  InitSigP();
   InitWave();  InitAudio();
   InitVQ();    InitModel();
   if(InitParm()<SUCCESS)  
      HError(2000,"HCompV: InitParm failed");

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(0);
   SetConfParms();

   CreateHMMSet(&hset,&gstack,FALSE);
   CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   CreateHMMSet(&synrhset,&hmmStack,TRUE);
   pathPattern[0]='\0';

   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      /*if (strlen(s)!=1)
         HError(2019,"HCompV: Bad switch %s; must be single letter",s);*/
      switch(s[0]){
      case 'f':
    	 if (!strcmp(s, "fc")){
    		forceFc=TRUE;
    		break;
    	 }
    	 if (!strcmp(s, "fl")){
    		forceLog=TRUE;
    		break;
    	 }
         if (NextArg() != FLOATARG)
            HError(2019,"HCompV: Variance floor scale expected");
         vFloorScale = GetChkedFlt(0.0,100.0,s);
         break;
      case 'n':
    	  if (!strcmp(s,"nf")){
    		  if(NextArg()!=INTARG){
    			  HError(2019,"HCompV: Number of files expected");
    		  }
    		  numFiles=GetIntArg();
    		  fileList=(char**)New(&gstack,sizeof(char*)*numFiles);
    		  fileList2=(char**)New(&gstack,sizeof(char*)*numFiles);
    	  }
    	  if (!strcmp(s,"nt")){
    		  if(NextArg()!=INTARG){
    			  HError(2019,"HCompV: Number of threads expected");
    		  }
    		  numThreads=GetIntArg();
    	  }
    	  break;
      case 'l':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: Segment label expected");
         segLab = GetStrArg();
         break;
      case 'm':
         meanUpdate = TRUE;
         break;
      case 't':
    	 if (!strcmp(s, "trim")){
    		 DoTrim=TRUE;
        	 if (NextArg()!=INTARG)
        		 HError(2319,"HCompV: to calculate noise model using head&tail frames, expecting trimming frame number (e.g.20)");
        	 trimNum = GetIntArg();
    	 }
    	 else if(!strcmp(s, "tiedCov")){
    		 tiedCov=TRUE;
    	 }
    	 else{
    		 HError(2019,"HCompV: Unknown switch %s",s);
    	 }
    	 break;
      case 'g':
    	 if (!strcmp(s,"gmm")){
    		 plist=CreatePoint(&gstack,NULL,NULL,NULL);
    		 if(NextArg()!=INTARG)
    			 HError(2319,"HCompV: to calculate a GMM noise model, expecting number of mixtures");
    		 gmm=GetIntArg();
    	 }
    	 else{
    		 HError(2019,"HCompV: Unknown switch %s",s);
    	 }
    	 break;
      case 'o':
    	  if (!strcmp(s,"os"))
    		  onlyStatic=TRUE;
    	  else if(strlen(s)==1)
    		  outfn = GetStrArg();
    	  break;
      case 'v':
         if (NextArg() != FLOATARG)
            HError(2019,"HCompV: Minimum variance level expected");
         minVar = GetChkedFlt(0.0,100.0,s);
         break;
      case 'k':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: speaker pattern expected");
         strcpy(spPattern,GetStrArg());
         if (strchr(spPattern,'%')==NULL)
            HError(2019,"HCompV: Speaker mask invalid");
         break;
      case 'c':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: CMV output dir expected");
         strcpy(cmDir,GetStrArg());
         DoCMV = TRUE;
         break;
      case 'p':
    	 if (!strcmp(s,"pca")){
    		 if(NextArg()!=INTARG){
    			 HError(1,"HCompV: target dimension of PCA transformation expected");
    		 }
    		 pca=GetIntArg();
    		 break;
    	 }
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: path pattern expected");
         strcpy(pathPattern,GetStrArg());
         if (strchr(pathPattern,'%')==NULL)
            HError(2019,"HCompV: Path mask invalid");
         break;
      case 'q':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: output flags (nmv)");
         strcpy(oflags,GetStrArg());
         break;
      case 'B':
         saveBinary = TRUE;
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: Data File format expected");
         if((dff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2089,"HCompV: Warning ALIEN Data file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: Label File format expected");
         if((lff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2089,"HCompV: Warning ALIEN Label file format set");
         break;
      case 'H':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: HMM macro file name expected");
         AddMMF(&hset,GetStrArg());
         break;
      case 'I':
         if (NextArg() != STRINGARG)
            HError(2019,"HCompV: MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'L':
         if (NextArg()!=STRINGARG)
            HError(2019,"HCompV: Label file directory expected");
         labDir = GetStrArg();
         break;
      case 'M':
         if (NextArg()!=STRINGARG)
            HError(2019,"HCompV: Output macro file directory expected");
         outDir = GetStrArg();
         break;
      case 'T':
         if (NextArg() != INTARG)
            HError(2019,"HCompV: Trace value expected");
         trace = GetChkedInt(0,077,s); 
         break;
      case 'X':
         if (NextArg()!=STRINGARG)
            HError(2019,"HCompV: Label file extension expected");
         labExt = GetStrArg();
         break;
      case 'E':
    	  if (NextArg()!=INTARG)
    		HError(1019,"HCompV: Feature expansion window size expected");
    	  xpwin=GetIntArg();break;
      case 'W':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCompV: -W: gmmcluster definition file expected");
    	 AddMMF(&synrhset,GetStrArg());
    	 gotgmmcluster=TRUE;
    	 break;
      case 'e':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCompV: -p: gmmcluster class list file expected");
    	 strcpy(gmmclusterListFn,GetStrArg());
    	 break;
      case 's':
    	  stereo=TRUE;break;
      default:
         HError(2019,"HCompV: Unknown switch %s",s);
         break;
      }
   }
   if(gotgmmcluster){/*load gmmcluster */
	   gmmcluster=LoadPosSynth(&synrhset,gmmclusterListFn);
   }
   /* if not doing CMV, do standard HCompV */
   if (DoCMV == FALSE){
      if (NextArg()!=STRINGARG)
         HError(2019,"HCompV: Source HMM file name expected");
      hmmfn = GetStrArg();
      Initialise();
      i=0;
      char ext[12];
      do {
         if (NextArg()!=STRINGARG)
            HError(2019,"HCompV: Training data file name expected");
         datafn = GetStrArg();
         if(numThreads==0){
        	 ExtnOf(datafn,ext);
        	 printf("%s\n",datafn);fflush(stdout);
        	 if(!strcmp(ext,"spa")){
        		 LoadSpaf(datafn,hset.vecSize);
        	 }
        	 else{
        		 LoadFile(datafn, plist, &pntCount);
        	 }
         }
         else{
        	 fileList[i]=(char*)New(&gstack,sizeof(char)*(strlen(datafn)+1));/*last char for terminal null character*/
        	 strcpy(fileList[i],datafn);
        	 if(stereo){
        		 datafn = GetStrArg();
            	 fileList2[i]=(char*)New(&gstack,sizeof(char)*(strlen(datafn)+1));/*last char for terminal null character*/
            	 strcpy(fileList2[i],datafn);
        	 }
        	 i++;
        	 if(i==numFiles){
        		 break;
        	 }
         }
      } while (NumArgs()>0);
      if(numThreads>0){
    	  if(pca){
    		  AccParallel(fileList, numFiles, hset.vecSize, numThreads);
    		  AddPCA2InputXForm();
    		  SetCovs();
    	  }
    	  else{
    		  Group *grps;
    		  if(gotgmmcluster){
    			  grps=GMMClusterFileX(fileList, numFiles, gmm, hset.vecSize, numThreads, forceFc, tiedCov);
    		  }
    		  else{
    			  grps=GMMClusterFile(fileList, numFiles, gmm, hset.vecSize, numThreads, forceFc, tiedCov);
    		  }
    		  AssignGMM(hmmLink,grps,gmm);
    	  }
    	  SaveModel(outfn);
      }
      else{
          if(gmm>1){
        	  points=(Point*)New(&gstack,sizeof(Point)*pntCount);
        	  iter=plist->next;
        	  for(i=0;i<pntCount;i++,iter=iter->next){
        		  points[i]=iter;
        	  }
        	  Group *grps=InitKmeansMU(points, gmm, hset.vecSize, pntCount, forceFc, tiedCov);
        	  AssignGMM(hmmLink,grps,gmm);
        	  SaveModel(outfn);
          }
          else{
              if(pca){/*save PCA transformation*/
            	  AddPCA2InputXForm();
              }
              SetCovs();
              if(!pca){
            	  FixGConsts(hmmLink);
              }
              SaveModel(outfn);
          }
      }
      if (trace&T_TOP)
         printf("Output written to directory %s\n",(outDir==NULL)?"./":outDir);
      if (vFloorScale>0.0)
         PutVFloor();
   }
   else {
      /* report export data type */
      ReportOutput();
      /* init input buffer mem heap */
      CreateHeap(&iStack,"BufferIn",MSTAK,1,0.5,100000,5000000);
      do {
         if (NextArg()!=STRINGARG){
            HError(2019,"HCompV: Training data file name expected");
         }
         datafn = GetStrArg();       
         /* accumulate stats for current utterance file and update speaker list */
         sa = AccGenUtt(spPattern,datafn,sa);
         salist = UpdateSpkrAccList(salist,sa);
         /* reset for next utterance */
         ClrSpkrAcc(sa);
      } while (NumArgs()>0);
      /* compute the means and variances for each speaker */
      UpdateMeanVar(salist);
      /* export NMV for each speaker */
      ExportNMV(salist,cmDir,TargetPKStr);
   }

   Exit(0);
   return (0);          /* never reached -- make compiler happy */
}

typedef struct _WorkLoad{
	int id;
	int vSize;
	int numFiles;
	int grpCount;/*current estimated number of clusters*/
	int numClust;/*total clusters*/
	float measure;
	Group* grps;/*a clone of current grps*/
	char** fileList;/*subset of filelist*/
	char** fileList2;/*subset of filelist*/
	MemHeap iStack;/*local memheap*/
	MemHeap iStack2;/*local memheap*/
}*WorkLoad;


void CopyGroups(Group* from, Group* to, int grpCount){
	int i;
	for(i=0;i<grpCount;i++){
		to[i]->fc=from[i]->fc;
		to[i]->gConst=from[i]->gConst;
		to[i]->weight=from[i]->weight;
		CopyVector(from[i]->mean,to[i]->mean);
		if(to[i]->fc){
			CopyTriMat(from[i]->inv,to[i]->inv);
		}
		else{
			CopyVector(from[i]->var,to[i]->var);
		}
	}
}
WorkLoad CreateWorkLoad(int id, int numFiles, int numClust, char** fileList, int vSize, Boolean fc, char** fileList2){
	int i;
	WorkLoad wl=New(&gstack,sizeof(struct _WorkLoad));
	wl->grps=(Group*)New(&gstack,sizeof(Group)*numClust);
	for(i=0;i<numClust;i++){
		wl->grps[i]=CreateGroup(&gstack,vSize,fc);
	}
	wl->id=id;
	wl->numFiles=numFiles;
	wl->vSize=vSize;
	wl->numClust=numClust;
	wl->fileList=fileList;
	wl->fileList2=fileList2;
	wl->measure=0;
	CreateHeap(&wl->iStack,"InBuf", MSTAK, 1, 0.5, 100000, LONG_MAX);
	if(stereo)
		CreateHeap(&wl->iStack2,"InBuf", MSTAK, 1, 0.5, 100000, LONG_MAX);
	return wl;
}
void* WorkThread(void* arg){
	WorkLoad wl=(WorkLoad)arg;
	ParmBuf pbuf, pbuf2;
	BufferInfo info, info2;
	int i,nObs,fidx,j,k,idx;
	char* fn,*fn2;
	int obsvsize,lobsvsize;
	int vecSize=wl->vSize;
	Vector xpObs,tempObs;
	int halfwin;
	Boolean offsetObs;
	int swidth0=1;
	wl->measure=0;
	Boolean eSep,eSep2;
	short swidth[5],swidth2[5];
	swidth[0]=1;/*only support one stream*/
	swidth2[0]=1;/*only support one stream*/
	Observation obs,obs2;
	for(fidx=0;fidx<wl->numFiles;fidx++){
		fn=wl->fileList[fidx];
		/*printf("Processing %s by thread %d\n",fn,wl->id);*/
		if((pbuf=OpenBuffer(&wl->iStack, fn, 0, dff, FALSE_dup, FALSE_dup))==NULL)
				HError(2050,"LoadFile: Config parameters invalid");
		GetBufferInfo(pbuf,&info);
		if(stereo){
			fn2=wl->fileList2[fidx];
			if((pbuf2=OpenBuffer(&wl->iStack2, fn2, 0, dff, FALSE_dup, FALSE_dup))==NULL)
					HError(2050,"LoadFile: Config parameters invalid");
			GetBufferInfo(pbuf2,&info2);
		}
		swidth[1]=info.tgtVecSize;
		SetStreamWidths(info.tgtPK,info.tgtVecSize,swidth,&eSep);
		obs=MakeObservation(&wl->iStack,swidth,info.tgtPK,FALSE,eSep);
		if(stereo){
			swidth2[1]=info2.tgtVecSize;
			SetStreamWidths(info2.tgtPK,info2.tgtVecSize,swidth2,&eSep2);
			obs2=MakeObservation(&wl->iStack2,swidth2,info2.tgtPK,FALSE,eSep2);
		}
		/*CheckData(fn,info);*/
		if(onlyStatic){
			obsvsize=NumStatic(info.tgtVecSize,info.tgtPK);
			if (info.tgtPK&HASENERGY) ++obsvsize;
			if (info.tgtPK&HASZEROC) ++obsvsize;
			lobsvsize=obsvsize*xpwin;/*cover only static part*/
		}
		else{
			obsvsize=info.tgtVecSize;/*cover all*/
			lobsvsize=obsvsize*xpwin;
		}
		xpObs=CreateVector(&wl->iStack,lobsvsize*(stereo?2:1));
		halfwin=xpwin/2;
		offsetObs=(xpwin-1)%2;/*e.g. win=9->no offset, win=8->offset*/
		tempObs=CreateVector(&wl->iStack,info.tgtVecSize);
		nObs = ObsInBuffer(pbuf);
		for(i=0; i < nObs; i++){
			if(offsetObs){
				ReadAsTable(pbuf,i,&obs);
				if(stereo)
					ReadAsTable(pbuf2,i,&obs2);
				CopyVector(obs.fv[swidth0],tempObs);
			}
			for(j=-halfwin;j<=halfwin;j++){
				if(offsetObs&&j==0){/*no need to read again the current frame*/
					continue;
				}
				k=i+j;
				k=k<0?0:(k>=nObs?nObs-1:k);
				ReadAsTable(pbuf,k,&obs);
				for(idx=1;idx<=obsvsize;idx++){/*tgtvsize refers to the part of features we care about*/
					if(offsetObs){
						xpObs[idx+obsvsize*((j>0?j-1:j)+halfwin)]=obs.fv[swidth0][idx];/*-tempObs[idx];/*offset the current frame*/
					}
					else{
						xpObs[idx+obsvsize*(j+halfwin)]=obs.fv[swidth0][idx];
					}
				}
				if(stereo){
					ReadAsTable(pbuf2,k,&obs2);
					for(idx=1;idx<=obsvsize;idx++){/*tgtvsize refers to the part of features we care about*/
						if(offsetObs){
							xpObs[lobsvsize+idx+obsvsize*((j>0?j-1:j)+halfwin)]=obs2.fv[swidth0][idx];/*-tempObs[idx];/*offset the current frame*/
						}
						else{
							xpObs[lobsvsize+idx+obsvsize*(j+halfwin)]=obs2.fv[swidth0][idx];
						}
					}
				}
			}
			Vector vec=CreateVector(&wl->iStack,vecSize);
			CopyVector(xpObs,vec);
			Point pnt=CreatePoint(&wl->iStack,vec,NULL,NULL);
			wl->measure+=SearchAndInsert(&wl->iStack,wl->grps, pnt, wl->grpCount, Gaussi, FALSE, FALSE);
		}
        CloseBuffer(pbuf);
        ResetHeap(&wl->iStack);
        if(stereo){
        	CloseBuffer(pbuf2);
        	ResetHeap(&wl->iStack2);
        }
	}
	return NULL;
}



void* AccThread(void* arg){
	AccLoad al=(AccLoad)arg;
	int vSize=VectorSize(al->meanSum);
	int idx,i,j,dnn_row,dnn_col;
	char ch;
	char dnn_buf[32];
	float value;
	int fidx;
	FILE* dnn_f;
	int x,y;
	float valx,valy;
	Vector v=al->tempV;
	for(fidx=0;fidx<al->numFiles;fidx++){
		if((dnn_f=fopen(al->fileList[fidx],"r"))==NULL){
			HError(1,"Cannot open file: %s",al->fileList[fidx]);
		}
		else{
			printf("Process %s\n",al->fileList[fidx]);
		}
		fscanf(dnn_f,"%d %d\n",&dnn_row,&dnn_col);
		assert(dnn_col==vSize);
		for(i=1;i<=dnn_row;i++){
			for(idx=1;idx<=vSize;idx++){/*set all small posterior to be floor rather than zero to avoid LZERO*/
				v[idx]=SMALLP;
			}
			j=0;
			while((ch=getc(dnn_f))!=EOF){
				if(ch==' '||ch=='\n'){
					dnn_buf[j++]='\0';
					if(j>1){
						sscanf(dnn_buf,"%d:%e",&idx,&value);
						v[idx]=value;
					}
					j=0;
					if(ch=='\n')
						break;
				}
				else{
					dnn_buf[j++]=ch;
				}
			}
			al->numFrames++;
			for (x=1;x<=vSize;x++) {
				valx=forceLog?(v[x]>0?log(v[x]):LSMALLP):v[x];
				al->meanSum[x] += valx;     /* accumulate mean */
				for (y=1;y<=x;y++){
					valy=forceLog?(v[y]>0?log(v[y]):LSMALLP):v[y];
					al->squareSum[x][y] += valx*valy;
				}
			}
		}
		fclose(dnn_f);
	}
	return NULL;
}

void AccParallel(char** fileList, int numFiles, int vSize, int numThreads){
	int i,j;
	AccLoad *al;
	int nt;
	al=(AccLoad*)New(&gstack,sizeof(AccLoad)*numThreads);
	int seg=numFiles/numThreads+1;
	int comsumedFiles=0;
	int procSeg=0;

	pthread_t *pths=(pthread_t*)New(&gstack,sizeof(pthread_t)*numThreads);
	for(nt=0;nt<numThreads;nt++){
		procSeg=seg;
		if(comsumedFiles+procSeg>numFiles)
			procSeg=numFiles-comsumedFiles;
		al[nt]=(AccLoad)New(&gstack,sizeof(struct _AccLoad));
		al[nt]->fileList=fileList+comsumedFiles;
		al[nt]->meanSum=CreateVector(&gstack,vSize);
		al[nt]->squareSum=CreateTriMat(&gstack,vSize);
		al[nt]->numFiles=procSeg;
		al[nt]->numFrames=0;
		al[nt]->tempV=CreateVector(&gstack,vSize);
		ZeroVector(al[nt]->meanSum);
		ZeroTriMat(al[nt]->squareSum);
		comsumedFiles+=procSeg;
	}
	for(nt=0;nt<numThreads;nt++){
		pthread_create(&pths[nt],NULL,AccThread,al[nt]);
	}
	for(nt=0;nt<numThreads;nt++){
		pthread_join(pths[nt], NULL /* void ** return value could go here */);
	}
	totalCount=0;
	for(nt=0;nt<numThreads;nt++){
		totalCount+=al[nt]->numFrames;
		for(i=1;i<=vSize;i++){
			accs[1].meanSum[i]+=al[nt]->meanSum[i];
			for(j=1;j<=i;j++){
				accs[1].squareSum.inv[i][j]+=al[nt]->squareSum[i][j];
			}
		}
	}
}

Group* GMMClusterFile(char** fileList, int numFiles, int numClust, int vSize, int numThreads, Boolean fc, Boolean tiedCov){
	int i,j,grpCount=1,k,iter;
	double x;
	double measure;
	const float pertDepth = 0.2;
	Group *grps,*pGrps;
	WorkLoad *wl;
	TriMat fcvar;
	int nt;
	float numPnts;
	float normWgt;
	float* floorV;
	pGrps=(Group*)New(&gstack,sizeof(Group)*numClust);
	grps=(Group*)New(&gstack,sizeof(Group)*numClust);
	wl=(WorkLoad*)New(&gstack,sizeof(WorkLoad)*numThreads);
	int seg=numFiles/numThreads+1;
	int comsumedFiles=0;
	int procSeg=0;
	int toSplitCount;
	if(stereo)
		vSize*=2;
	pthread_t *pths=(pthread_t*)New(&gstack,sizeof(pthread_t)*numThreads);
	for(nt=0;nt<numThreads;nt++){
		procSeg=seg;
		if(comsumedFiles+procSeg>numFiles)
			procSeg=numFiles-comsumedFiles;
		wl[nt]=CreateWorkLoad(nt,procSeg, numClust, fileList+comsumedFiles, vSize, fc, fileList2+comsumedFiles);
		comsumedFiles+=procSeg;
	}
	fcvar=CreateTriMat(&gstack,vSize);
	for(i=0;i<numClust;i++){
		grps[i]=CreateGroup(&gstack,vSize,fc);
		pGrps[i]=CreateGroup(&gstack,vSize,fc);
	}
	pGrps[0]->weight=1;
	floorV=pGrps[0]->floorV;
	for(nt=0;nt<numThreads;nt++){
		wl[nt]->grpCount=grpCount;
		CopyGroups(pGrps, wl[nt]->grps, grpCount);
		ResetGroups(wl[nt]->grps,grpCount);
		pthread_create(&pths[nt],NULL,WorkThread,wl[nt]);
	}
	for(nt=0;nt<numThreads;nt++){
		pthread_join(pths[nt], NULL /* void ** return value could go here */);
	}
	measure=0;
	ResetGroups(pGrps,grpCount);
	for(nt=0;nt<numThreads;nt++){
		measure+=wl[nt]->measure;
		for(i=0;i<grpCount;i++){
			pGrps[i]->occ+=wl[nt]->grps[i]->occ;
			for(j=1;j<=vSize;j++){
				pGrps[i]->acc1[j]+=wl[nt]->grps[i]->acc1[j];
			}
			if(fc){
				for(j=1;j<=vSize;j++){
					for(k=1;k<=j;k++){
						pGrps[i]->acc2x[j][k]+=wl[nt]->grps[i]->acc2x[j][k];
					}
				}
			}
			else{
				for(j=1;j<=vSize;j++)
					pGrps[i]->acc2[j]+=wl[nt]->grps[i]->acc2[j];
			}
		}
	}
	numPnts=pGrps[0]->occ;
	UpdateGroups(&gstack, pGrps, grpCount, numPnts, TRUE);
	if(fc){
		CovInvert(&gstack, pGrps[0]->inv,fcvar);/*Note that you can't just invert the diagonal element to get the variance*/
	}
	for(j=1;j<=vSize;j++){
		floorV[j]=fc?(0.01*fcvar[j][j]):(0.01*pGrps[0]->var[j]);
	}
	printf("Starting from measure %f with one mixture [mu=0,var=1]\n",measure/numPnts);
	while(grpCount<numClust){
		toSplitCount=numClust-grpCount;
		toSplitCount=((toSplitCount<grpCount)?toSplitCount:grpCount);
		if(toSplitCount)
			normWgt=1.0/(grpCount+toSplitCount);
		else
			normWgt=0;
		for(i=0;i<toSplitCount;i++){/*mix up*/
			grps[i*2]->weight=normWgt;/*pGrps[i]->weight/2.0;/*split the weight*/
			grps[i*2+1]->weight=normWgt;/*pGrps[i]->weight/2.0;*/
			grps[i*2]->gConst=pGrps[i]->gConst;/*keep the variance*/
			grps[i*2+1]->gConst=pGrps[i]->gConst;
			if(fc){
				CovInvert(&gstack, pGrps[i]->inv,fcvar);/*Note that you can't just invert the diagonal element to get the variance*/
			}
			for(j=1;j<=vSize;j++){/*Although the split algorithm is the same, the reestimated mean could be very different using diagonal or full covariance matrix*/
				x = sqrt(fc?fcvar[j][j]:pGrps[i]->var[j])*pertDepth;
				grps[i*2]->mean[j]=pGrps[i]->mean[j]+x;/*pertub the mean*/
				grps[i*2+1]->mean[j]=pGrps[i]->mean[j]-x;
			}
			if(fc){
				CopyMatrix(pGrps[i]->inv,grps[i*2]->inv);
				CopyMatrix(pGrps[i]->inv,grps[i*2+1]->inv);
			}
			else{
				CopyVector(pGrps[i]->var,grps[i*2]->var);
				CopyVector(pGrps[i]->var,grps[i*2+1]->var);
			}
		}
        for(i=toSplitCount;i<grpCount;i++){/*Do hard copy instead of pointer copy because CopyGroups(grps, pGrps, grpCount) call will mess it up*/
        	grps[2*toSplitCount+i-toSplitCount]->weight=normWgt?normWgt:pGrps[i]->weight;
        	grps[2*toSplitCount+i-toSplitCount]->gConst=pGrps[i]->gConst;
        	CopyVector(pGrps[i]->mean,grps[2*toSplitCount+i-toSplitCount]->mean);
        	if(fc){
        		CopyMatrix(pGrps[i]->inv,grps[2*toSplitCount+i-toSplitCount]->inv);
        		CopyMatrix(pGrps[i]->inv,grps[2*toSplitCount+i-toSplitCount]->inv);
        	}
        	else{
        		CopyVector(pGrps[i]->var,grps[2*toSplitCount+i-toSplitCount]->var);
        		CopyVector(pGrps[i]->var,grps[2*toSplitCount+i-toSplitCount]->var);
        	}
        }
		grpCount+=toSplitCount;
		printf("\tIncreased to %d mixes\n",grpCount);
		/*printVector(grps[0]->mean,"m1");
		printVector(grps[1]->mean,"m1");*/
		for(iter=0;iter<4;iter++){/*do 4 iteration EM estimation*/
			for(nt=0;nt<numThreads;nt++){
				wl[nt]->grpCount=grpCount;
				CopyGroups(grps, wl[nt]->grps, grpCount);
				ResetGroups(wl[nt]->grps,grpCount);
				pthread_create(&pths[nt],NULL,WorkThread,wl[nt]);
			}
			for(nt=0;nt<numThreads;nt++){
				pthread_join(pths[nt], NULL /* void ** return value could go here */);
			}
			measure=0;
			ResetGroups(grps,grpCount);/*reset the accumulators*/
			for(nt=0;nt<numThreads;nt++){
				measure+=wl[nt]->measure;
				for(i=0;i<grpCount;i++){
					grps[i]->occ+=wl[nt]->grps[i]->occ;
					for(j=1;j<=vSize;j++){
						grps[i]->acc1[j]+=wl[nt]->grps[i]->acc1[j];
					}
					if(fc){
						for(j=1;j<=vSize;j++){
							for(k=1;k<=j;k++){
								grps[i]->acc2x[j][k]+=wl[nt]->grps[i]->acc2x[j][k];
							}
						}
					}
					else{
						for(j=1;j<=vSize;j++){
							grps[i]->acc2[j]+=wl[nt]->grps[i]->acc2[j];
						}
					}
				}
			}
			UpdateGroups(&gstack, grps, grpCount, numPnts, tiedCov?FALSE:TRUE);/*update mean&variance&weight of the cluster*/
			printf("\t\tIteration %d with measure %f\n",iter,measure/numPnts);
			/*printVector(grps[0]->mean,"m1");
			printVector(grps[1]->mean,"m1");*/
		}
		CopyGroups(grps, pGrps, grpCount);
	}
	return grps;
}

Group* GMMClusterFileX(char** fileList, int numFiles, int numClust, int vSize, int numThreads, Boolean fc, Boolean tiedCov){
	int i,j,grpCount=1,k,iter;
	double x;
	double measure;
	const float pertDepth = 0.2;
	Group *grps,*pGrps;
	WorkLoad *wl;
	TriMat fcvar;
	int nt;
	float numPnts;
	float* floorV;
	grps=(Group*)New(&gstack,sizeof(Group)*numClust);
	pGrps=(Group*)New(&gstack,sizeof(Group)*numClust);
	wl=(WorkLoad*)New(&gstack,sizeof(WorkLoad)*numThreads);
	int seg=numFiles/numThreads+1;
	int comsumedFiles=0;
	int procSeg=0;
	int toSplitCount;
	float normWgt;
	pthread_t *pths=(pthread_t*)New(&gstack,sizeof(pthread_t)*numThreads);
	for(nt=0;nt<numThreads;nt++){
		procSeg=seg;
		if(comsumedFiles+procSeg>numFiles)
			procSeg=numFiles-comsumedFiles;
		wl[nt]=CreateWorkLoad(nt,procSeg, numClust, fileList+comsumedFiles, vSize, fc, NULL);
		comsumedFiles+=procSeg;
	}
	fcvar=CreateTriMat(&gstack,vSize);
	for(i=0;i<numClust;i++){
		grps[i]=CreateGroup(&gstack,vSize,fc);
		pGrps[i]=CreateGroup(&gstack,vSize,fc);
	}
	grpCount=synrhset.numPhyHMM;
	for(i=0;i<grpCount;i++){
		MixPDF* mp=(((gmmcluster[i+1]->svec+2)->info->pdf+1)->spdf.cpdf+1)->mpdf;
		pGrps[i]->weight=1.0/grpCount;
		pGrps[i]->gConst=-0.5*mp->gConst;
		CopyVector(mp->mean,pGrps[i]->mean);
		if(fc){
			CopyMatrix(mp->cov.inv,pGrps[i]->inv);
		}
		else{
			CopyVector(mp->cov.var,pGrps[i]->var);
		}
	}
	floorV=pGrps[0]->floorV;
	for(nt=0;nt<numThreads;nt++){
		wl[nt]->grpCount=grpCount;
		CopyGroups(pGrps, wl[nt]->grps, grpCount);
		ResetGroups(wl[nt]->grps,grpCount);
		pthread_create(&pths[nt],NULL,WorkThread,wl[nt]);
	}
	for(nt=0;nt<numThreads;nt++){
		pthread_join(pths[nt], NULL /* void ** return value could go here */);
	}
	measure=0;
	ResetGroups(pGrps,grpCount);
	for(nt=0;nt<numThreads;nt++){
		measure+=wl[nt]->measure;
		for(i=0;i<grpCount;i++){
			pGrps[i]->occ+=wl[nt]->grps[i]->occ;
			for(j=1;j<=vSize;j++){
				pGrps[i]->acc1[j]+=wl[nt]->grps[i]->acc1[j];
			}
			if(fc){
				for(j=1;j<=vSize;j++){
					for(k=1;k<=j;k++){
						pGrps[i]->acc2x[j][k]+=wl[nt]->grps[i]->acc2x[j][k];
					}
				}
			}
			else{
				for(j=1;j<=vSize;j++)
					pGrps[i]->acc2[j]+=wl[nt]->grps[i]->acc2[j];
			}
		}
	}
	numPnts=0;
	ZeroVector(floorV);
	for(i=0;i<grpCount;i++){
		numPnts+=pGrps[i]->occ;
		if(fc){
			CovInvert(&gstack, pGrps[i]->inv,fcvar);/*Note that you can't just invert the diagonal element to get the variance*/
		}
		for(j=1;j<=vSize;j++){
			floorV[j]+=(fc?(0.01*fcvar[j][j]):(0.01*pGrps[i]->var[j]))/grpCount;
		}
	}
	printf("Starting from measure %f with %d mixture \n",measure/numPnts,grpCount);
	UpdateGroups(&gstack, pGrps, grpCount, numPnts, TRUE);
	while(grpCount<=numClust){
		toSplitCount=numClust-grpCount;
		toSplitCount=((toSplitCount<grpCount)?toSplitCount:grpCount);
		if(toSplitCount)
			normWgt=1.0/(grpCount+toSplitCount);
		else
			normWgt=0;
		for(i=0;i<toSplitCount;i++){/*mix up*/
			grps[i*2]->weight=normWgt;/*pGrps[i]->weight/2.0;/*split the weight*/
			grps[i*2+1]->weight=normWgt;/*pGrps[i]->weight/2.0;*/
			grps[i*2]->gConst=pGrps[i]->gConst;/*keep the variance*/
			grps[i*2+1]->gConst=pGrps[i]->gConst;
			if(fc){
				CovInvert(&gstack, pGrps[i]->inv,fcvar);/*Note that you can't just invert the diagonal element to get the variance*/
			}
			for(j=1;j<=vSize;j++){/*Although the split algorithm is the same, the reestimated mean could be very different using diagonal or full covariance matrix*/
				x = sqrt(fc?fcvar[j][j]:pGrps[i]->var[j])*pertDepth;
				grps[i*2]->mean[j]=pGrps[i]->mean[j]+x;/*pertub the mean*/
				grps[i*2+1]->mean[j]=pGrps[i]->mean[j]-x;
			}
			if(fc){
				CopyMatrix(pGrps[i]->inv,grps[i*2]->inv);
				CopyMatrix(pGrps[i]->inv,grps[i*2+1]->inv);
			}
			else{
				CopyVector(pGrps[i]->var,grps[i*2]->var);
				CopyVector(pGrps[i]->var,grps[i*2+1]->var);
			}
		}
        for(i=toSplitCount;i<grpCount;i++){/*Do hard copy instead of pointer copy because CopyGroups(grps, pGrps, grpCount) call will mess it up*/
        	grps[2*toSplitCount+i-toSplitCount]->weight=normWgt?normWgt:pGrps[i]->weight;
        	grps[2*toSplitCount+i-toSplitCount]->gConst=pGrps[i]->gConst;
        	CopyVector(pGrps[i]->mean,grps[2*toSplitCount+i-toSplitCount]->mean);
        	if(fc){
        		CopyMatrix(pGrps[i]->inv,grps[2*toSplitCount+i-toSplitCount]->inv);
        		CopyMatrix(pGrps[i]->inv,grps[2*toSplitCount+i-toSplitCount]->inv);
        	}
        	else{
        		CopyVector(pGrps[i]->var,grps[2*toSplitCount+i-toSplitCount]->var);
        		CopyVector(pGrps[i]->var,grps[2*toSplitCount+i-toSplitCount]->var);
        	}
        }
		grpCount+=toSplitCount;
		printf("\tIncreased to %d mixes\n",grpCount);
		/*printVector(grps[0]->mean,"m1");
		printVector(grps[1]->mean,"m1");*/
		for(iter=0;iter<4;iter++){/*do 4 iteration EM estimation*/
			for(nt=0;nt<numThreads;nt++){
				wl[nt]->grpCount=grpCount;
				CopyGroups(grps, wl[nt]->grps, grpCount);
				ResetGroups(wl[nt]->grps,grpCount);
				pthread_create(&pths[nt],NULL,WorkThread,wl[nt]);
			}
			for(nt=0;nt<numThreads;nt++){
				pthread_join(pths[nt], NULL /* void ** return value could go here */);
			}
			measure=0;
			ResetGroups(grps,grpCount);/*reset the accumulators*/
			for(nt=0;nt<numThreads;nt++){
				measure+=wl[nt]->measure;
				for(i=0;i<grpCount;i++){
					grps[i]->occ+=wl[nt]->grps[i]->occ;
					for(j=1;j<=vSize;j++){
						grps[i]->acc1[j]+=wl[nt]->grps[i]->acc1[j];
					}
					if(fc){
						for(j=1;j<=vSize;j++){
							for(k=1;k<=j;k++){
								grps[i]->acc2x[j][k]+=wl[nt]->grps[i]->acc2x[j][k];
							}
						}
					}
					else{
						grps[i]->acc2[j]+=wl[nt]->grps[i]->acc2[j];
					}
				}
			}
			UpdateGroups(&gstack, grps, grpCount, numPnts, tiedCov?FALSE:TRUE);/*update mean&variance&weight of the cluster*/
			printf("\t\tIteration %d with measure %f\n",iter,measure/numPnts);
			/*printVector(grps[0]->mean,"m1");
			printVector(grps[1]->mean,"m1");*/
		}
		if(toSplitCount==0)break;
		CopyGroups(grps, pGrps, grpCount);
	}
	return grps;
}

/* ----------------------------------------------------------- */
/*                      END:  HCompV.c                         */
/* ----------------------------------------------------------- */
