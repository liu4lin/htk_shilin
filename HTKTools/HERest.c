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
/*          2002-2004 Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HERest.c: Embedded B-W ReEstimation           */
/* ----------------------------------------------------------- */

char *herest_version = "!HVER!HERest:   3.4.1 [CUED 12/03/09]";
char *herest_vc_id = "$Id: HERest.c,v 1.2 2006/12/07 11:09:08 mjfg Exp $";

/*
   This program is used to perform a single reestimation of
   the parameters of a set of HMMs using Baum-Welch.  Training
   data consists of one or more utterances each of which has a 
   transcription in the form of a standard label file (segment
   boundaries are ignored).  For each training utterance, a
   composite model is effectively synthesized by concatenating
   the phoneme models given by the transcription.  Each phone
   model has the usual set of accumulators allocated to it,
   these are updated by performing a standard B-W pass over
   each training utterance using the composite model. This program
   supports arbitrary parameter tying and multiple data streams.
   
   Added in V1.4 - support for tee-Models ie HMMs with a non-
   zero transition from entry to exit states.

   In v2.2 most of the core functionality has been moved to the
   library module HFB
*/

#include "HShell.h"     /* HMM ToolKit Modules */
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HAudio.h"
#include "HWave.h"
#include "HVQ.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HTrain.h"
#include "HUtil.h"
#include "HAdapt.h"
#include "HMap.h"
#include "HFB.h"
#include "pthread.h"

/* Trace Flags */
#define T_TOP   0001    /* Top level tracing */
#define T_MAP   0002    /* logical/physical hmm map */
#define T_UPD   0004    /* Model updates */
#define E	2.7182818284 /*natural logrigthm*/


/* possible values of updateMode */
#define UPMODE_DUMP 1
#define UPMODE_UPDATE 2
#define UPMODE_BOTH 3

#define PJMMODE_UNIFORM 11		/*uniform distribution*/
#define PJMMODE_GSMPL   12		/*use global samples*/
#define PJMMODE_MSMPL	13		/*use component dependent samples*/

/* Global Settings */

static char * labDir = NULL;     /* label (transcription) file directory */
static char * labExt = "lab";    /* label file extension */
static char * hmmDir = NULL;     /* directory to look for hmm def files */
static char * hmmExt = NULL;     /* hmm def file extension */
static char * newDir = NULL;     /* directory to store new hmm def files */
static char * newExt = NULL;     /* extension of new reestimated hmm files */
static char * newMacExt =NULL;
static char * statFN;            /* stats file, if any */
static char * posDir[5];	 /* directory to store posterior files */
static int numPosDirs=0;
static char * posExt = "pos";	 /* extension of posterior files */

static float minVar  = 0.0;      /* minimum variance (diagonal only) */
static float mixWeightFloor=0.0; /* Floor for mixture weights */
static int minEgs    = 3;        /* min examples to train a model */
static UPDSet uFlags = (UPDSet) (UPMEANS|UPVARS|UPTRANS|UPMIXES); /* update flags */
static int parMode   = -1;       /* enable one of the // modes */
static Boolean stats = FALSE;    /* enable statistics reports */
static char * mmfFn  = NULL;     /* output MMF file, if any */
static int trace     = 0;        /* Trace level */
static Boolean saveBinary = FALSE;  /* save output in binary  */
static Boolean ldBinary = TRUE;        /* load/dump in binary */
static FileFormat dff=UNDEFF;       /* data file format */
static FileFormat lff=UNDEFF;       /* label file format */
static int updateMode = UPMODE_UPDATE; /* dump summed accs, update models or do both? */


static ConfParam *cParm[MAXGLOBS];   /* configuration parameters */
static int nParm = 0;               /* total num params */

static Boolean al_hmmUsed = FALSE;   /* Set for 2-model ReEstimation */
static char al_hmmDir[MAXFNAMELEN];  /* dir to look for alignment hmm defs */
static char al_hmmExt[MAXSTRLEN];  	 /* alignment hmm def file extension */
static char al_hmmMMF[MAXFNAMELEN];  /* alignment hmm MMF */
static char al_hmmLst[MAXFNAMELEN];  /* alignment hmm list */
static char up_hmmMMF[MAXFNAMELEN];  /* alignment hmm list */

static HMMSet al_hset ;      	 /* Option 2nd set of models for alignment */

/* Global Data Structures - valid for all training utterances */
static LogDouble pruneInit = NOPRUNE;    /* pruning threshold initially */
static LogDouble pruneInc = 0.0;         /* pruning threshold increment */
static LogDouble pruneLim = NOPRUNE;     /* pruning threshold limit */
static float minFrwdP = NOPRUNE;         /* mix prune threshold */


static Boolean firstTime = TRUE;    /* Flag used to enable creation of ot */
static Boolean twoDataFiles = FALSE; /* Enables creation of ot2 for FB
                                        training using two data files */
static long totalT=0;       /* total number of frames in training data */
static LogDouble totalPr=0;   /* total log prob upto current utterance */
static Vector vFloor[SMAX]; /* variance floor - default is all zero */

static MemHeap hmmStack;   /*For Storage of all dynamic structures created...*/
static MemHeap uttStack;
static MemHeap fbInfoStack;
static MemHeap accStack;

/* information about transforms */
static XFInfo xfInfo;
static int maxSpUtt = 0;
static float varFloorPercent = 0;

static char *labFileMask = NULL;
static Boolean fMPE=FALSE;
/*TVWR*/
static short posWidth[5];
static short posize=0;
/*static Boolean tied=FALSE;*/

static char *dnn_label_list=NULL;

static int numClust=0;

static Boolean trueRwgtUpOnly=FALSE;
static int numThreads=1;
static Boolean NAT=FALSE;
static Boolean xNdir=FALSE;

/*store a copy of original MMF*/
static Matrix x_mu_bk=NULL;
static Matrix x_cv_bk=NULL;
static int* negs_bk=NULL;

static Matrix lmus_x=NULL;
static Matrix* lcvs_x=NULL;
static DMatrix* lGs=NULL;

static char* nlist=NULL;
static FILE* nlistFile=NULL;
static NoiseModel* nlistHead=NULL;
static NoiseModel* nlistTail=NULL;
static Boolean nlistWrite=FALSE;

static SpkrModel* slistHead=NULL;
static SpkrModel* slistTail=NULL;


static char* lnlist=NULL;
static FILE* lnlistFile=NULL;
static NoiseModel* lnlistHead=NULL;
static NoiseModel* lnlistTail=NULL;
static Boolean lnlistWrite=FALSE;
static Boolean upNoiseModel=FALSE;/*default update noise model*/
static Boolean upAcousModel=FALSE;
static int trimNum=20;
static Boolean tiedTVWR=FALSE;
static int NatEMiters=4;
static float glrate=1.0;
static int acc_subK=99999;
static char* pr_file=NULL;
Boolean NATUpdateNoiseModel(HMMSet *hset, NoiseModel* nm);
Boolean NATUpdateLNoiseModel(HMMSet *hset, NoiseModel* lnm);
void NATUpdateAcousticModel(HMMSet *hset, NoiseModel *nlistHead);
void NATUpdateLAcousticModel(HMMSet *hset);

/*used for SNAT training, track current speaker MLLR mean transform*/
static LinXForm* mllrMean=NULL;
/* ------------------ Process Command Line -------------------------- */
   
/* SetConfParms: set conf parms relevant to HCompV  */
void SetConfParms(void)
{
   int i;
   Boolean b;
   double f;
   char buf[MAXSTRLEN];
   
   nParm = GetConfig("HEREST", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfFlt(cParm,nParm,"VARFLOORPERCENTILE",&f)) varFloorPercent = f;
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&b)) saveBinary = b;
      if (GetConfBool(cParm,nParm,"BINARYACCFORMAT",&b)) ldBinary = b;
      /* 2-model reestimation alignment model set */
      if (GetConfStr(cParm,nParm,"ALIGNMODELMMF",buf)) {
          strcpy(al_hmmMMF,buf); al_hmmUsed = TRUE;
      }
      if (GetConfStr(cParm,nParm,"ALIGNHMMLIST",buf)) {
          strcpy(al_hmmLst,buf); al_hmmUsed = TRUE;
      }
      /* allow multiple individual model files */
      if (GetConfStr(cParm,nParm,"ALIGNMODELDIR",buf)) {
          strcpy(al_hmmDir,buf); al_hmmUsed = TRUE;
      }
      if (GetConfStr(cParm,nParm,"ALIGNMODELEXT",buf)) {
          strcpy(al_hmmExt,buf); al_hmmUsed = TRUE;
      }
      if (GetConfStr(cParm,nParm,"ALIGNXFORMEXT",buf)) {
         xfInfo.alXFormExt = CopyString(&hmmStack,buf);
      }
      if (GetConfStr(cParm,nParm,"ALIGNXFORMDIR",buf)) {
         xfInfo.alXFormDir = CopyString(&hmmStack,buf);
      }
      if (GetConfStr(cParm,nParm,"INXFORMMASK",buf)) {
         xfInfo.inSpkrPat = CopyString(&hmmStack,buf);
      }
      if (GetConfStr(cParm,nParm,"PAXFORMMASK",buf)) {
         xfInfo.paSpkrPat = CopyString(&hmmStack,buf);
      }
      if (GetConfStr(cParm,nParm,"LABFILEMASK",buf)) {
         labFileMask = (char*)malloc(strlen(buf)+1); 
         strcpy(labFileMask, buf);
      }

      if (GetConfStr(cParm,nParm,"UPDATEMODE",buf)) {
         if (!strcmp (buf, "DUMP")) updateMode = UPMODE_DUMP;
         else if (!strcmp (buf, "UPDATE")) updateMode = UPMODE_UPDATE;
         else if (!strcmp (buf, "BOTH")) updateMode = UPMODE_BOTH;
         else HError(2319, "Unknown UPDATEMODE specified (must be DUMP, UPDATE or BOTH)");
      }
      if (GetConfBool(cParm,nParm,"FMPE",&b)){
    	  fMPE=b;
      }
      if (GetConfInt(cParm,nParm,"POSIZE",&i)){
    	  posize=i;
    	  posWidth[1]=i;
    	  posWidth[0]=1;
      }
      if (GetConfInt(cParm,nParm,"POSIZE1",&i)){
    	  posize+=i;
    	  posWidth[1]=i;
    	  posWidth[0]=1;
      }
      if (GetConfInt(cParm,nParm,"POSIZE2",&i)){
    	  posize+=i;
    	  posWidth[2]=i;
    	  posWidth[0]=2;
      }
      if (GetConfInt(cParm,nParm,"POSIZE3",&i)){
    	  posize+=i;
    	  posWidth[3]=i;
    	  posWidth[0]=3;
      }
      if (GetConfInt(cParm,nParm,"POSIZE4",&i)){
    	  posize+=i;
    	  posWidth[4]=i;
    	  posWidth[0]=4;
      }
   }
}

void ReportUsage(void)
{
   printf("\nUSAGE: HERest [options] hmmList dataFiles...\n\n");
   printf(" Option                                       Default\n\n");
   printf(" -a      Use an input linear transform        off\n");
   printf(" -c f    Mixture pruning threshold            10.0\n");
   printf(" -d s    dir to find hmm definitions          current\n");
   printf(" -h s    set output speaker name pattern   *.%%%%%%\n");
   printf("         to s, optionally set input and parent patterns\n");
   printf(" -l N    set max files per speaker            off\n");
   printf(" -m N    set min examples needed per model    3\n");
   printf(" -o s    extension for new hmm files          as src\n");
   printf(" -p N    set parallel mode to N               off\n");
   printf(" -r      Enable Single Pass Training...       \n");
   printf("         ...using two parameterisations       off\n");
   printf(" -s s    print statistics to file s           off\n");
   printf(" -t f [i l] set pruning to f [inc limit]      inf\n");
   printf(" -u tmvwapsr  update t)rans m)eans v)ars w)ghts tmvw\n");
   printf("                a)daptation xform p)rior used     \n");
   printf("                s)semi-tied xform r)rwght         \n");
   printf(" -v f    set minimum variance to f            0.0\n");
   printf(" -w f    set mix weight floor to f*MINMIX     0.0\n");
   printf(" -x s    extension for hmm files              none\n");
   printf(" -z s    Save all xforms to TMF file s        TMF\n");
   printf(" -nomixup default '-u r' will also update mixwgts \n");
   PrintStdOpts("BEFGHIJKLMSTX");
   printf("\n\n");
}

void ResetTVWRCache(HMMSet *hset){
	HMMScanState hss;
	MixtureElem* me;
	NewHMMScan(hset,&hss);
	while(GoNextMix(&hss,FALSE)){
		me=hss.me;
		if(me->rwgt)
			me->rwgt->cached_t=-1;
	}
	EndHMMScan(&hss);
}

void SetuFlags(void)
{
   char *s;
   
   s=GetStrArg();
   uFlags=(UPDSet) 0;        
   while (*s != '\0')
      switch (*s++) {
      case 't': uFlags = (UPDSet) (uFlags+UPTRANS); break;
      case 'm': uFlags = (UPDSet) (uFlags+UPMEANS); break;
      case 'v': uFlags = (UPDSet) (uFlags+UPVARS); break;
      case 'w': uFlags = (UPDSet) (uFlags+UPMIXES); break;
      case 's': uFlags = (UPDSet) (uFlags+UPSEMIT); break;
      case 'a': uFlags = (UPDSet) (uFlags+UPXFORM); break;
      case 'p': uFlags = (UPDSet) (uFlags+UPMAP); break;
      case 'r': uFlags = (UPDSet) (uFlags+UPRWGTS); break;
      default: HError(2320,"SetuFlags: Unknown update flag %c",*s);
         break;
      }
}

/* ScriptWord: return next word from script */
char *ScriptWord(FILE *script, char *scriptBuf)
{
   int ch,qch,i;
   
   i=0; ch=' ';
   while (isspace(ch)) ch = fgetc(script);
   if (ch==EOF) {
      scriptBuf=NULL;
      return NULL;
   }
   if (ch=='\'' || ch=='"'){
      qch = ch;
      ch = fgetc(script);
      while (ch != qch && ch != EOF) {
         scriptBuf[i++] = ch; 
         ch = fgetc(script);
      }
      if (ch==EOF)
         HError(5051,"ScriptWord: Closing quote missing in script file");
   } else {
      do {
         scriptBuf[i++] = ch; 
         ch = fgetc(script);
      }while (!isspace(ch) && ch != EOF);
   }
   scriptBuf[i] = '\0';

   return scriptBuf;
}

void CheckUpdateSetUp()
{
  AdaptXForm *xf;

  xf = xfInfo.paXForm;
  if ((xfInfo.paXForm != NULL) && !(uFlags&UPXFORM)) {
    while (xf != NULL) {
       if ((xf->xformSet->xkind != CMLLR) && (xf->xformSet->xkind != SEMIT))
	HError(999,"SAT only supported with SEMIT/CMLLR transforms");
      xf = xf->parentXForm;
    }
  }
}

int main(int argc, char *argv[]){
   char *datafn=NULL;
   char *datafn2=NULL;
   char *s;
   char *scriptFile;
   char datafn1[MAXSTRLEN];
   char newFn[MAXSTRLEN];
   char posfn[MAXSTRLEN];
   FILE *f;
   UttInfo *utt;            /* utterance information storage */
   FBInfo *fbInfo;          /* forward-backward information storage */
   HMMSet hset;             /* Set of HMMs to be re-estimated */
   Source src;
   float tmpFlt;
   int tmpInt;
   int numUtt,spUtt=0;
   char* hmmListFn;
   HMMSet synrhset;
   char synrListFn[MAXFNAMELEN];
   Boolean useGmmSynr=FALSE;
   Boolean useTiedStatePos=FALSE;
   HLink* phnSynr=NULL;
   TreeNode tspSynr=NULL;
   long preskip;
   int posStream;
   int vSize=0, lvSize=0;
   posDir[1]=NULL;

   void Initialise(FBInfo *fbInfo, MemHeap *x, HMMSet *hset, char *hmmListFn);
   float DoForwardBackward(FBInfo *fbInfo, UttInfo *utt, char *datafn, char *datafn2);
   void UpdateModels(HMMSet *hset, ParmBuf pbuf2);
   void StatReport(HMMSet *hset);

   if(InitShell(argc,argv,herest_version,herest_vc_id)<SUCCESS)
      HError(2300,"HERest: InitShell failed");
   InitMem();    InitMath();
   InitSigP();   InitAudio();
   InitWave();   InitVQ();
   InitLabel();  InitModel();
   if(InitParm()<SUCCESS)  
      HError(2300,"HERest: InitParm failed");
   InitTrain();
   InitUtil();   InitFB();
   InitAdapt(&xfInfo); InitMap();

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(0);
   al_hmmDir[0] = '\0'; al_hmmExt[0] = '\0'; 
   al_hmmMMF[0] = '\0'; al_hmmLst[0] = '\0'; 
   up_hmmMMF[0] = '\0';
   CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   SetConfParms(); 
   CreateHMMSet(&hset,&hmmStack,TRUE);
   CreateHMMSet(&synrhset,&hmmStack,TRUE);
   CreateHeap(&uttStack,   "uttStore",    MSTAK, 1, 0.5, 100,   1000);
   utt = (UttInfo *) New(&uttStack, sizeof(UttInfo));
   CreateHeap(&fbInfoStack,   "FBInfoStore",  MSTAK, 1, 0.5, 100 ,  1000 );
   fbInfo = (FBInfo *) New(&fbInfoStack, sizeof(FBInfo));
   CreateHeap(&accStack,   "accStore",    MSTAK, 1, 1.0, 50000,   500000);

   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      /*if (strlen(s)!=1)
         HError(2319,"HERest: Bad switch %s; must be single letter",s);*/
      switch(s[0]){
      case 'b':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: script file expected");
         scriptFile = GetStrArg(); break;
      case 'c':
    	 if(!strcmp(s,"clr")){
    		 if(NextArg()!=INTARG)
    			 HError(1,"HERest: num of clusters expected");
    		 numClust=GetChkedInt(1,10000,s); break;
    		 break;
    	 }
         minFrwdP = GetChkedFlt(0.0,1000.0,s);
         break;
      case 'd':
    	  if (!strcmp(s, "dll")){
    		  dnn_label_list=GetStrArg();
    		  break;
    	  }
    	  if (NextArg()!=STRINGARG)
    		  HError(2319,"HERest: HMM definition directory expected");
    	  hmmDir = GetStrArg();
    	  break;
      case 'm':
         minEgs = GetChkedInt(0,1000,s); break;
      case 'o':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: HMM file extension expected");
         newExt = GetStrArg(); break;
      case 'p':
    	 if(!strcmp(s,"pr")){
    		 pr_file=GetStrArg();
    		 break;
    	 }
         parMode = GetChkedInt(0,500,s); break;
      case 'r':
    	 if(!strcmp(s,"ronly")){
    		 trueRwgtUpOnly=TRUE;break;
    	 }
         twoDataFiles = TRUE; break;
      case 's':
         stats = TRUE;
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: Stats file name expected");
         statFN = GetStrArg(); break;
      case 't':
    	  if (!strcmp(s,"tiedStatePos")){
    		  useTiedStatePos=TRUE;
    		  break;
    	  }
    	  if (!strcmp(s,"tiedTVWR")){
    		  tiedTVWR=TRUE;
    		  break;
    	  }
    	  pruneInit =  GetChkedFlt(0.0,1.0E20,s);
    	  if (NextArg()==FLOATARG || NextArg()==INTARG)
    	  {
    		  pruneInc = GetChkedFlt(0.0,1.0E20,s);
    		  pruneLim = GetChkedFlt(0.0,1.0E20,s);
    	  }
    	  else
    	  {
    		  pruneInc = 0.0;
    		  pruneLim = pruneInit  ;
    	  }
         break;
      case 'g':
    	  glrate = GetChkedFlt(0.0,1.0,s);
    	  break;
      case 'n':
    	  if (!strcmp(s,"nt")){
    		  if(NextArg()!=INTARG){
    			  HError(2019,"HERest: Number of threads expected");
    		  }
    		  numThreads=GetIntArg();
    	  }
    	  if (!strcmp(s,"nList")){/*noise per utterance*/
    		  /*nlistFile=GetStrArg();*/
    		  nlist=GetStrArg();
    	  }
    	  break;
      case 'k':
		  if(NextArg()!=INTARG){
			  HError(2019,"HERest: Number of examples expected");
		  }
    	  acc_subK=GetIntArg();
    	  break;
      case 'u':
    	  if(!strcmp(s,"upNmodel")){
    		  upNoiseModel=TRUE;
    		  upAcousModel=FALSE;
    		  break;
    	  }
    	  if(!strcmp(s,"upAmodel")){
    		  upAcousModel=TRUE;
    		  upNoiseModel=FALSE;
    		  break;
    	  }
    	  if(!strcmp(s, "u"))
    		  SetuFlags();
    	  else{
    		  HError(1,"HERest: Unknown switch %s",s);
    	  }
    	  break;
      case 'v':
         minVar = GetChkedFlt(0.0,10.0,s); break;
      case 'w':
         mixWeightFloor = MINMIX * GetChkedFlt(0.0,10000.0,s); 
         break;
      case 'x':
    	 if (!strcmp(s,"xNdir")){
    		 xNdir=TRUE;break;
    	 }
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: HMM file extension expected");
         hmmExt = GetStrArg(); break;
      case 'B':
         saveBinary=TRUE;
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(2319,"HERest: Data File format expected");
         if((dff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2389,"HERest: Warning ALIEN Data file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(2319,"HERest: Label File format expected");
         if((lff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2389,"HERest: Warning ALIEN Label file format set");
         break;
      case 'H':
         if (NextArg() != STRINGARG)
            HError(2319,"HERest: HMM macro file name expected");
         strcpy(up_hmmMMF,GetStrArg());
         AddMMF(&hset,up_hmmMMF);
         break;     
      case 'I':
         if (NextArg() != STRINGARG)
            HError(2319,"HERest: MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'L':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: Label file directory expected");
         labDir = GetStrArg(); break;
      case 'M':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: Output macro file directory expected");
         newDir = GetStrArg();
         break;     
      case 'T':
         trace = GetChkedInt(0,0100000,s);
         break;
      case 'X':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: Label file extension expected");
         labExt = GetStrArg(); break;
	 /* additional options for transform support */
      case 'a':
    	  xfInfo.useInXForm = TRUE; break;
      case 'h':
    	  if (NextArg()!=STRINGARG)
    		  HError(1,"Speaker name pattern expected");
    	  xfInfo.outSpkrPat = GetStrArg();
    	  break;
      case 'l':
    	  if (!strcmp(s,"lnList")){
    		  lnlist=GetStrArg();
    		  break;
    	  }
    	  maxSpUtt = GetChkedInt(0,0100000,s);
    	  break;
      case 'E':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: parent transform directory expected");
         xfInfo.usePaXForm = TRUE;
         xfInfo.paXFormDir = GetStrArg(); 
         if (NextArg()==STRINGARG)
        	 xfInfo.paXFormExt = GetStrArg();
         if (NextArg() != SWITCHARG)
        	 HError(2319,"HERest: cannot have -E as the last option");
         break;              
      case 'J':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: input transform directory expected");
         AddInXFormDir(&hset,GetStrArg());
         if (NextArg()==STRINGARG) {
            if (xfInfo.inXFormExt == NULL)
               xfInfo.inXFormExt = GetStrArg(); 
            else
               HError(2319,"HERest: only one input transform extension may be specified");
         }
         if (NextArg() != SWITCHARG)
        	 HError(2319,"HERest: cannot have -J as the last option");
         break;              
      case 'K':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: output transform directory expected");
         xfInfo.outXFormDir = GetStrArg(); 
         if (NextArg()==STRINGARG)
        	 xfInfo.outXFormExt = GetStrArg();
         if (NextArg() != SWITCHARG)
        	 HError(2319,"HERest: cannot have -K as the last option");
         break;
      case 'Y':
    	  if (NextArg()!=STRINGARG)
    		  HError(2319,"HERest: -Y: Posterior directory expected");
    	  posDir[++numPosDirs] = GetStrArg();
    	  if (NextArg()==STRINGARG)
    		  posExt = GetStrArg();
    	 break;
      case 'W':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HERest: -W: Posterior synthesizer model file expected");
    	 AddMMF(&synrhset,GetStrArg());
    	 useGmmSynr=TRUE;
    	 numPosDirs=1;
    	 break;
      case 'e':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HERest: -e: Posterior synthesizer list file expected");
    	 strcpy(synrListFn,GetStrArg());
    	 break;
      case 'N':
    	 if(!strcmp(s,"NAT")){
    		 xfInfo.hasVTS=NAT=TRUE;
    		 if(NextArg()==INTARG){
    			 trimNum=GetIntArg();
    		 }
    	 }
    	 break;
      case 'i':
    	  if(NextArg()!=INTARG){
    		  HError(2319,"HERest: NAT-EM iteration number expected");
    	  }
    	  NatEMiters=GetIntArg();
    	  break;
      case 'z':
         if (NextArg() != STRINGARG)
            HError(2319,"HERest: output TMF file expected");
         xfInfo.xformTMF = GetStrArg(); break;
      default:
         HError(2319,"HERest: Unknown switch %s",s);break;
      }
   } 
   if (NextArg() != STRINGARG)
      HError(2319,"HERest: file name of vocabulary list expected");

   hmmListFn=GetStrArg();
   if(parMode>0)
	   assert(numPosDirs==posWidth[0]);
   if(useGmmSynr){/*load posterior synthesizer */
	   if(useTiedStatePos)/*only work for fMPE */
		   tspSynr=BuildClusteredGMMTree(&synrhset, synrListFn);
	   else/*mono phone synthesizer*/
		   phnSynr=LoadPosSynth(&synrhset,synrListFn);
   }
   Initialise(fbInfo, &fbInfoStack, &hset, hmmListFn);
   NoiseModel* iter=NULL;
   NoiseModel* liter=NULL;

   if(NAT){/*backup the clean back-end and front-end GMM models*/
	   HMMScanState hss;
	   MixPDF* mp;
	   int mixIdx=0, hmmIdx=0;
	   if(nlist){
		   /*********Back up the "clean" acoustic model for compensation*********/
		   vSize=hset.vecSize;
		   x_mu_bk=CreateMatrix(&gstack,hset.numMix,vSize);
		   x_cv_bk=CreateMatrix(&gstack,hset.numMix,vSize);
		   negs_bk=New(&gstack,sizeof(int)*hset.numPhyHMM);
		   NewHMMScan(&hset,&hss);
		   while(GoNextMix(&hss,FALSE)){
				mp=hss.mp;
				mixIdx=mp->mIdx;
				if(mixIdx<=0)continue;
				CopyVector(mp->mean,x_mu_bk[mixIdx]);
				if(hss.mp->ckind==DIAGC)
					CopyVector(mp->cov.var,x_cv_bk[mixIdx]);
				if(hss.mp->ckind==INVDIAGC)
					CopyiVector(mp->cov.var,x_cv_bk[mixIdx]);
		   }
		   EndHMMScan(&hss);

		   hmmIdx=0;
		   NewHMMScan(&hset,&hss);
		   do {
		         negs_bk[hmmIdx++]=0;/*reset the number of examples*/
		   }while(GoNextHMM(&hss));
		   EndHMMScan(&hss);
		   /***********Load the noise model if available*********/
		   nlistFile=fopen(nlist,"r");
		   nlistHead=(NoiseModel*)New(&gstack,sizeof(NoiseModel));
		   nlistHead->next=NULL;
		   if(nlistFile){/*read the stored noise models*/
			   nlistWrite=FALSE;
			   fclose(nlistFile);
			   ReadNoiseModels(nlist, nlistHead, FALSE);
		   }
		   else{
			   nlistWrite=TRUE;
		   }
		   nlistTail=nlistHead;
		   iter=nlistHead->next;
	   }
	   if(lnlist){
		   /*********Back up the "clean" posterior synthesizer model*********/
		   lvSize=synrhset.vecSize;
		   lmus_x=CreateMatrix(&gstack,synrhset.numMix,lvSize);
		   lcvs_x=(Matrix*)New(&gstack,(synrhset.numMix+1)*sizeof(Matrix));
		   lGs=(DMatrix*)New(&gstack,(synrhset.numMix+1)*sizeof(DMatrix));
		   NewHMMScan(&synrhset,&hss);
		   while(GoNextMix(&hss,FALSE)){
			   mp=hss.mp;
			   mixIdx=mp->mIdx;
			   if(mixIdx<=0)continue;
			   CopyVector(mp->mean,lmus_x[mixIdx]);
			   lcvs_x[mixIdx]=CreateMatrix(&gstack,lvSize,lvSize);
			   CovInvert(&gstack, mp->cov.inv,lcvs_x[mixIdx]);/*switch to variance*/
			   lGs[mixIdx]=CreateDMatrix(&gstack,lvSize,lvSize);
		   }
		   EndHMMScan(&hss);
		   /***********Load the noise model if available*********/
		   lnlistFile=fopen(lnlist,"r");
		   lnlistHead=(NoiseModel*)New(&gstack,sizeof(NoiseModel));
		   lnlistHead->next=NULL;
		   if(lnlistFile){/*read the stored noise models*/
			   lnlistWrite=FALSE;
			   fclose(lnlistFile);
			   ReadNoiseModels(lnlist, lnlistHead, FALSE);
		   }
		   else{
			   lnlistWrite=TRUE;
		   }
		   lnlistTail=lnlistHead;
		   liter=lnlistHead->next;
		   /*****************Allocate the noise accumulate space*******/
		   synrhset.upVTS_lAM=upAcousModel;
		   AttachAccs(&synrhset, &accStack, uFlags);
	   }
   }

   InitUttInfo(utt, twoDataFiles);
   numUtt = 1;
   if (trace&T_TOP) 
      SetTraceFB(); /* allows HFB to do top-level tracing */

   do {
      if (NextArg()!=STRINGARG)
         HError(2319,"HERest: data file name expected");
      if (twoDataFiles && (parMode!=0)){
         if ((NumArgs() % 2) != 0)
            HError(2319,"HERest: Must be even num of training files for single pass training");
         strcpy(datafn1,GetStrArg());
         datafn = datafn1;
         
         datafn2 = GetStrArg();
      }else
         datafn = GetStrArg();
      if (parMode==0){
         src=LoadAccs(&hset, datafn,uFlags);
         ReadFloat(&src,&tmpFlt,1,ldBinary);
         totalPr += (LogDouble)tmpFlt;
         ReadInt(&src,&tmpInt,1,ldBinary);
         totalT += tmpInt;
         CloseSource( &src );
         if(lnlist){
        	 datafn[strlen(datafn)-1]='l';/*switch to .acl extension for synthesizer*/
        	 src=LoadAccs(&synrhset, datafn, uFlags);
        	 CloseSource(&src);
         }
      }
      else {
         /* track speakers */	 
    	 if (numUtt==1){
    		 if(NAT&&xfInfo.useOutXForm){
    			 SpkrVTSAccInit(&hset, &xfInfo, datafn);
    		 }
    	 }
         if (UpdateSpkrStats(&hset,&xfInfo, datafn)) {
        	 spUtt=0;
        	 if(NAT&&xfInfo.useOutXForm){
        		 SpkrVTSAccInit(&hset, &xfInfo, datafn);
        	 }
         }
         if(NAT&&xfInfo.useOutXForm){
    		 ZeroAdaptAccs(&hset,xfInfo.outXForm);
         }

         /* Check to see whether set-up is valid */
         CheckUpdateSetUp();
         fbInfo->inXForm = xfInfo.inXForm;
         fbInfo->al_inXForm = xfInfo.al_inXForm;
         fbInfo->paXForm = xfInfo.paXForm;

    	 if(NAT){
        	 char id[MAXSTRLEN];
        	 BaseOf(datafn, id);
        	 if(nlistWrite||lnlistWrite){/*estimate noise model*/
        		iter=nlistWrite?CreateNoiseModel(&gstack, DIAGC, vSize, id, &nlistTail):iter;
        		liter=lnlistWrite?CreateNoiseModel(&gstack, FULLC, lvSize, id, &lnlistTail):liter;
        		EstimateNoiseModel(nlistWrite?iter:NULL, lnlistWrite?liter:NULL, vSize, lvSize, datafn, trimNum);
        	 }
        	 else{/*check name consistency*/
        		 if((iter&&strcmp(iter->id,id))||(liter&&strcmp(liter->id,id))){
        			 HError(1,"Inconsistent noise model name: disk=%s, mem=%s",iter->id,id);
        		 }
        	 }
        	 int iter0,hmmIdx=0,negs=0;
        	 HMMScanState hss;
        	 float totalT_bk, totalPr_bk;
        	 float prv_like=LZERO,cur_like=LZERO,min_like_inc=0.02;
        	 totalT_bk=totalT;
        	 totalPr_bk= totalPr;

        	 for(iter0=0;iter0<NatEMiters;iter0++){
        		 printf("#######################VTS noise estimation EM loop iteration %d#######################\n",iter0);
        		 totalT=totalT_bk;
        		 totalPr=totalPr_bk;
            	 UPDSet uFlags_bk = uFlags;
            	 if(upAcousModel||upNoiseModel){/*if we simply want the adaptation function, there is no need to reset statistics*/
            		 uFlags = (UPDSet) (UPMEANS|UPVARS); /* only conventional mean and variance's occ need to be reset */
            		 ZeroAccs(&hset, uFlags);/*Note that NAT related occ (AtPA,gradient..) is accumulated over utterances, weights and trans occ issue is fixed now 16/01/2013, 19:45:00*/
            		 if(lnlist)
            			 ZeroAccs(&synrhset,uFlags);
            	 }
            	 if(hset.TVWR)
            		 ResetTVWRCache(&hset);
            	 uFlags = uFlags_bk;
            	 if(nlist){/*compensate HMM parameters*/
            		 VtsAdapt(&hset, iter->n_mean, iter->h_mean, iter->cov.var, x_mu_bk, x_cv_bk);
            	 }
            	 if(lnlist){/*compensate synthesizer parameters*/
           			 lVtsAdapt(&synrhset, liter->n_mean, liter->h_mean, liter->cov.inv, lmus_x, lcvs_x, lGs);
            	 }
            	 if (xfInfo.useInXForm ) {/*after the VTS compensation*/
            		 SetXForm(&hset,xfInfo.inXForm,xfInfo.hasVTS);
            		 ApplyHMMSetXForm(&hset,xfInfo.inXForm);
            	 }
                 if(hset.TVWR&&lnlist){
                	 for(posStream=1;posStream<=numPosDirs;posStream++){
                    	 if(posDir[posStream]!=NULL){
                        	 MakeFN(datafn,posDir[posStream],posExt,posfn);
                        	 printf("Reading posterior from file: %s\n", posfn); fflush(stdout);
                    	 }
                    	 if(phnSynr!=NULL){
                    		 printf("Synthesize posterior for file: %s\n", datafn); fflush(stdout);
                    	 }
                    	 if(iter0==1||phnSynr)/*in case of NN posterior, no need to load multiple times*/
                    		 LoadHDF4TVWR(&(fbInfo->hdfInfo), dff, posDir[posStream], posExt, datafn, &synrhset, phnSynr, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL, posStream);
                	 }
                 }
                 if ((maxSpUtt==0) || (spUtt<maxSpUtt)){
                     if(hset.TVWR){
                    	 for(posStream=1;posStream<=numPosDirs;posStream++){
                    		 if(posDir[posStream]!=NULL){
                    			 MakeFN(datafn,posDir[posStream],posExt,posfn);
                    			 printf("Reading posterior from file: %s\n", posfn); fflush(stdout);
                    		 }
                    		 if(phnSynr!=NULL){
                    			 printf("Synthesize posterior for file: %s\n", datafn); fflush(stdout);
                    		 }
                    		 LoadHDF4TVWR(&(fbInfo->hdfInfo), dff, posDir[posStream], posExt, datafn, &synrhset, phnSynr, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL, posStream);
                    	 }
                     }
                    cur_like=DoForwardBackward(fbInfo, utt, datafn, datafn2);
                 }
                 if(iter0==0){
                	 hmmIdx=0;
                	 NewHMMScan(&hset,&hss);
                	 do {
                		 negs=negs_bk[hmmIdx]+(int)hss.hmm->hook;
                		 hss.hmm->hook=(void *)negs;
                		 negs_bk[hmmIdx++]=negs;
                	 }while(GoNextHMM(&hss));
                	 EndHMMScan(&hss);
                 }
                 if(cur_like>=prv_like+min_like_inc){
                	 prv_like=cur_like;
                 }
                 else{
                	 break;
                 }
                 mllrMean=xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL;/*assume speaker transform for joint adaptation is mllr mean*/
                 if(upAcousModel||upNoiseModel){
                     if(nlist){
                    	 if(!NATUpdateNoiseModel(&hset, iter))/*if upAcousModel, then inc is zero and it will stop trying up_noise*/
                    		 break;
                     }
                     if(lnlist){
                    	 NATUpdateLNoiseModel(&synrhset, liter);
                     }
                 }
                 else{/*we simply want the adaptation function*/
                	 break;
                 }
        	 }
    	 }
    	 else{
             if(hset.TVWR){
            	 for(posStream=1;posStream<=numPosDirs;posStream++){
            		 if(posDir[posStream]!=NULL){
            			 MakeFN(datafn,posDir[posStream],posExt,posfn);
            			 printf("Reading posterior from file: %s\n", posfn); fflush(stdout);
            		 }
            		 if(phnSynr!=NULL){
            			 printf("Synthesize posterior for file: %s\n", datafn); fflush(stdout);
            		 }
            		 LoadHDF4TVWR(&(fbInfo->hdfInfo), dff, posDir[posStream], posExt, datafn, &synrhset, phnSynr, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL, posStream);
            	 }
             }
             if(hset.fMPE){
            	 if(useTiedStatePos){
            		 LoadHDF4fMPE2(&(fbInfo->hdfInfo), dff, datafn, tspSynr, hset.TVWR);
            	 }
            	 else{
            		 LoadHDF4fMPE(&(fbInfo->hdfInfo), dff, posfn, datafn, &synrhset, phnSynr, hset.TVWR, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL);
            	 }
             }
             if ((maxSpUtt==0) || (spUtt<maxSpUtt))
                DoForwardBackward(fbInfo, utt, datafn, datafn2);
    	 }
    	 if(!nlistWrite){
    		 iter=nlist?iter->next:NULL;
    		 liter=lnlist?liter->next:NULL;
    	 }
    	 if(hset.TVWR){
    		 ClearRwgtCache(&hset);
    	 }
    	 if(NAT&&xfInfo.useOutXForm){
    		 ResetXFormHMMSet(&hset);
    		 SpkrVTSAccStat(&xfInfo);
    	 }
         numUtt += 1; spUtt++;
      }
   } while (NumArgs()>0);
   if (parMode>0&&(trace&T_TOP)) {
      printf("Reestimation complete - average log prob per frame = %e\n",
             totalPr/totalT);
      printf("     - total frames seen          = %e\n", (double)totalT);
   }
   if(nlist&&upNoiseModel){
	   char fname[256];
	   MakeFN(nlist,newDir,NULL,fname);
	   WriteNoiseModels(xNdir?fname:nlist, nlistHead, FALSE);
   }
   if(lnlist&&upNoiseModel){
	   char fname[256];
	   MakeFN(lnlist,newDir,NULL,fname);
	   WriteNoiseModels(xNdir?fname:lnlist, lnlistHead, FALSE);
   }
   if (uFlags&UPXFORM) {/* ensure final speaker correctly handled */
      UpdateSpkrStats(&hset,&xfInfo, NULL);
      if (trace&T_TOP) {
         printf("Reestimation complete - average log prob per frame = %e (%d frames)\n",
                totalPr/totalT, totalT);
      }
   } else {
      if ((parMode>0 || (parMode==0 && (updateMode&UPMODE_DUMP)))&&!upNoiseModel){
         MakeFN("HER$.acc",newDir,NULL,newFn);
         f=DumpAccs(&hset,newFn,uFlags,parMode);
         tmpFlt = (float)totalPr;
         WriteFloat(f,&tmpFlt,1,ldBinary);
         WriteInt(f,(int*)&totalT,1,ldBinary);
         fclose( f );
         if(lnlist){
        	 MakeFN("HER$.acl",newDir,NULL,newFn);
        	 f=DumpAccs(&synrhset,newFn,uFlags,parMode);
        	 fclose(f);
         }
      }
      if (parMode <= 0) {
    	  if (stats) {
    		  StatReport(&hset);
    	  }
    	  if (nlist&&upAcousModel){/*update speech models, disable mean&variance update flag*/
    		  uFlags = uFlags & ~(UPMEANS|UPVARS);
    		  if(xfInfo.useInXForm){/*add association of noise and speaker models*/
    			  NoiseModel *nmIter=nlistHead;
    			  slistHead=(SpkrModel*)New(&gstack,sizeof(struct _SpkrModel));
    			  slistHead->next=NULL;
    			  slistTail=slistHead;
    			  SpkrModel *smIter=NULL;
    			  char spkr[MAXSTRLEN];
    			  char newMn[MAXSTRLEN];
    			  Boolean maskMatch = FALSE;
    			  while((nmIter=nmIter->next)!=NULL){
    				  maskMatch=MaskMatch(xfInfo.outSpkrPat,spkr,nmIter->id);
    				  if(!maskMatch){
    					  HError(1,"speaker is not recognized from ", nmIter->id);
    				  }
    				  smIter=slistHead;
    				  while((smIter=smIter->next)!=NULL){/*searching existing speaker model*/
    					  if(!strcmp(smIter->id,spkr)){/*found speaker*/
    						  break;
    					  }
    				  }
    				  if(smIter==NULL){/*load a new speaker model*/
    					  smIter=(SpkrModel*)New(&gstack,sizeof(struct _SpkrModel));
    					  strcpy(smIter->id,spkr);
    					  MakeFN(spkr,NULL,xfInfo.inXFormExt,newMn);
    					  smIter->mllrMean = LoadOneXForm(&hset,newMn,NULL)->xformSet->xforms[1];
    					  slistTail->next=smIter;
    					  smIter->next=NULL;
    					  slistTail=smIter;
    				  }
    				  /*add a soft link to noise model*/
    				  nmIter->spkr=smIter;
    			  }
    		  }
    		  NATUpdateAcousticModel(&hset, nlistHead);
    		  UpdateModels(&hset,utt->pbuf2);
    	  }
    	  if (lnlist&&upAcousModel){
    		  uFlags = UPMIXES;
    		  char* synfn=synrhset.mmfNames->fName;
    		  if(strcmp(synfn+strlen(synfn)-3,"syn"))
    			  newMacExt = "syn";
    		  NATUpdateLAcousticModel(&synrhset);
    		  UpdateModels(&synrhset,utt->pbuf2);
    	  }
    	  if (!nlist&&(updateMode&UPMODE_UPDATE))
    		  UpdateModels(&hset,utt->pbuf2);/* update all the model parameters after collecting the sufficient statistics */
      }
   }
   ResetHeap(&uttStack);
   ResetHeap(&fbInfoStack);
   ResetHeap(&hmmStack);
   ResetHeap(&accStack);
   Exit(0);
   return (0);          /* never reached -- make compiler happy */
}



/* -------------------------- Initialisation ----------------------- */

void Initialise(FBInfo *fbInfo, MemHeap *x, HMMSet *hset, char *hmmListFn)
{   
   HSetKind hsKind;
   int L,P,S,vSize,maxM,s,regwsize=0;

   /* Load HMMs and init HMMSet related global variables */
   if(MakeHMMSet( hset, hmmListFn, dnn_label_list)<SUCCESS) /* make the containers for each HMM */
      HError(2321,"Initialise: MakeHMMSet failed");
   regwsize=posize*hset->temp_ctx_size;
   if(regwsize){
	   hset->TVWR=TRUE;
	   hset->regwSize=regwsize;
	   hset->pwidth[0]=posWidth[0];
	   for(s=1;s<=posWidth[0];s++){
		   hset->pwidth[s]=posWidth[s];
	   }
   }
   hset->tiedTVWR=tiedTVWR;
   hset->upVTS_AM=upAcousModel;
   hset->upVTS_NM=upNoiseModel;
   if(LoadHMMSet( hset,hmmDir,hmmExt)<SUCCESS) /* load the parameters into the containers from the specified macro model file */
      HError(2321,"Initialise: LoadHMMSet failed");
   /*SetPrior(hset, pr_file);*/
   SetParmHMMSet(hset);
   SetVFloor(hset, vFloor, 1E-5);
   if (uFlags&UPSEMIT) uFlags = uFlags|UPMEANS|UPVARS;
   AttachAccs(hset, &accStack, uFlags);
   ZeroAccs(hset, uFlags);
   P = hset->numPhyHMM;
   L = hset->numLogHMM;
   vSize = hset->vecSize;
   S = hset->swidth[0];
   maxM = MaxMixInSet(hset);

   hsKind = hset->hsKind;
   if (hsKind==DISCRETEHS)
     uFlags = (UPDSet) (uFlags & (~(UPMEANS|UPVARS|UPXFORM|UPSEMIT)));

   hset->fMPE=fMPE;
   if (parMode != 0) {
      ConvDiagC(hset,TRUE);
   }
   if (trace&T_TOP) {
      if (uFlags&UPMAP)  printf("HERest  MAP Updating: ");
      else printf("HERest  ML Updating: ");
      if (uFlags&UPTRANS) printf("Transitions "); 
      if (uFlags&UPMEANS) printf("Means "); 
      if (uFlags&UPVARS)  printf("Variances "); 
      if (uFlags&UPSEMIT)  printf("SemiTied "); 
      if (uFlags&UPXFORM)  printf("XForms "); 
      if ((uFlags&UPMIXES) && maxM>1)  printf("MixWeights ");
      if ((uFlags&UPRWGTS) && maxM>1)  printf("RegWeights ");
      printf("\n\n ");
    
      if (parMode>=0) printf("Parallel-Mode[%d] ",parMode);

      printf("System is ");
      switch (hsKind){
      case PLAINHS:  printf("PLAIN\n");  break;
      case SHAREDHS: printf("SHARED\n"); break;
      case TIEDHS:   printf("TIED\n"); break;
      case DISCRETEHS: printf("DISCRETE\n"); break;
      }

      printf("%d Logical/%d Physical Models Loaded, VecSize=%d\n",L,P,vSize);
      if (hset->numFiles>0)
         printf("%d MMF input files\n",hset->numFiles);
      if (mmfFn != NULL)
         printf("Output to MMF file:  %s\n",mmfFn); 
      fflush(stdout);
   }
   SetVFloor( hset, vFloor, minVar);
   totalPr = 0.0;

   if (xfInfo.inSpkrPat == NULL) xfInfo.inSpkrPat = xfInfo.outSpkrPat; 
   if (xfInfo.paSpkrPat == NULL) xfInfo.paSpkrPat = xfInfo.outSpkrPat; 
   if (uFlags&UPXFORM) {
      if ((hsKind != PLAINHS) && (hsKind != SHAREDHS))
         HError(999,"Can only estimated transforms with PLAINHS and SHAREDHS!");
      if (uFlags != UPXFORM)
         HError(999,"Can only update linear transforms OR model parameters!");
      xfInfo.useOutXForm = TRUE;
      /* This initialises things - temporary hack - THINK!! */
      CreateAdaptXForm(hset, "tmp");
   } 

   
   /* initialise and  pass information to the forward backward library */
   InitialiseForBack(fbInfo, x, hset, uFlags, pruneInit, pruneInc,
                     pruneLim, minFrwdP);

   if (parMode != 0) {
      ConvLogWt(hset);
   }
   /* 2-model reestimation */
   if (al_hmmUsed){
       if (trace&T_TOP)
           printf("2-model re-estimation enabled\n");
       /* load alignment HMM set */
       CreateHMMSet(&al_hset,&hmmStack,TRUE);
       xfInfo.al_hset = &al_hset;
       if (xfInfo.alXFormExt == NULL) xfInfo.alXFormExt = xfInfo.inXFormExt;
       /* load multiple MMFs */
       if (strlen(al_hmmMMF) > 0 ) {
           char *p,*q;
           Boolean eos;
           p=q=al_hmmMMF;
           for(;;) {
               eos = (*p=='\0');
               if ( ( isspace((int) *p) || *p == '\0' ) && (q!=p) ) {
                   *p='\0';
                   if (trace&T_TOP) { 
                       printf("Loading alignment HMM set %s\n",q);
                   }
                   AddMMF(&al_hset,q);
                   if (eos)
                       break;
                   q=p+1;
               }
               p++;
           }
       }
       if (strlen(al_hmmLst) > 0 ) 
           MakeHMMSet(&al_hset, al_hmmLst,NULL );
       else /* use same hmmList */
           MakeHMMSet(&al_hset, hmmListFn,NULL );
       if (strlen(al_hmmDir) > 0 )
           LoadHMMSet(&al_hset,al_hmmDir,al_hmmExt);
       else
           LoadHMMSet(&al_hset,NULL,NULL);

       /* switch model set */
       UseAlignHMMSet(fbInfo,x,&al_hset);
       if (parMode != 0) {
    	   ConvDiagC(&al_hset,TRUE);
    	   ConvLogWt(&al_hset);
       }

       /* and echo status */
       if (trace&T_TOP) { 
           if (strlen(al_hmmDir) > 0 )
               printf(" HMM Dir %s",al_hmmDir);
           if (strlen(al_hmmExt) > 0 )
               printf(" Ext %s",al_hmmExt);
           printf("\n");
           if (strlen(al_hmmLst) > 0 )
               printf("HMM List %s\n",al_hmmLst);
           printf(" %d Logical/%d Physical Models Loaded, VecSize=%d\n",
                  al_hset.numLogHMM,al_hset.numPhyHMM,al_hset.vecSize);
       }
   }
}

/* ------------------- Statistics Reporting  -------------------- */

/* PrintStats: for given hmm */
void PrintStats(HMMSet *hset,FILE *f, int n, HLink hmm, int numEgs)
{
   WtAcc *wa;
   char buf[MAXSTRLEN];
   StateInfo *si;
   int i,N;
    
   N = hmm->numStates;
   ReWriteString(HMMPhysName(hset,hmm),buf,DBL_QUOTE);
   fprintf(f,"%4d %14s %4d ",n,buf,numEgs);/* HMM access number, s2-occ, s3-occ, s4-occ */
   for (i=2;i<N;i++) {
      si = hmm->svec[i].info;
      wa = (WtAcc *)((si->pdf+1)->hook);/* state.info.pdf+1: the first stream element of this state */
      fprintf(f," %10f",wa->occ);
   }
   fprintf(f,"\n");
}

/* StatReport: print statistics report */
void StatReport(HMMSet *hset)
{
   HMMScanState hss;
   HLink hmm;
   FILE *f;
   int px;

   if ((f = fopen(statFN,"w")) == NULL){
      HError(2311,"StatReport: Unable to open stats file %s",statFN);
      return;
   }
   NewHMMScan(hset,&hss);
   px=1;
   do {
      hmm = hss.hmm;
      PrintStats(hset,f,px,hmm,(int)hmm->hook);
      px++;
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);
   fclose(f);
}

/* -------------------- Top Level of F-B Updating ---------------- */


/* Load data and call FBFile: apply forward-backward to given utterance */
float DoForwardBackward(FBInfo *fbInfo, UttInfo *utt, char * datafn, char * datafn2)
{
   char datafn_lab[MAXFNAMELEN];

   utt->twoDataFiles = twoDataFiles ;
   utt->S = fbInfo->al_hset->swidth[0];

   /* Load the labels - support for label masks */
   if (labFileMask) {
      if (!MaskMatch (labFileMask, datafn_lab, datafn))
         HError(2319,"HERest: LABFILEMASK %s has no match with segemnt %s", labFileMask, datafn);
   }
   else
      strcpy (datafn_lab, datafn);
   LoadLabs(utt, lff, datafn_lab, labDir, labExt);
   /* Load the data */
   LoadData(fbInfo->al_hset, utt, dff, datafn, datafn2);

   if (firstTime) {
      InitUttObservations(utt, fbInfo->al_hset, datafn, fbInfo->maxMixInS); /* Initialise the observation structures within UttInfo */
      firstTime = FALSE;
   }
  
   /* fill the alpha beta and otprobs (held in fbInfo) */
   if (FBFile(fbInfo, utt, datafn)) {
      /* update totals */
      totalT += utt->T ;
      totalPr += utt->pr ;
      /* Handle the input xform Jacobian if necssary */
      if (fbInfo->al_hset->xf != NULL) {
         totalPr += utt->T*0.5*fbInfo->al_hset->xf->xform->det;
      }
   }
   return utt->pr/utt->T;
}

/* --------------------------- Model Update --------------------- */

static int nFloorVar = 0;     /* # of floored variance comps */
static int nFloorVarMix = 0;  /* # of mix comps with floored vars */

/* UpdateTrans: use acc values to calc new estimate for transP */
void UpdateTrans(HMMSet *hset, int px, HLink hmm)
{
   int i,j,N;
   float x,occi;
   TrAcc *ta;
   
   ta = (TrAcc *) GetHook(hmm->transP);
   if (ta==NULL) return;   /* already done */
   N = hmm->numStates;
   for (i=1;i<N;i++) {
      occi = ta->occ[i];
      if (occi > 0.0) 
         for (j=2;j<=N;j++) {
            x = ta->tran[i][j]/occi;/* the proportion of ij within occuring of i */
            hmm->transP[i][j] = (x>MINLARG)?log(x):LZERO;
         }
      else
         HError(-2326,"UpdateTrans: Model %d[%s]: no transitions out of state %d",
                px,HMMPhysName(hset,hmm),i);
   }
   SetHook(hmm->transP,NULL);
}

/* FloorMixes: apply floor to given mix set */
void FloorMixes(HMMSet *hset, MixtureElem *mixes, int M, float floor)
{
   float sum,fsum,scale;
   MixtureElem *me;
   int m;
   
   if (hset->logWt == TRUE) HError(999,"FloorMixes requires linear weights");
   sum = fsum = 0.0;
   for (m=1,me=mixes; m<=M; m++,me++) {
      if (MixWeight(hset,me->weight)>floor)
         sum += me->weight;
      else {
         fsum += floor; me->weight = floor;
      }
   }
   if (fsum>1.0) HError(2327,"FloorMixes: Floor sum too large");
   if (fsum == 0.0) return;
   if (sum == 0.0) HError(2328,"FloorMixes: No mixture weights above floor");
   scale = (1.0-fsum)/sum;
   for (m=1,me=mixes; m<=M; m++,me++)
      if (me->weight>floor) me->weight *= scale;
}

/* FloorTMMixes: apply floor to given tied mix set */
void FloorTMMixes(Vector mixes, int M, float floor)
{
   float sum,fsum,scale,fltWt;
   int m;
   
   sum = fsum = 0.0;
   for (m=1; m<=M; m++) {
      fltWt = mixes[m];
      if (fltWt>floor)
         sum += fltWt;
      else {
         fsum += floor;
         mixes[m] = floor;
      }
   }
   if (fsum>1.0) HError(2327,"FloorTMMixes: Floor sum too large");
   if (fsum == 0.0) return;
   if (sum == 0.0) HError(2328,"FloorTMMixes: No mixture weights above floor");
   scale = (1.0-fsum)/sum;
   for (m=1; m<=M; m++){
      fltWt = mixes[m];
      if (fltWt>floor)
         mixes[m] = fltWt*scale;
   }
}

/* FloorDProbs: apply floor to given discrete prob set */
void FloorDProbs(ShortVec mixes, int M, float floor)
{
   float sum,fsum,scale,fltWt;
   int m;
   
   sum = fsum = 0.0;
   for (m=1; m<=M; m++) {
      fltWt = Short2DProb(mixes[m]);
      if (fltWt>floor)
         sum += fltWt;
      else {
         fsum += floor;
         mixes[m] = DProb2Short(floor);
      }
   }
   if (fsum>1.0) HError(2327,"FloorDProbs: Floor sum too large");
   if (fsum == 0.0) return;
   if (sum == 0.0) HError(2328,"FloorDProbs: No probabilities above floor");
   scale = (1.0-fsum)/sum;
   for (m=1; m<=M; m++){
      fltWt = Short2DProb(mixes[m]);
      if (fltWt>floor)
         mixes[m] = DProb2Short(fltWt*scale);
   }
}



void UpdateRWeights(HMMSet *hset, int px, HLink hmm)
{
	int i,s,m,M=0,N,S,j,rSize;
	float occi;
	WtAcc *wa;
	StateElem *se;
	StreamElem *ste;
	MixtureElem *me;
	RWtAcc *rwa;
	float chksum=0;
	N = hmm->numStates;
	se = hmm->svec+2;
	S = hset->swidth[0];
	char *name=HMMPhysName(hset,hmm);
	Boolean silst=FALSE;
	Vector hj=CreateVector(&gstack,hset->regwSize);
	float occj=0;
	/*if(!strcmp(name,"sil")||!strcmp(name,"sp")){
		silst=TRUE;
	}*/
	for (j=2; j<N; j++,se++){
		ste = se->info->pdf+1;
		for (s=1;s<=S; s++,ste++){
			if(silst){
				ZeroVector(hj);
				occj=0;
			}
			wa = (WtAcc *)ste->hook;
			M=ste->nMix;
			if (wa != NULL) {
				occi = wa->occ;/* state occupancy */
				if (occi>0) {
					for(m=1;m<=M;m++){
						me=ste->spdf.cpdf+m;
						rSize=hset->rkind>=HIGH_SPA_REG?hset->confuTbl[me->rwgt->spIdx][0][0]:hset->regwSize;
						if(!trueRwgtUpOnly)
							me->weight=me->sweight=wa->c[m]/occi;/*update static weight*/
						rwa=(RWtAcc *) (me->rwgt->hook);
						if(rwa&&rwa->occ>1){
							occj+=rwa->occ;
							chksum=0;
							for(i=1;i<=rSize;i++){
								me->rwgt->regparm[i]=rwa->hjm[i]/rwa->occ;
								chksum+=me->rwgt->regparm[i];
								hj[i]+=rwa->hjm[i];
							}
						}
						me->rwgt->hook=NULL;
					}
					if(silst&&occj>1){
						for(m=1;m<=M;m++){
							me=ste->spdf.cpdf+m;
							rSize=hset->rkind>=HIGH_SPA_REG?hset->confuTbl[me->rwgt->spIdx][0][0]:hset->regwSize;
							for(i=1;i<=rSize;i++){
								me->rwgt->regparm[i]=hj[i]/occj;
							}
							if(hset->temp_ctx_size>1)
								me->rwgt->regparm[1]=-1;/*marker for silst*/
						}
					}
				}
				ste->hook = NULL;
			}
		}
	}
	FreeVector(&gstack,hj);
}

void UpdateTiedTVWR(HMMSet *hset, int px, HLink hmm)
{
	int i,m,M=0,N,j,s;
	TVWRElem *ze;
	RWtAcc *rwa;
	RegWgt **rwgts;
	N = hmm->numStates;
	ze = hmm->zvec+2;
	Vector hj=CreateVector(&gstack,hset->regwSize);
	float occj=0;
	Boolean silst=FALSE;
	char *name=HMMPhysName(hset,hmm);
	/*if(!strcmp(name,"sil")||!strcmp(name,"sp")){
		silst=TRUE;
	}*/
	for (j=2; j<N; j++,ze++){
		for(s=1;s<=hset->swidth[0];s++){
			if(silst){
				ZeroVector(hj);
				occj=0;
			}
			rwgts=ze->info->tstr[s].rwgts;
			M=ze->info->tstr[s].nMix;
			for(m=1;m<=M;m++){
				rwa=(RWtAcc *) (rwgts[m]->hook);
				if(rwa&&rwa->occ>1){
					occj+=rwa->occ;
					for(i=1;i<=hset->regwSize;i++){
						rwgts[m]->regparm[i]=rwa->hjm[i]/rwa->occ;
						hj[i]+=rwa->hjm[i];
					}
				}
				rwgts[m]->hook=NULL;
			}
			if(silst&&occj>1){
				for(m=1;m<=M;m++){
					for(i=1;i<=hset->regwSize;i++){
						rwgts[m]->regparm[i]=hj[i]/occj;
					}
					if(hset->temp_ctx_size>1)
						rwgts[m]->regparm[1]=-1;/*marker for silst*/
				}
			}
		}
	}
	FreeVector(&gstack,hj);
}

Group* ClusterMixParms(Point* pnts, int numClust, int vSize, int numPnts, int maxIter){
	int i;
	int iterCount;
	double measure, lastmeasure, change;
	Group* grps;
	printf("Begin clustering initialisation using mix up algorithm..\n");fflush(stdout);
	grps=InitKmeansMU(pnts, numClust, vSize, numPnts, FALSE, FALSE);
	lastmeasure=0;
	ResetGroups(grps,numClust);
	for(i=0;i<numPnts;i++){
		lastmeasure+=SearchAndInsert(&gstack, grps, pnts[i], numClust, Userdt, TRUE, TRUE);
	}
	UpdateGroups(&gstack, grps, numClust, numPnts, FALSE);
	lastmeasure/=numPnts;
	printf("Begin K-Means clustering with measure %f\n", lastmeasure);
	for(iterCount=0;iterCount<maxIter;iterCount++){
		ResetGroups(grps, numClust);
		measure=0;
		for(i=0;i<numPnts;i++){
			measure+=SearchAndInsert(&gstack, grps, pnts[i], numClust, Userdt, TRUE, TRUE);
		}
		UpdateGroups(&gstack, grps, numClust, numPnts, FALSE);
		measure/=numPnts;
		change=lastmeasure-measure;
		printf("\tDistortion measure change %f (%f -> %f) at iteration %d\n",change,lastmeasure,measure,iterCount);fflush(stdout);
		lastmeasure=measure;
		if(change==0)break;
	}
	return grps;
}

void RegWgtCluster(HMMSet* hset, int numClust){
	int i,j,k;
	char type = 'q';
	char* macName="Q_rw_";
	char buf[255];
	int numPoints;
	int maxIter=100;
	Point* points;
	Group* grps;
	Point iter;
	LabId labid;
	MLink mac;
	Point plist;
	RegWgt* rwgt0;
	float chksum=0;
	float initMeasure=0;
	plist=CreatePoint(&gstack,NULL,NULL,NULL);
	numPoints=LoadRegWgtList(hset,plist);
	if (numPoints<=numClust){
		HError(1,"RegWgtCluster: requested cluster num %d is >= the total availabe regweight num %d",numClust, numPoints);
	}
	else{
		printf("Expecting %d clusters out of %d regweights\n",numClust,numPoints);fflush(stdout);
	}
	points=(Point*)New(&gstack,sizeof(Point)*numPoints);
	iter=plist->next;
	for(i=0;i<numPoints;i++,iter=iter->next){
		points[i]=iter;
		Vector hjm=(Vector) (iter->hook);
		for(j=1;j<=hset->regwSize;j++){
			initMeasure-=hjm[j]*(iter->value[j]<=0?LZERO:log(iter->value[j]));
		}
	}
	initMeasure/=numPoints;
	printf("Initial measure %f\n",initMeasure);
	grps=ClusterMixParms(points,numClust,hset->regwSize,numPoints,maxIter);
	printf("Begin cluster mean using maximising the objective function\n");
	for(i=0;i<numClust;i++){
		printf("\tCluster %d has %d items\n",i,(int)grps[i]->hardOcc);fflush(stdout);
		if(grps[i]->hardOcc==0)continue;
		sprintf(buf,"%s%x",macName,i);
		labid=GetLabId(buf,TRUE);
        RegWgt* rwgt=CreateRegWgt(hset,1,-1);
        chksum=0;
        for(k=1;k<=hset->regwSize;k++){
        	rwgt->regparm[k]=grps[i]->mean[k];
        	chksum+=rwgt->regparm[k];
        }
		IncUse(rwgt->regparm);
		mac=NewMacro(hset,0,type,labid,rwgt);
		SetMacroUse(mac,1);
		iter=grps[i]->plist->next;
		for(;iter!=NULL;iter=iter->next){
			((MixtureElem*)iter->owner)->rwgt=rwgt;
		}
	}
}


/* UpdateWeights: use acc values to calc new estimate of mix weights */
void UpdateWeights(HMMSet *hset, int px, HLink hmm)
{
   int i,s,m,M=0,N,S;
   float x,occi;
   WtAcc *wa;
   StateElem *se;
   StreamElem *ste;
   MixtureElem *me;
   HSetKind hsKind;
   N = hmm->numStates;
   se = hmm->svec+2;
   hsKind = hset->hsKind;
   S = hset->swidth[0];
   for (i=2; i<N; i++,se++){
      ste = se->info->pdf+1;
      for (s=1;s<=S; s++,ste++){
         wa = (WtAcc *)ste->hook;
         switch (hsKind){
         case TIEDHS:
            M=hset->tmRecs[s].nMix;
            break;
         case DISCRETEHS:
         case PLAINHS:
         case SHAREDHS:
            M=ste->nMix;
            break;
         }
         if (wa != NULL) {
            occi = wa->occ;/* state occupancy */
            if (occi>0) {
               for (m=1; m<=M; m++){
                  x = wa->c[m]/occi;/* occ(state=i,component=m)/occ(state=i) */
                  if (x>1.0){
                     if (x>1.001)
                        HError(2393,"UpdateWeights: Model %d[%s]: mix too big in %d.%d.%d %5.5f",
                               px,HMMPhysName(hset,hmm),i,s,m,x);
                     x = 1.0;
                  }
                  switch (hsKind){
                  case TIEDHS:
                     ste->spdf.tpdf[m] = (x>MINMIX) ? x : 0.0;
                     break;
                  case DISCRETEHS:
                     ste->spdf.dpdf[m]=(x>MINMIX) ? DProb2Short(x) : DLOGZERO;
                     break;
                  case PLAINHS:
                  case SHAREDHS:
                     me=ste->spdf.cpdf+m;
                     me->weight = (x>MINMIX) ? x : 0.0;
                     break;
                  }
               }
               if (mixWeightFloor>0.0){
                  switch (hsKind){
                  case DISCRETEHS:
                     FloorDProbs(ste->spdf.dpdf,M,mixWeightFloor);
                     break;
                  case TIEDHS:
                     FloorTMMixes(ste->spdf.tpdf,M,mixWeightFloor);
                     break;
                  case PLAINHS:
                  case SHAREDHS:
                     FloorMixes(hset,ste->spdf.cpdf+1,M,mixWeightFloor);
                     break;
                  }
               }
            }else
               HError(-2330,"UpdateWeights: Model %d[%s]: no use of mixtures in %d.%d",
                      px,HMMPhysName(hset,hmm),i,s);
            ste->hook = NULL;
         }
      }
   }
}
      
/* UpdateMeans: use acc values to calc new estimate of means */
void UpdateMeans(HMMSet *hset, int px, HLink hmm)
{
   int i,s,m,k,M,N,S,vSize;
   float occim;
   MuAcc *ma;
   StateElem *se;
   StreamElem *ste;
   MixtureElem *me;
   Vector mean;
   
   N = hmm->numStates;
   se = hmm->svec+2;
   S = hset->swidth[0];
   for (i=2; i<N; i++,se++){
      ste = se->info->pdf+1;
      for (s=1;s<=S;s++,ste++){
         /* nuisance dimensions not updated */
         vSize = hset->swidth[s]-hset->projSize;
         me = ste->spdf.cpdf + 1; M = ste->nMix;
         for (m=1;m<=M;m++,me++){
            if (me->weight > MINMIX){
               mean = me->mpdf->mean;
               ma = (MuAcc *) GetHook(mean);
               if (ma != NULL){
                  occim = ma->occ;
                  if (occim > 0.0)
                     for (k=1; k<=vSize; k++) 
                       mean[k] += ma->mu[k]/occim;
                  else{
                     M = ste->nMix;
                     HError(-2330,"UpdateMeans: Model %d[%s]: no use of mean %d.%d.%d",
                            px,HMMPhysName(hset,hmm),i,s,m);
                  }
                  SetHook(mean,NULL);
               }
            }
         }
      }
   }
}

/* UpdateTMMeans: use acc values to calc new estimate of means for TIEDHS */
void UpdateTMMeans(HMMSet *hset)
{
   int s,m,k,M,S,vSize;
   float occim;
   MuAcc *ma;
   MixPDF *mpdf;
   Vector mean;
   
   S = hset->swidth[0];
   for (s=1;s<=S;s++){
      vSize = hset->swidth[s];
      M = hset->tmRecs[s].nMix;
      for (m=1;m<=M;m++){
         mpdf = hset->tmRecs[s].mixes[m];
         mean = mpdf->mean;
         ma = (MuAcc *) GetHook(mean);
         if (ma != NULL){
            occim = ma->occ;
            if (occim > 0.0)
               for (k=1; k<=vSize; k++)
                  mean[k] += ma->mu[k]/occim;
            else
               HError(-2330,"UpdateMeans: No use of mean %d in stream %d",m,s);
            SetHook(mean,NULL);
         }
      }
   }
}

/* UpdateVars: use acc values to calc new estimate of variances */
void UpdateVars(HMMSet *hset, int px, HLink hmm)
{
	int i,s,m,k,l,M,N,S,vSize;
	float occim,x,muDiffk,muDiffl;
	Vector minV,mpV;
	VaAcc *va;
	MuAcc *ma;
	StateElem *se;
	StreamElem *ste;
	MixtureElem *me;
	Vector mean;
	Covariance cov;
	Boolean mixFloored,shared;
	N = hmm->numStates;
	se = hmm->svec+2;
	S = hset->swidth[0];
	/*char *name=HMMPhysName(hset,hmm);*/
	for (i=2; i<N; i++,se++){
		ste = se->info->pdf+1;
		for (s=1;s<=S;s++,ste++){
			vSize = hset->swidth[s]-hset->projSize;
			minV = vFloor[s];
			me = ste->spdf.cpdf + 1; M = ste->nMix;
			for (m=1;m<=M;m++,me++){
				if (me->weight > MINMIX){
					if (me->mpdf->vFloor == NULL) mpV=minV;
					else mpV=me->mpdf->vFloor;
					cov = me->mpdf->cov;
					va = (VaAcc *) GetHook(cov.var);
					mean = me->mpdf->mean;
					ma = (MuAcc *) GetHook(mean);
					if (va != NULL){
						occim = va->occ;
						mixFloored = FALSE;
						if (occim > 0.0){
							shared=(GetUse(cov.var)>1 || ma==NULL || ma->occ<=0.0);
							if (me->mpdf->ckind==DIAGC) {
								for (k=1; k<=vSize; k++){
									muDiffk=(shared)?0.0:ma->mu[k]/ma->occ;
									x = va->cov.var[k]/occim - muDiffk*muDiffk;
									if (x<mpV[k]) {
										x = mpV[k];
										nFloorVar++;
										mixFloored = TRUE;
									}
									cov.var[k] = x;
								}
							}
							else { /* FULLC */
		                        for (k=1; k<=vSize; k++){
		                           muDiffk=(shared)?0.0:ma->mu[k]/ma->occ;
		                           for (l=1; l<=k; l++){
		                              muDiffl=(shared)?0.0:ma->mu[l]/ma->occ;
		                              x = va->cov.inv[k][l]/occim - muDiffk*muDiffl;
		                              if (k==l && x<mpV[k]) {
		                                 x = mpV[k];
		                                 nFloorVar++;
		                                 mixFloored = TRUE;
		                              }
		                              cov.inv[k][l] = x;
		                           }
		                        }
		                        if(IsPositiveTriM(cov.inv))
		                        	CovInvert(&gstack, cov.inv,cov.inv);
		                        else{
			                        for (k=1; k<=vSize; k++){
			                           muDiffk=(shared)?0.0:ma->mu[k]/ma->occ;
			                           for (l=1; l<=k; l++){
			                              muDiffl=(shared)?0.0:ma->mu[l]/ma->occ;
			                              x = va->cov.inv[k][l]/occim - muDiffk*muDiffl;
			                              if (k==l && x<mpV[k]) {
			                                 x = mpV[k];
			                                 nFloorVar++;
			                                 mixFloored = TRUE;
			                              }
			                              cov.inv[k][l] = k==l?1/x:0;
			                           }
			                        }
		                        }
		                     }
						}
						else{
							MixtureElem *me2;
							me2 = ste->spdf.cpdf + 1; M = ste->nMix;
							HError(-2330,"UpdateVars: Model %d[%s]: no use of variance %d.%d.%d",
									px,HMMPhysName(hset,hmm),i,s,m);
						}
						if (mixFloored == TRUE) nFloorVarMix++;
						SetHook(cov.var,NULL);
					}
				}
			}
		}
	}
}

/* UpdateTMVars: use acc values to calc new estimate of vars for TIEDHS */
void UpdateTMVars(HMMSet *hset)
{
   int s,m,k,l,M,S,vSize;
   float occim,x,muDiffk,muDiffl;
   Vector minV;
   VaAcc *va;
   MuAcc *ma;
   MixPDF *mpdf;
   Vector mean;
   Covariance cov;
   Boolean mixFloored,shared;
   
   S = hset->swidth[0];
   for (s=1;s<=S;s++){
      vSize = hset->swidth[s];
      minV = vFloor[s];
      M = hset->tmRecs[s].nMix;
      for (m=1;m<=M;m++){
         mpdf = hset->tmRecs[s].mixes[m];
         cov = mpdf->cov;
         va = (VaAcc *) GetHook(cov.var);
         mean = mpdf->mean;
         ma = (MuAcc *) GetHook(mean);
         if (va != NULL){
            occim = va->occ;
            mixFloored = FALSE;
            if (occim > 0.0){
               shared=(GetUse(cov.var)>1 || ma==NULL || ma->occ<=0.0);
               if (mpdf->ckind==DIAGC) {
                  for (k=1; k<=vSize; k++){
                     muDiffk=(shared)?0.0:ma->mu[k]/ma->occ;
                     x = va->cov.var[k]/occim - muDiffk*muDiffk;
                     if (x<minV[k]) {
                        x = minV[k];
                        nFloorVar++;
                        mixFloored = TRUE;
                     }
                     cov.var[k] = x;
                  }
               }
               else { /* FULLC */
                  for (k=1; k<=vSize; k++){
                     muDiffk=(shared)?0.0:ma->mu[k]/ma->occ;
                     for (l=1; l<=k; l++){
                        muDiffl=(shared)?0.0:ma->mu[l]/ma->occ;
                        x = va->cov.inv[k][l]/occim - muDiffk*muDiffl;
                        if (k==l && x<minV[k]) {
                           x = minV[k];
                           nFloorVar++;
                           mixFloored = TRUE;
                        }
                        cov.inv[k][l] = x;
                     }
                  }
                  CovInvert(&gstack, cov.inv,cov.inv);
               }
            }
            else
               HError(-2330,"UpdateTMVars: No use of var %d in stream %d",m,s);
            if (mixFloored == TRUE) nFloorVarMix++;
            SetHook(cov.var,NULL);
         }
      }
   }
}

static  int fltcompare(const void *_i, const void *_j)
{
  const float *i = (const float*)_i;
  const float *j = (const float*)_j;
  if (*i > *j)
    return (1);
  if (*i < *j)
    return (-1);
  return (0);
}


void FloorVars(HMMSet *hset1, int s){
  HMMScanState hss1;
  int vsize;
  int i;
  if(!(hset1->hsKind==PLAINHS || hset1->hsKind==SHAREDHS)){
     HError(1, "Percentile var flooring not supported for this kind of hmm set. (e.g. tied.) should be easy.");
  } else { 
     float **varray;
     int M=0,m=0,floored=0;
     vsize = hset1->swidth[s];
     
     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
        if (hss1.s == s) M++;
     }
     EndHMMScan(&hss1); 

     varray = New(&gstack, sizeof(float*) * (vsize+1));
     for(i=1;i<=vsize;i++) varray[i] = New(&gstack, sizeof(float) * M);

     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
        if (hss1.s == s) {
           int k;
           if(hss1.mp->ckind != DIAGC ) HError(1, "FloorVars expects DIAGC covariances. ");
           
           for(k=1;k<=vsize;k++){
              varray[k][m] = hss1.mp->cov.var[k];
           }
           m++;
        }
     }
     EndHMMScan(&hss1); 
     
     for(i=1;i<=vsize;i++){
        qsort((char *) varray[i], M, sizeof(float), fltcompare);
     }

     if(varFloorPercent <=0 || varFloorPercent >= 100) HError(1, "varFloorPercent should be <100 and >0..");
     

     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
        if (hss1.s == s) {
           int k, Pos = (int)(varFloorPercent*0.01*M);
           for(k=1;k<=vsize;k++){
              if(hss1.mp->cov.var[k] < varray[k][Pos]){
                 hss1.mp->cov.var[k] =  varray[k][Pos];
                 floored++;
              }
           }
        }
     }
     EndHMMScan(&hss1); 
     printf("Floored %d (expected to floor %d)\n", floored, (int)( varFloorPercent * 0.01 * M * vsize));
  }
  FixAllGConsts(hset1);
}


void MLUpdateModels(HMMSet *hset, UPDSet uFlags)
{
   HSetKind hsKind;
   HMMScanState hss;
   HLink hmm;
   int px,n,maxM;
   hsKind = hset->hsKind;
   maxM = MaxMixInSet(hset);

   if (hsKind == TIEDHS){ /* TIEDHS - update mu & var once per HMMSet */
      if (uFlags & UPVARS)
         UpdateTMVars(hset);
      if (uFlags & UPMEANS)
         UpdateTMMeans(hset);
      if (uFlags & (UPMEANS|UPVARS))
         FixAllGConsts(hset);
   }

   NewHMMScan(hset,&hss);
   px=1;
   do {   
      hmm = hss.hmm;
      n = (int)hmm->hook;
      /*if (n<minEgs && !(trace&T_UPD))
         HError(-2331,"UpdateModels: %s[%d] copied: only %d egs\n",
                HMMPhysName(hset,hmm),px,n);*/
      if (trace&T_UPD) {
         if (n<minEgs)
            printf("Model %s[%d] copied: only %d examples\n",
                   HMMPhysName(hset,hmm),px,n);
         else
            printf("Model %s[%d] to be updated with %d examples\n",
                   HMMPhysName(hset,hmm),px,n);
         fflush(stdout);
      }
      if (n>=minEgs && n>0) {/*update per state*/
         if (uFlags & UPTRANS)
            UpdateTrans(hset,px,hmm);
         if (uFlags&UPRWGTS){
        	if(hset->tiedTVWR){
        		if(maxM>1&&!trueRwgtUpOnly)
        			UpdateWeights(hset,px,hmm);
        		UpdateTiedTVWR(hset, px, hmm);
        	}
        	else{
        		UpdateRWeights(hset, px, hmm);
        	}
         }
         if (maxM>1 && (uFlags&UPMIXES)){
        	 UpdateWeights(hset,px,hmm);
         }
         if (hsKind != TIEDHS){
            if (uFlags & UPVARS)
               UpdateVars(hset,px,hmm);
            if (uFlags & UPMEANS)
               UpdateMeans(hset,px,hmm);
            if (uFlags & (UPMEANS|UPVARS)) 
               FixGConsts(hmm);
         }  
      }
      px++;
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);
   /*printf("TVWR update stats:ML update count=%d, skip update count=%d, cur average obj=%f, prv average obj=%f\n",ml_up_count,no_up_count,curAvgObj/ml_up_count,prvAvgObj/ml_up_count);*/
   if (trace&T_TOP) {
      if (nFloorVar > 0)
         printf("Total %d floored variance elements in %d different mixes\n",
                nFloorVar,nFloorVarMix);
      fflush(stdout);
   }
}


/* UpdateModels: update all models and save them in newDir if set,
   new files have newExt if set */
void UpdateModels(HMMSet *hset, ParmBuf pbuf2)
{
   int maxM;
   static char str[100];
   BufferInfo info2;
   char macroname[MAXSTRLEN];

   if (trace&T_UPD){
      printf("Starting Model Update\n"); fflush(stdout);
   }
   if (parMode == -1) {
      ForceDiagC(hset); /* invert the variances again to DIAGC from INVDIAGC */
      ConvExpWt(hset); /* invert all mixture log-weights into weights */
   }
   maxM = MaxMixInSet(hset);

   /* 
      This routine tidies up the semi-tied transform and 
      tidies the model up. The transition and priors are 
      not updated 
   */

   if(numClust>0){
	   RegWgtCluster(hset, numClust);
	   uFlags=0;
   }
   if (uFlags & UPSEMIT) {
      UpdateSemiTiedModels(hset, &xfInfo);
      uFlags = uFlags & ~(UPMEANS|UPVARS);
   }

   if (uFlags & UPMAP)
     MAPUpdateModels(hset, uFlags);
   else {
     MLUpdateModels(hset, uFlags);     
   }
   
   if(varFloorPercent){
      int s;
      printf("Flooring all vars to the %f'th percentile of distribution... ", varFloorPercent);
      for(s=1;s<=hset->swidth[0];s++)
         FloorVars(hset,s);
   }

   if (trace&T_TOP){
      if (mmfFn == NULL)
         printf("Saving hmm's to dir %s\n",(newDir==NULL)?"Current":newDir); 
      else
         printf("Saving hmm's to MMF %s\n",mmfFn);
      fflush(stdout);
   }
   ClearSeenFlags(hset,CLR_ALL);
   if (twoDataFiles){
      if (parMode == 0){
         SetChannel("HPARM2");
         nParm = GetConfig("HPARM2", TRUE, cParm, MAXGLOBS);
         if (GetConfStr(cParm,nParm,"TARGETKIND",str))
            hset->pkind = Str2ParmKind(str);
	 if (GetConfStr(cParm,nParm,"MATTRANFN",str)) {
            /* The transform needs to be set-up */
            hset->xf = LoadInputXForm(hset,NameOf(str,macroname),str);
         }
      } else {
         GetBufferInfo(pbuf2,&info2);
         hset->pkind = info2.tgtPK;
	 hset->xf = info2.xform;
      }
   }
   SaveHMMSet(hset,newDir,newExt,newMacExt,saveBinary);
   if (trace&T_TOP) {
      printf("Reestimation complete - average log prob per frame = %e\n",
             totalPr/totalT);
      printf("     - total frames seen          = %e\n", (double)totalT);
   }
}
void accMu(DMatrix J, DMatrix Jt, DVector y_icv, DMatrix _JtP, DMatrix _JtPJ, TriMat JtPJ, Vector JtPacc1, DVector acc1, double occ){
	double temp;
	int i,j;
	int vSize=NumDRows(J);
	DMatDiag(Jt,y_icv,_JtP);
	DMatMult(_JtP,J,_JtPJ);
	for(i=1;i<=vSize;i++){
		DDotProduct(_JtP[i],acc1,&temp);
		JtPacc1[i]+=temp;
		for(j=1;j<=i;j++){
			JtPJ[i][j]+=occ*_JtPJ[i][j];
		}
	}
}

void laccMu(DMatrix J, DMatrix Jt, DMatrix y_icv, DMatrix _JtP, DMatrix _JtPJ, TriMat JtPJ, Vector JtPacc1, DVector acc1, double occ){
	double temp;
	int i,j;
	int vSize=NumDRows(J);
	DMatMult(Jt,y_icv,_JtP);
	DMatMult(_JtP,J,_JtPJ);
	for(i=1;i<=vSize;i++){
		DDotProduct(_JtP[i],acc1,&temp);
		JtPacc1[i]+=temp;
		for(j=1;j<=i;j++){
			JtPJ[i][j]+=occ*_JtPJ[i][j];
		}
	}
}

void jaccMu(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DVector P, DMatrix XX, DMatrix XXX, DVector acc1, DVector mu_n0, DVector mu_h0, double occ,
		DMatrix A, DMatrix B, DMatrix D, DMatrix K, DVector g, DVector h){
	double temp;
	int i,j;
	int vSize=NumDRows(G);
	for(i=1;i<=vSize;i++){
		DDotProduct(G[i],mu_h0,&temp);
		acc1[i]+=occ*temp;
		DDotProduct(F[i],mu_n0,&temp);
		acc1[i]+=occ*temp;
	}
	DMatDiag(Ft,P,XX);
	DMatMult(XX,G,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			A[i][j]+=occ*XXX[i][j];/*A=\sum_{jm} occ_{jm}*F'PG*/
			K[j][i]+=occ*XXX[i][j];/*K=\sum_{jm} occ_{jm}*G'PF*/
		}
	}
	DMatMult(XX,F,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			B[i][j]+=occ*XXX[i][j];/*B=\sum_{jm} occ_{jm}*F'PF*/
		}
	}
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],acc1,&temp);
		g[i]+=temp;/*g=\sum_{jm} F'P{acc1+occ_{jm}*(G*mu_h0+F*mu_n0)}*/
	}

	DMatDiag(Gt,P,XX);
	DMatMult(XX,G,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			D[i][j]+=occ*XXX[i][j];/*D=\sum_{jm} occ_{jm}*G'PG*/
		}
	}
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],acc1,&temp);
		h[i]+=temp;/*g=\sum_{jm} G'P{acc1+occ_{jm}*(G*mu_h0+F*mu_n0)}*/
	}
}

void jlaccMu(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DMatrix P, DMatrix XX, DMatrix XXX, DVector acc1, DVector mu_n0, DVector mu_h0, double occ,
		DMatrix A, DMatrix B, DMatrix D, DMatrix K, DVector g, DVector h){
	double temp;
	int i,j;
	int vSize=NumDRows(G);
	for(i=1;i<=vSize;i++){
		DDotProduct(G[i],mu_h0,&temp);
		acc1[i]+=occ*temp;
		DDotProduct(F[i],mu_n0,&temp);
		acc1[i]+=occ*temp;
	}
	DMatMult(Ft,P,XX);
	DMatMult(XX,G,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			A[i][j]+=occ*XXX[i][j];/*A=\sum_{jm} occ_{jm}*F'PG*/
			K[j][i]+=occ*XXX[i][j];/*K=\sum_{jm} occ_{jm}*G'PF*/
		}
	}
	DMatMult(XX,F,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			B[i][j]+=occ*XXX[i][j];/*B=\sum_{jm} occ_{jm}*F'PF*/
		}
	}
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],acc1,&temp);
		g[i]+=temp;/*g=\sum_{jm} F'P{acc1+occ_{jm}*(G*mu_h0+F*mu_n0)}*/
	}

	DMatMult(Gt,P,XX);
	DMatMult(XX,G,XXX);
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			D[i][j]+=occ*XXX[i][j];/*D=\sum_{jm} occ_{jm}*G'PG*/
		}
	}
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],acc1,&temp);
		h[i]+=temp;/*g=\sum_{jm} G'P{acc1+occ_{jm}*(G*mu_h0+F*mu_n0)}*/
	}
}

void accCv(DMatrix J, DVector x_cv, DVector y_icv, DVector y_acc2, Vector gradient, TriMat hessian, double occ){
	int p,l,d;
	double temp;
	double _gradient_p;
	int vSize=NumDRows(J);
	for(p=1;p<=vSize;p++){
		_gradient_p=0;
		for(d=1;d<=vSize;d++){
			_gradient_p+=(J[d][p]*J[d][p]*x_cv[p]*y_icv[d])*(occ-y_acc2[d]*y_icv[d]);
		}
		gradient[p]-=0.5*_gradient_p;
		if(hessian){
			for(l=1;l<=p;l++){
				temp=0;
				for(d=1;d<=vSize;d++){
					temp+=(J[d][p]*J[d][p]*x_cv[p]*J[d][l]*J[d][l]*x_cv[l]*y_icv[d]*y_icv[d])*(occ-2*(y_acc2[d]*y_icv[d]));
				}
				hessian[p][l]+=0.5*(temp-(p==l?_gradient_p:0));
			}
		}
	}
}

void laccCv(DMatrix A, DMatrix y_icv, DMatrix y_acc2, TriMat gradient, double occ){
	int i,j,p,q,start4p,start4q;
	double temp,chain0_ij;
	int vSize=NumDRows(A);
	DMatrix PAcc2=CreateDMatrix(&gstack,vSize,vSize);
	DMatrix PAcc2P=CreateDMatrix(&gstack,vSize,vSize);
	DMatMult(y_icv,y_acc2,PAcc2);
	DMatMult(PAcc2,y_icv,PAcc2P);
	for(p=1;p<=vSize;p++){
		for(q=1;q<=p;q++){
			start4p=((p-1)/13)*13;/*nonzero chain1_pq starting row index*/
			start4q=((q-1)/13)*13;/*nonzero chain1_pq starting col index*/
			temp=0;
			for(i=start4q+1;i<=start4q+13;i++){/*calculate trance, only care about nonzero elements, trace=\sum_i \sum_j chain0[i][j]*A[j][p]*A[i][q] */
				for(j=start4p+1;j<=start4p+13;j++){
					chain0_ij=occ*y_icv[i][j]-PAcc2P[i][j];
					temp+=chain0_ij*A[j][p]*A[i][q];/*d[i][i]=\sum_j chain0[i][j]*chain1_pq[j][i], chain1_pq[j][i]=A[j][p]*A[i][q]*/
				}
			}
			gradient[p][q]-=0.5*temp;
		}
	}
	FreeDMatrix(&gstack,PAcc2);
}
void laccCv_pn(DMatrix A, DMatrix y_icv, DMatrix y_acc2, TriMat gradient, TriMat pgradient, TriMat ngradient, double occ){
	int i,j,p,q,start4p,start4q;
	double temp,chain0_ij;
	int vSize=NumDRows(A);
	DMatrix PAcc2=CreateDMatrix(&gstack,vSize,vSize);
	DMatrix PAcc2P=CreateDMatrix(&gstack,vSize,vSize);
	DMatMult(y_icv,y_acc2,PAcc2);
	DMatMult(PAcc2,y_icv,PAcc2P);
	for(p=1;p<=vSize;p++){
		for(q=1;q<=p;q++){
			start4p=((p-1)/13)*13;/*nonzero chain1_pq starting row index*/
			start4q=((q-1)/13)*13;/*nonzero chain1_pq starting col index*/
			temp=0;
			for(i=start4q+1;i<=start4q+13;i++){/*calculate trance, only care about nonzero elements, trace=\sum_i \sum_j chain0[i][j]*A[j][p]*A[i][q] */
				for(j=start4p+1;j<=start4p+13;j++){
					chain0_ij=occ*y_icv[i][j]-PAcc2P[i][j];
					temp+=chain0_ij*A[j][p]*A[i][q];/*d[i][i]=\sum_j chain0[i][j]*chain1_pq[j][i], chain1_pq[j][i]=A[j][p]*A[i][q]*/
				}
			}
			if(gradient){
				gradient[p][q]-=0.5*temp;
			}
			if(pgradient&&ngradient){
				if(temp<0){
					pgradient[p][q]-=0.5*temp;
				}
				else{
					ngradient[p][q]+=0.5*temp;
				}
			}
		}
	}
	FreeDMatrix(&gstack,PAcc2);
}
/*matrix A is either G or F, depending on updating additive noise or channel*/
void upMu(MemHeap* x, TriMat JtPJ, Matrix iJtPJ, Vector JtPacc, DVector cep_mu){
	float step_i;
	int i;
	ModDiagHessian(x, JtPJ,iJtPJ,TRUE,NULL,FALSE);
	/*MatInvert(JtPJ,iJtPJ);*/
	int vSize=TriMatSize(JtPJ);
	for(i=1;i<=vSize;i++){
		DotProduct(iJtPJ[i],JtPacc,&step_i);
		cep_mu[i]=cep_mu[i]+step_i;
	}
}

void jcalMu(DMatrix A, DMatrix B, DMatrix D, DMatrix K, DVector g, DVector h, DVector up_h, DVector up_n){
	int i,j;
	double temp;
	int vSize=DVectorSize(g);
	DMatrix X=CreateDMatrix(&gstack,vSize,vSize);
	DMatrix XX=CreateDMatrix(&gstack,vSize,vSize);
	DMatrix XXX=CreateDMatrix(&gstack,vSize,vSize);
	DVector x=CreateDVector(&gstack,vSize);
	DMatInvert(&gstack, K,X);
	DMatMult(B,X,XX);
	DMatMult(XX,D,XXX);
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],h,&temp);
		x[i]=g[i]-temp;
		for(j=1;j<=vSize;j++){
			XXX[i][j]=A[i][j]-XXX[i][j];
		}
	}
	DMatInvert(&gstack, XXX,XXX);
	for(i=1;i<=vSize;i++){
		DDotProduct(XXX[i],x,&temp);
		up_h[i]=temp;/*mu_h=(A-BK^D)^(g-BK^h), ^:inverse*/
	}

	DMatInvert(&gstack, D,X);
	DMatMult(A,X,XX);
	DMatMult(XX,K,XXX);
	for(i=1;i<=vSize;i++){
		DDotProduct(XX[i],h,&temp);
		x[i]=g[i]-temp;
		for(j=1;j<=vSize;j++){
			XXX[i][j]=B[i][j]-XXX[i][j];
		}
	}
	DMatInvert(&gstack, XXX,XXX);
	for(i=1;i<=vSize;i++){
		DDotProduct(XXX[i],x,&temp);
		up_n[i]=temp;/*mu_n=(B-AD^K)^(g-AD^h)*/
	}
	FreeDMatrix(&gstack,X);
}

void backoffUpdate(DVector work, DVector update, float lrate){
	int i;
	int vSize=DVectorSize(work);
	for(i=1;i<=vSize;i++){
		work[i]=lrate*update[i]+(1.0-lrate)*work[i];
	}
}
void jcalCv(TriMat hessian, Matrix ihessian, Vector gradient, DVector cep_cv, DVector up_cv){
	int i;
	float step_i,lrate=1.0,bound=1;/*we let this bound wider for noise update*/
	int vSize=VectorSize(gradient);
	if(hessian){
		CopyNTriMat(hessian,hessian);/*convert to negative*/
		ModDiagHessian(&gstack, hessian, ihessian, TRUE, NULL, FALSE);
		CopyNMatrix(ihessian,ihessian);
		CopyNTriMat(hessian,hessian);/*remember convert it back, important*/
		/*for(i=1;i<=vSize;i++){
			hessian[i][i]-=1;
		}
		MatInvert(hessian,ihessian);*/
	}
	for(i=1;i<=vSize;i++){
		if(hessian)
			DotProduct(ihessian[i],gradient,&step_i);/*for Newton's update*/
		else
			step_i=-gradient[i];/* for gradient update*/
		step_i=lrate*step_i;
		step_i=step_i>bound?bound:(step_i<-bound?-bound:step_i);
		up_cv[i]=exp(log(cep_cv[i])-step_i);
	}
}

/*used for acoustic model variance update with stabilized bound*/
void upCv(MemHeap* x, TriMat hessian, Matrix ihessian, Vector gradient, DVector cep_cv){
	int i;
	float step_i,bound=1;
	int vSize=TriMatSize(hessian);
	CopyNTriMat(hessian,hessian);
	ModDiagHessian(x, hessian, ihessian, TRUE, NULL, FALSE);
	CopyNMatrix(ihessian,ihessian);
	CopyNTriMat(hessian,hessian);/*remember convert it back, important*/
	/*for(i=1;i<=vSize;i++){
		hessian[i][i]-=1;
	}
	MatInvert(hessian,ihessian);*/
	for(i=1;i<=vSize;i++){
		DotProduct(ihessian[i],gradient,&step_i);
		step_i=step_i>bound?bound:(step_i<-bound?-bound:step_i);
		cep_cv[i]=exp(log(cep_cv[i])-step_i);
	}
}


void lupCv_pn(TriMat pgradient, TriMat ngradient, DMatrix cep_cv, float *lrate){
	int i,j;
	float minlrate;
	int vSize=TriMatSize(pgradient);
	TriMat _cv=CreateTriMat(&gstack,vSize);
    minlrate=1E-4;
    do{
		for(i=1;i<=vSize;i++){
			for(j=1;j<=i;j++){
				_cv[i][j]=cep_cv[i][j]+*lrate*sqrt(cep_cv[i][i])*sqrt(cep_cv[j][j])*(pgradient[i][j]-ngradient[i][j])/(pgradient[i][j]+ngradient[i][j]);
				/*tempX[i][j]=cep_cv[i][j]+curlrate*(pgradient[i][j]-ngradient[i][j]);*/
			}
		}
		if(!IsPositiveTriM(_cv)){
			*lrate/=2;
		}
		else{
			break;
		}
    }while(*lrate>=minlrate);
	if(*lrate>=minlrate)
		Tri2DMat(_cv,cep_cv);
	FreeTriMat(&gstack, _cv);
}



/**
 * Given a new noise model, calculate the auxiliary function
 */
float ReAdaptAndCalcAux(HMMSet *hset, DMatrix dct, DMatrix idct, DVector n_mu_work, DVector h_mu_work, DVector* n_cv_work){
	HMMScanState hss;
	MixPDF* mp;
	int mixIdx;
	int i,j;
	float obj=0;
	float gconst=numCepCoef*3*log(TPI);
	DVector x_mu[3];
	DVector x_cv[3];
	DVector y_mu[3];
	DVector y_cv[3];
	for(j=0;j<3;j++){
		x_mu[j]=CreateDVector(&gstack,numCepCoef);
		x_cv[j]=CreateDVector(&gstack,numCepCoef);
		y_mu[j]=CreateDVector(&gstack,numCepCoef);
		y_cv[j]=CreateDVector(&gstack,numCepCoef);
	}
	DVector exPoint=CreateDVector(&gstack,numCepCoef);
	DVector g0=CreateDVector(&gstack,numCepCoef);
	DMatrix G=CreateDMatrix(&gstack,numCepCoef,numCepCoef);DMatrix Gt=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	DMatrix F=CreateDMatrix(&gstack,numCepCoef,numCepCoef);DMatrix Ft=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	NewHMMScan(hset,&hss);
	while(GoNextMix(&hss,FALSE)){
		mp=hss.mp;
		mixIdx=mp->mIdx;
		MuAcc* ma = (MuAcc *) GetHook(mp->mean);
		VaAcc* va = (VaAcc *) GetHook(mp->cov.var);
		float occ=ma->occ;
		if(occ==0||mixIdx<=0){
			continue;
		}
		for(i=1;i<=numCepCoef;i++){
			for(j=0;j<3;j++){
				x_mu[j][i]=x_mu_bk[mixIdx][numCepCoef*j+i];
				x_cv[j][i]=x_cv_bk[mixIdx][numCepCoef*j+i];
			}
		}
		/********************************Jacobian re-calculation************************************/
		for(i=1;i<=numCepCoef;i++){/*mu=mu_n-mu_x-mu_h*/
			exPoint[i]=n_mu_work[i]-x_mu[0][i]-h_mu_work[i];
		}
		VTS(&gstack, dct, idct, exPoint, G, Gt, F, Ft, g0);
		/*******************************Re-adaptation****************************************/
		for(i=1;i<=numCepCoef;i++){
			y_mu[0][i]=x_mu[0][i]+h_mu_work[i]+g0[i];
			DDotProduct(G[i],x_mu[1],&y_mu[1][i]);
			DDotProduct(G[i],x_mu[2],&y_mu[2][i]);
		}
		for(j=0;j<3;j++){
			dVtsUpVar(G, Gt, F, Ft, x_cv[j], n_cv_work[j], y_cv[j]);
		}
		for(i=1;i<=numCepCoef;i++){
			for(j=0;j<3;j++){
				mp->mean[numCepCoef*j+i]=y_mu[j][i];
				mp->cov.var[numCepCoef*j+i]=mp->ckind==DIAGC?y_cv[j][i]:1/y_cv[j][i];
			}
		}
		if(mllrMean){
			ApplyXForm2Vector(mllrMean,mp->mean);
		}
		mp->ckind==DIAGC?FixDiagGConst(mp):FixInvDiagGConst(mp);
		obj-=occ*(mp->gConst-gconst);
		for(i=1;i<=numCepCoef*3;i++){
			obj-=(va->cov.var[i]-2*ma->mu[i]*mp->mean[i]+occ*mp->mean[i]*mp->mean[i])/(mp->ckind==DIAGC?mp->cov.var[i]:1/mp->cov.var[i]);
		}
	}
	EndHMMScan(&hss);
	FreeDVector(&gstack,x_mu[0]);
	return obj;
}

/**
 * Given the noise model and compensated acoustic model, calculate the new derivatives
 * make sure the model has been compensated
 */
float CalcAuxAndDiff(HMMSet *hset, DMatrix dct, DMatrix idct, DVector n_mu, DVector h_mu, DVector* n_cv,
		DMatrix A, DMatrix B, DMatrix D, DMatrix K, DVector g, DVector h, Vector* gradient, TriMat* hessian, char* accMode, NoiseModel* nm){
	HMMScanState hss;
	MixPDF* mp;
	int mixIdx;
	int i,j;
	float obj=0;
	float gconst=numCepCoef*3*log(TPI);
	int vSize=numCepCoef*3;
	DVector y_acc1=CreateDVector(&gstack,numCepCoef);
	DVector y_acc2[3];
	DVector y_icv[3];
	DVector exPoint=CreateDVector(&gstack,numCepCoef);
	DMatrix G=CreateDMatrix(&gstack,numCepCoef,numCepCoef);DMatrix Gt=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	DMatrix F=CreateDMatrix(&gstack,numCepCoef,numCepCoef);DMatrix Ft=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	DMatrix XX=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	DMatrix XXX=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	DMatrix As=CreateDMatrix(&gstack,numCepCoef,numCepCoef);
	for(j=0;j<3;j++){
		y_icv[j]=CreateDVector(&gstack,numCepCoef);
		y_acc2[j]=CreateDVector(&gstack,numCepCoef);
	}
	if(strchr(accMode,'m')){
		ZeroDMatrix(A);ZeroDMatrix(B);ZeroDMatrix(D);ZeroDMatrix(K);ZeroDVector(g);ZeroDVector(h);
	}
	if(strchr(accMode,'v')){
		for(j=0;j<3;j++){
			ZeroVector(gradient[j]);
			ZeroTriMat(hessian[j]);
		}
	}
	if(mllrMean){/*Copy the Jacobian of static parameters*/
		Mat2DMat(mllrMean->xform[1],As);
	}
	NewHMMScan(hset,&hss);
	while(GoNextMix(&hss,FALSE)){
		mp=hss.mp;
		mixIdx=mp->mIdx;
		SVector mean = mp->mean;/*this is the mean after both noise and speaker adaptation*/
		SVector var=mp->cov.var;
		MuAcc* ma = (MuAcc *) GetHook(mean);
		VaAcc* va = (VaAcc *) GetHook(var);
		float occ=ma->occ;
		if(occ==0||mixIdx<=0){
			continue;
		}
		/*****************************calc obj***********************/
		obj-=occ*(mp->gConst-gconst);
		for(i=1;i<=vSize;i++){/*based on compensated acoustic model*/
			float temp=(va->cov.var[i]-2*ma->mu[i]*mean[i]+occ*mean[i]*mean[i])/(mp->ckind==DIAGC?var[i]:1/var[i]);
			obj-=temp;
		}
		/*****************************calc Jacobian*****************************************/
		for(i=1;i<=numCepCoef;i++){/*mu=mu_n-mu_x-mu_h*/
			exPoint[i]=n_mu[i]-x_mu_bk[mixIdx][i]-h_mu[i];
		}
		VTS(&gstack, dct, idct, exPoint, G, Gt, F, Ft, NULL);
		/*****************************calc NAT statistics for noise model***********************************/
		if(upNoiseModel){
			for(i=1;i<=numCepCoef;i++){
				for(j=0;j<3;j++){
					y_icv[j][i]=(mp->ckind==DIAGC?1/var[numCepCoef*j+i]:var[numCepCoef*j+i]);
					y_acc2[j][i]=va->cov.var[numCepCoef*j+i]-2*ma->mu[numCepCoef*j+i]*mean[numCepCoef*j+i]+occ*mean[numCepCoef*j+i]*mean[numCepCoef*j+i];/*sum_t gamma_jmt^k (o_t^k-u_jmk)*(o_t^k-u_jmk)*/
				}
				y_acc1[i]=ma->mu[i]-occ*mean[i];/*sum_t gamma_jmt^k (o_t^k-u_jmk)*/
			}
			if(strchr(accMode,'v')){
				for(j=0;j<3;j++){
					accCv(F, n_cv[j], y_icv[j], y_acc2[j], gradient[j], hessian[j], occ);
				}
			}
			if(strchr(accMode,'m')){
				if(mllrMean){
					/*G=A*G, F=A*F, note that mp->mean has been jointly adapted*/
					DMatMult(As,G,XX); CopyDMatrix(XX,G); DTranspose(G,Gt);
					DMatMult(As,F,XX); CopyDMatrix(XX,F); DTranspose(F,Ft);
				}
				jaccMu(G, Gt, F, Ft, y_icv[0], XX, XXX, y_acc1, n_mu, h_mu, occ, A, B, D, K, g, h);
			}
		}
		/**********************************calc NAT statistics for speech model************************************/
		if(upAcousModel){
			/*update the example list*/
			if(ma->acc_subK<acc_subK){
				NATAcc* newAcc=(NATAcc*)New(&accStack,sizeof(NATAcc));
				newAcc->occ_k=occ;
				newAcc->mu_k=CreateVector(&accStack,vSize);
				newAcc->cov_k.var=CreateVector(&accStack,vSize);
				CopyVector(ma->mu,newAcc->mu_k);
				CopyVector(va->cov.var,newAcc->cov_k.var);
				newAcc->noise_k=nm;
				NATAcc* temp=ma->acc_k->next;
				ma->acc_k->next=newAcc;
				newAcc->next=temp;
				ma->acc_subK++;
			}
			else{
				NATAcc* accIter=ma->acc_k->next;
				NATAcc* minAcc=accIter;
				while((accIter=accIter->next)!=NULL){
					if(accIter->occ_k<minAcc->occ_k){
						minAcc=accIter;
					}
				}
				if(occ>minAcc->occ_k){
					minAcc->occ_k=occ;
					CopyVector(ma->mu,minAcc->mu_k);
					CopyVector(va->cov.var,minAcc->cov_k.var);
					minAcc->noise_k=nm;
				}
			}
			ma->acc_K++;
		}
	}
	EndHMMScan(&hss);
	FreeDVector(&gstack,y_acc1);
	return obj;
}

/*
 * additive and channel means are jointly updated
 * ApplyXForm2Vector(xformSet->xforms[numXf],mp->mean);
 */

Boolean NATUpdateNoiseModel(HMMSet *hset, NoiseModel* nm){
	Vector memMarker=CreateVector(&gstack,1);
	/*go through acc of each mix*/
	int i,j;
	DMatrix dct,idct;
	DVector n_mu_work, h_mu_work;
	DVector n_cv_work[3];
	DMatrix A, B, D, K;
	DVector g,h,h_mu_update,n_mu_update;
	DVector n_cv_up[3];
	Vector gradient[3];
	TriMat hessian[3];
	Matrix ihessian;
	RestoreAccs(hset);/*offset static statistics*/
	for(j=0;j<3;j++){
		n_cv_work[j]=CreateDVector(&gstack, numCepCoef);
		gradient[j]=CreateVector(&gstack,numCepCoef);
		hessian[j]=CreateTriMat(&gstack,numCepCoef);
		n_cv_up[j]=CreateDVector(&gstack, numCepCoef);
	}
	ihessian=CreateMatrix(&gstack,numCepCoef,numCepCoef);
	A=CreateDMatrix(&gstack,numCepCoef, numCepCoef);B=CreateDMatrix(&gstack,numCepCoef, numCepCoef);D=CreateDMatrix(&gstack,numCepCoef, numCepCoef);
	K=CreateDMatrix(&gstack,numCepCoef, numCepCoef);g=CreateDVector(&gstack,numCepCoef);h=CreateDVector(&gstack,numCepCoef);
	h_mu_update=CreateDVector(&gstack,numCepCoef);n_mu_update=CreateDVector(&gstack,numCepCoef);
	dct=liftDCT(&gstack, numChans, numCepCoef, cepLifter);idct=iLiftDCT(&gstack, numChans,numCepCoef, cepLifter);
	n_mu_work=CreateDVector(&gstack,numCepCoef);h_mu_work=CreateDVector(&gstack,numCepCoef);
	float origobj0=0, origobj=0, curobj=0;
	int k,upNoiseLoop=upAcousModel?1:4;
	int k1,bkOffLoop=6;
	int k0,upVarLoop=5;
	float lrate=glrate;
	float min_obj_inc=1;
	Boolean isUpdated=FALSE;
	for(k=0;k<upNoiseLoop;k++){
		/**********************Load noise model********************************/
		for(i=1;i<=numCepCoef;i++){
			n_mu_work[i]=nm->n_mean[i];
			h_mu_work[i]=nm->h_mean[i];
			for(j=0;j<3;j++)
				n_cv_work[j][i]=nm->cov.var[numCepCoef*j+i];
		}
		/************accumulate stats for NAT noise update**************/
		origobj0=origobj=CalcAuxAndDiff(hset, dct, idct, n_mu_work, h_mu_work, n_cv_work,
				A, B, D, K, g, h, gradient, hessian, "m", nm);
		if(upNoiseModel){
			printf("\tUpdate noise parameter loop iteration %d\n",k);
			/*update means once due to close form update*/
			jcalMu(A, B, D, K, g, h, h_mu_update, n_mu_update);
			lrate=glrate;
			for(k1=0;k1<=bkOffLoop;k1++){
				backoffUpdate(n_mu_work, n_mu_update, lrate);
				backoffUpdate(h_mu_work, h_mu_update, lrate);
				curobj=ReAdaptAndCalcAux(hset, dct, idct, n_mu_work, h_mu_work, n_cv_work);
				if(curobj>origobj+min_obj_inc){
					for(i=1;i<=numCepCoef;i++){
						nm->n_mean[i]=n_mu_work[i];
						nm->h_mean[i]=h_mu_work[i];
					}
					isUpdated=TRUE;
					printf("\t\tUpdate noise means: origobj=%e, curobj=%e, lrate=%e\n",origobj,curobj,lrate);
					origobj=curobj;
					break;
				}
				else{/*reset the working means*/
					for(i=1;i<=numCepCoef;i++){
						n_mu_work[i]=nm->n_mean[i];
						h_mu_work[i]=nm->h_mean[i];
					}
					if(k1==bkOffLoop-1){/*reset and force re-adapt for obj calculation at final iteration, since CalcAuxAndDiff() will assume models are compensated*/
						printf("\t\tNo Update for noise means: origobj=%e, curobj=%e, lrate=%e\n",origobj,curobj,lrate);
						lrate=0;
					}
					lrate*=(k1==0?0.5:lrate);
				}
			}
			/*update variances 10 times due to line search method*/
			for(k0=0;k0<upVarLoop;k0++){
				lrate=glrate;
				origobj=CalcAuxAndDiff(hset, dct, idct, n_mu_work, h_mu_work, n_cv_work,
						A, B, D, K, g, h, gradient, hessian,"v",nm);
				for(j=0;j<3;j++){
					jcalCv(hessian[j], ihessian, gradient[j], n_cv_work[j], n_cv_up[j]);
				}
				for(k1=0;k1<=bkOffLoop;k1++){
					for(j=0;j<3;j++){
						backoffUpdate(n_cv_work[j], n_cv_up[j], lrate);
					}
					curobj=ReAdaptAndCalcAux(hset, dct, idct, n_mu_work, h_mu_work, n_cv_work);
					if(curobj>origobj+min_obj_inc){
						for(i=1;i<=numCepCoef;i++){
							for(j=0;j<3;j++){
								nm->cov.var[numCepCoef*j+i]=n_cv_work[j][i];
							}
						}
						isUpdated=TRUE;
						printf("\t\tUpdate noise variance loop iteration %d\t: origobj=%e, curobj=%e, lrate=%e\n",k0,origobj,curobj,lrate);
						origobj=curobj;
						break;
					}
					else{
						for(i=1;i<=numCepCoef;i++){
							for(j=0;j<3;j++){
								n_cv_work[j][i]=nm->cov.var[numCepCoef*j+i];
							}
						}
						if(k1==bkOffLoop-1){/*reset and force re-adapt for obj calculation at final iteration*/
							printf("\t\tNo Update for noise variance loop iteration %d\t: origobj=%e, curobj=%e, lrate=%e\n",k0,origobj,curobj,lrate);
							lrate=0;
						}
						lrate*=(k1==0?0.5:lrate);
					}
				}
				if(k1==bkOffLoop+1)break;
			}
		}
		if(curobj-origobj0<min_obj_inc){/*check upNoise Loop*/
			break;
		}
	}
	FreeVector(&gstack,memMarker);
	return isUpdated;
}


Boolean NATUpdateLNoiseModel(HMMSet *hset, NoiseModel* lnm){
	HMMScanState hss;
	Vector mu_n, mu_h;
	TriMat cv_n;
	Vector memMarker=CreateVector(&gstack,1);
	mu_n=lnm->n_mean;
	mu_h=lnm->h_mean;
	cv_n=lnm->cov.inv;
	/*go through acc of each mix*/
	int i,j,lnumCeps;
	lnumCeps=hset->vecSize;/*e.g. 13*8 or 13*9 */
	int mixIdx=0;
	int k,KK=10;
	DVector x_mu, n_mu, h_mu, y_mu, n_mu_up, h_mu_up;
	DMatrix n_cv, y_icv, x_cv,y_cv;
	DVector y_acc1;
	DMatrix y_acc2;
	TriMat gradient, pgradient, ngradient;
	DMatrix XXX,XX;
	DMatrix G,Gt,F,Ft,GCx,GCxGt,FCn,FCnFt, A, B, D, K;
	DVector exPoint,g0,g,h;
	MixPDF* mp;
	x_mu=CreateDVector(&gstack,lnumCeps);
	n_mu=CreateDVector(&gstack,lnumCeps);
	h_mu=CreateDVector(&gstack,lnumCeps);
	n_mu_up=CreateDVector(&gstack,lnumCeps);
	h_mu_up=CreateDVector(&gstack,lnumCeps);
	exPoint=CreateDVector(&gstack,lnumCeps);
	g0=CreateDVector(&gstack,lnumCeps);
	y_mu=CreateDVector(&gstack,lnumCeps);
	y_acc1=CreateDVector(&gstack,lnumCeps);
	y_acc2=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	XXX=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	XX=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	gradient=CreateTriMat(&gstack,lnumCeps);
	pgradient=CreateTriMat(&gstack,lnumCeps);
	ngradient=CreateTriMat(&gstack,lnumCeps);
	A=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	B=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	D=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	K=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	g=CreateDVector(&gstack,lnumCeps);
	h=CreateDVector(&gstack,lnumCeps);
	ZeroTriMat(gradient);
	ZeroTriMat(pgradient);
	ZeroTriMat(ngradient);
	y_icv=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	n_cv=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	x_cv=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	y_cv=CreateDMatrix(&gstack,lnumCeps, lnumCeps);
	int trajLen=lnumCeps/numCepCoef;
	DMatrix dct=lLiftDCT(&gstack, trajLen, numChans, numCepCoef, cepLifter);
	DMatrix idct=ilLiftDCT(&gstack, trajLen, numChans, numCepCoef, cepLifter);
	F=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	Gt=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	Ft=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	GCx=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	GCxGt=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	FCn=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	FCnFt=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	/**********************Load noise model********************************/
	for(i=1;i<=lnumCeps;i++){
		n_mu[i]=mu_n[i];
		h_mu[i]=mu_h?mu_h[i]:0;
		for(j=1;j<=i;j++){
			n_cv[j][i]=n_cv[i][j]=cv_n[i][j];
		}
	}
	float origobj0=0,origobj=0, curobj=0;
	float gconst=lnumCeps*log(TPI);
	NewHMMScan(hset,&hss);
	while(GoNextMix(&hss,FALSE)){
		mp=hss.mp;
		mixIdx=mp->mIdx;
		SVector mean = mp->mean;
		STriMat inv=mp->cov.inv;
		MuAcc* ma = (MuAcc *) GetHook(mean);
		VaAcc* va = (VaAcc *) GetHook(inv);
		double occ=ma->occ;
		if(occ<1){/*to speed up*/
			continue;
		}
		for(i=1;i<=lnumCeps;i++){
			y_acc1[i]=ma->mu[i];
			for(j=1;j<=i;j++){
				y_icv[i][j]=y_icv[j][i]=inv[i][j];
				y_acc2[i][j]=y_acc2[j][i]=va->cov.inv[i][j];
			}
		}
		G=lGs[mixIdx];/*Jacobian is calculated based on the previous expansion point*/
		for(i=1;i<=lnumCeps;i++){
			for(j=1;j<=lnumCeps;j++){
				Gt[j][i]=G[i][j];
				Ft[j][i]=F[i][j]=(i==j)?(1-G[i][j]):-G[i][j];
			}
		}
		float logdet=mp->gConst-gconst;
		origobj-=occ*logdet;
		float trace=0;
		for(i=1;i<=lnumCeps;i++){
			for(j=1;j<=lnumCeps;j++){
				trace+=y_acc2[i][j]*y_icv[j][i];
			}
		}
		origobj=origobj-trace;
		for(i=1;i<=lnumCeps;i++){
			ma->mu[i]=ma->mu[i]+occ*mean[i];/*\sum_t y_t-mu_y------------>\sum_t y_t, last bug found here, forget multiplying occ*/
		}
		for(i=1;i<=lnumCeps;i++){
			for(j=1;j<=i;j++){
				va->cov.inv[i][j]=va->cov.inv[i][j]+ma->mu[i]*mean[j]+ma->mu[j]*mean[i]-occ*mean[i]*mean[j];/*\sum_t (y_t-mu_y)^2--------------->\sum_t y_t^2*/
			}
		}
		if(upNoiseModel){
			jlaccMu(G, Gt, F, Ft, y_icv, XX, XXX, y_acc1, n_mu, h_mu, occ, A, B, D, K, g, h);
			laccCv_pn(F, y_icv, y_acc2, gradient, pgradient, ngradient, occ);
		}
		if(upAcousModel){
			laccMu(G, Gt, y_icv, XX, XXX, ma->lGtPG, ma->lGtPacc1, y_acc1, occ);
			laccCv_pn(G, y_icv, y_acc2, NULL, va->lpgradient, va->lngradient, occ);
		}
	}
	EndHMMScan(&hss);
	origobj0=origobj;
	float lrate=1;
	G=CreateDMatrix(&gstack,lnumCeps,lnumCeps);
	int mode;
	if(upNoiseModel){
		jcalMu(A, B, D, K, g, h, h_mu_up, n_mu_up);
		for(mode=0;mode<=1;mode++){
			lrate=1;
			for(k=1;k<=KK;k++){
				switch(mode){
				case 0:
					backoffUpdate(n_mu, n_mu_up, lrate);
					backoffUpdate(h_mu, h_mu_up, lrate);
					break;
				case 1:
					lupCv_pn(pgradient, ngradient, n_cv, &lrate);
					break;
				}
				curobj=0;
				NewHMMScan(hset,&hss);
				while(GoNextMix(&hss,FALSE)){
					mp=hss.mp;
					mixIdx=mp->mIdx;
					MuAcc* ma = (MuAcc *) GetHook(mp->mean);
					VaAcc* va = (VaAcc *) GetHook(mp->cov.inv);
					float occ=ma->occ;
					if(occ<1){
						continue;
					}
					for(i=1;i<=lnumCeps;i++){
						x_mu[i]=lmus_x[mixIdx][i];
						for(j=1;j<=i;j++){
							x_cv[i][j]=x_cv[j][i]=lcvs_x[mixIdx][i][j];
							y_icv[i][j]=y_icv[j][i]=mp->cov.inv[i][j];
						}
					}
					/********************************Jacobian re-calculation @ new expansion point************************************/
					for(i=1;i<=lnumCeps;i++){/*mu=mu_n-mu_x-mu_h*/
						exPoint[i]=n_mu[i]-x_mu[i]-h_mu[i];
					}
					VTS(&gstack, dct, idct, exPoint, G, Gt, F, Ft, g0);
					/*******************************Re-adaptation****************************************/
					for(i=1;i<=lnumCeps;i++){
						y_mu[i]=x_mu[i]+h_mu[i]+g0[i];
					}
					VtsUpVar(G, Gt, F, Ft, GCx, GCxGt, FCn, FCnFt, x_cv, n_cv, y_cv);
					DMatInvert(&gstack, y_cv,y_icv);
					curobj-=occ*(-DCovDet(y_icv));
					for(i=1;i<=lnumCeps;i++){
						y_acc1[i]=ma->mu[i]-occ*y_mu[i];
						for(j=1;j<=i;j++){
							y_acc2[i][j]=y_acc2[j][i]=va->cov.inv[i][j]-ma->mu[i]*y_mu[j]-ma->mu[j]*y_mu[i]+occ*y_mu[i]*y_mu[j];
						}
					}
					float trace=0;
					for(i=1;i<=lnumCeps;i++){
						for(j=1;j<=lnumCeps;j++){
							trace+=y_acc2[i][j]*y_icv[j][i];
						}
					}
					curobj=curobj-trace;
				}
				EndHMMScan(&hss);
				if(curobj>origobj){
					switch(mode){
					case 0:
						DVec2Vec(n_mu,lnm->n_mean);
						DVec2Vec(h_mu,lnm->h_mean);
						printf("\tUpdate lnoise means: curobj=%f, origobj=%f, iteration=%d, lrate=%f\n",curobj,origobj,k,lrate);
						break;
					case 1:
						DMat2Tri(n_cv,lnm->cov.inv);
						printf("\tUpdate lnoise variance: curobj=%f, origobj=%f, iteration=%d, lrate=%f\n",curobj,origobj,k,lrate);
						break;
					}
					origobj=curobj;
					break;
				}
				else{
					switch(mode){
					case 0:
						Vec2DVec(mu_n,n_mu);
						Vec2DVec(mu_h,h_mu);
						break;
					case 1:
						Tri2DMat(cv_n,n_cv);
						break;
					}
					lrate/=2;
				}
			}
		}
	}
	FreeVector(&gstack,memMarker);
	return curobj-origobj0;
}

void NATUpdateLAcousticModel(HMMSet *hset){
	 HMMScanState hss;
	 /*go through acc of each mix*/
	 int i,j;
	 int mixIdx=0;
	 DVector x_mu;
	 DMatrix x_cv;
	 Matrix iGtPG;
	 DMatrix iCov;
	 MixPDF* mp;
	 float lrate=1,temp;
	 int vSize=hset->vecSize;
	 x_mu=CreateDVector(&gstack,vSize);
	 x_cv=CreateDMatrix(&gstack,vSize, vSize);
	 iGtPG=CreateMatrix(&gstack,vSize,vSize);
	 iCov=CreateDMatrix(&gstack,vSize,vSize);
	 NewHMMScan(hset,&hss);
	 while(GoNextMix(&hss,FALSE)){
		 mp=hss.mp;
		 mixIdx=mp->mIdx;
		 MuAcc* ma=(MuAcc *) GetHook(mp->mean);
		 VaAcc* va=(VaAcc *) GetHook(mp->cov.var);
		 if(ma->occ==0)
			 continue;
		 float norm=0;
		 for(i=1;i<=vSize;i++){
			 x_mu[i]=lmus_x[mixIdx][i];
			 for(j=1;j<=vSize;j++){
				 x_cv[i][j]=lcvs_x[mixIdx][i][j];
			 }
			 for(j=1;j<=i;j++){
				 temp=va->lpgradient[i][j]-va->lngradient[i][j];
				 norm+=temp*temp;
			 }
		 }
		 lrate=1;
		 upMu(&gstack, ma->lGtPG, iGtPG, ma->lGtPacc1, x_mu);
		 lupCv_pn(va->lpgradient, va->lngradient, x_cv, &lrate);
		 for(i=1;i<=vSize;i++){
			 mp->mean[i]=x_mu[i];
		 }
		 DMatInvert(&gstack, x_cv,iCov);
		 DMat2Tri(iCov,mp->cov.inv);
		 FixFullGConst(mp,-CovDet(&gstack, mp->cov.inv));
	}
	EndHMMScan(&hss);
	FreeDVector(&gstack,x_mu);
}

void ApplyXForm2VectorX(MemHeap* x, LinXForm *linXForm, Vector mean)
{
   Vector vec, bias;
   int size,b,bsize;
   Matrix A;
   float tmp;
   int i,j;
   int cnt,cnti,cntj;

   /* Check dimensions */
   size = linXForm->vecSize;
   if (size != VectorSize(mean))
      HError(999,"Transform dimension (%d) does not match mean dimension (%d)",
             size,VectorSize(mean));
   vec = CreateVector(x,size);
   CopyVector(mean,vec); ZeroVector(mean);
   /* Transform mean */
   for (b=1,cnti=1,cnt=1;b<=IntVecSize(linXForm->blockSize);b++) {
      bsize = linXForm->blockSize[b];
      A = linXForm->xform[b];
      for (i=1;i<=bsize;i++,cnti++) {
         tmp = 0;
         for (j=1,cntj=cnt;j<=bsize;j++,cntj++)
            tmp += A[i][j] * vec[cntj];
         mean[cnti] = tmp;
      }
      cnt += bsize;
   }
   /* Apply bias if required */
   bias = linXForm->bias;
   if (bias != NULL) {
      for (i=1;i<=size;i++)
         mean[i] += bias[i];
   }
   FreeVector(x,vec);
}

float CalcMiniAuxAndStat(MemHeap* x, NATAcc* acc_head, DMatrix dct, DMatrix idct, DVector* x_mu_work, DVector* x_cv_work, MixPDF* mp, char* accMode){
	int i,j;
	float obj=0;
	DVector n_cv[3],y_mu[3],y_cv[3],y_acc1[3],y_acc2[3],y_icv[3];
	for(i=0;i<3;i++){
		n_cv[i]=CreateDVector(x,numCepCoef);
		y_mu[i]=CreateDVector(x,numCepCoef);
		y_cv[i]=CreateDVector(x,numCepCoef);
		y_acc1[i]=CreateDVector(x,numCepCoef);
		y_acc2[i]=CreateDVector(x,numCepCoef);
		y_icv[i]=CreateDVector(x,numCepCoef);
	}
	DVector exPoint=CreateDVector(x,numCepCoef);
	DVector g0=CreateDVector(x,numCepCoef);
	DMatrix G=CreateDMatrix(x,numCepCoef,numCepCoef);DMatrix Gt=CreateDMatrix(x,numCepCoef,numCepCoef);
	DMatrix F=CreateDMatrix(x,numCepCoef,numCepCoef);DMatrix Ft=CreateDMatrix(x,numCepCoef,numCepCoef);
	DMatrix XX=CreateDMatrix(x,numCepCoef,numCepCoef);DMatrix XXX=CreateDMatrix(x,numCepCoef,numCepCoef);
	DMatrix As=CreateDMatrix(x,numCepCoef,numCepCoef);
	Vector full_mean=CreateVector(&gstack,numCepCoef*3);
	DMatrix G_bk=CreateDMatrix(x,numCepCoef,numCepCoef);
	MuAcc* ma=(MuAcc *) GetHook(mp->mean);
	VaAcc* va=(VaAcc *) GetHook(mp->cov.var);
	/*reset the accumulators*/
	if(strchr(accMode,'m')){
		for(j=0;j<3;j++){
			ZeroTriMat(ma->GtPG[j]);
			ZeroVector(ma->GtPacc1[j]);
		}
	}
	if(strchr(accMode,'v')){
		for(j=0;j<3;j++){
			ZeroTriMat(va->hessian[j]);
			ZeroVector(va->gradient[j]);
		}
	}

	NATAcc* iter=acc_head;
	LinXForm* mllrMean=NULL;
	while((iter=iter->next)!=NULL){
		float occ=iter->occ_k;
		if(occ==0){
			continue;
		}
		if(iter->noise_k->spkr){/*Copy the Jacobian of static parameters*/
			mllrMean=iter->noise_k->spkr->mllrMean;
		}
		/********************************Jacobian re-calculation************************************/
		for(i=1;i<=numCepCoef;i++){/*mu=mu_n-mu_x-mu_h*/
			exPoint[i]=iter->noise_k->n_mean[i]-x_mu_work[0][i]-iter->noise_k->h_mean[i];
			for(j=0;j<3;j++)
				n_cv[j][i]=iter->noise_k->cov.var[numCepCoef*j+i];
		}
		VTS(x, dct, idct, exPoint, G, Gt, F, Ft, g0);
		/*******************************Re-adaptation****************************************/
		for(i=1;i<=numCepCoef;i++){
			y_mu[0][i]=x_mu_work[0][i]+iter->noise_k->h_mean[i]+g0[i];
			DDotProduct(G[i],x_mu_work[1],&y_mu[1][i]);
			DDotProduct(G[i],x_mu_work[2],&y_mu[2][i]);
		}
		for(j=0;j<3;j++){
			dVtsUpVar(G, Gt, F, Ft, x_cv_work[j], n_cv[j], y_cv[j]);
		}
		if(mllrMean){
			for(j=0;j<3;j++){
				for(i=1;i<=numCepCoef;i++){
					full_mean[j*numCepCoef+i]=y_mu[j][i];
				}
			}
			ApplyXForm2VectorX(x, mllrMean,full_mean);/*mllr mean on top of vts adaptation*/
			for(j=0;j<3;j++){
				for(i=1;i<=numCepCoef;i++){
					y_mu[j][i]=full_mean[j*numCepCoef+i];
				}
			}
		}
		for(i=1;i<=numCepCoef;i++){
			for(j=0;j<3;j++){
				y_acc1[j][i]=iter->mu_k[numCepCoef*j+i]-occ*y_mu[j][i];
				y_acc2[j][i]=iter->cov_k.var[numCepCoef*j+i]-2*iter->mu_k[numCepCoef*j+i]*y_mu[j][i]+occ*y_mu[j][i]*y_mu[j][i];
				y_icv[j][i]=1/y_cv[j][i];
				obj-=occ*log(y_cv[j][i]);
				obj-=y_acc2[j][i]/y_cv[j][i];
			}
		}
		if(strchr(accMode,'v')){
			for(j=0;j<3;j++){
				accCv(G, x_cv_work[j], y_icv[j], y_acc2[j], va->gradient[j], va->hessian[j], occ);
			}
		}
		if(strchr(accMode,'m')){
			if(mllrMean){
				CopyDMatrix(G,G_bk);
			}
			for(j=0;j<3;j++){
				if(mllrMean){/*Copy the Jacobian of static parameters, revise the joint Jacobian for mean update*/
					Mat2DMat(mllrMean->xform[1+j],As);
					DMatMult(As,G_bk,G); DTranspose(G,Gt);
				}
				accMu(G, Gt, y_icv[j], XX, XXX, ma->GtPG[j], ma->GtPacc1[j], y_acc1[j], occ);
			}
		}
	}
	FreeDVector(x,n_cv[0]);
	return obj;
}

typedef struct _NATObj{
	int mixIdx;
	float updateCnt;
	float totalSubK;
	float totalFulK;
	MixPDF* mp;
	struct _NATObj *next;
}*NATObj;

void* NATThread(void* arg){
	NATObj noHead=(NATObj)arg;
	int i,j;
	int mixIdx=0;
	int k, bkOffLoop=6;
	float totalSubK=0, totalFulK=0, updateCnt=0;
	float lrate,origobj,curobj;
	int mode,iter;
	int Iters=NatEMiters;
	DMatrix dct, idct;
	Matrix iGtPG;
	Matrix ihessian;
	MemHeap x;/*only header has*/
	char storeName[16];
	sprintf(storeName,"NATStore_%d",noHead->mixIdx);
	CreateHeap(&x,   storeName,    MSTAK, 1, 1.0, 50000,   500000);
	DVector x_mu_work[3],x_cv_work[3],x_mu_up[3],x_cv_up[3];
	for(j=0;j<3;j++){
		x_mu_work[j]=CreateDVector(&x,numCepCoef);
		x_cv_work[j]=CreateDVector(&x,numCepCoef);
		x_mu_up[j]=CreateDVector(&x,numCepCoef);
		x_cv_up[j]=CreateDVector(&x,numCepCoef);
	}
	iGtPG=CreateMatrix(&x,numCepCoef,numCepCoef);
	ihessian=CreateMatrix(&x,numCepCoef,numCepCoef);
	dct=liftDCT(&x, numChans, numCepCoef, cepLifter);
	idct=iLiftDCT(&x, numChans,numCepCoef, cepLifter);
	NATObj noIter=noHead;
	while((noIter=noIter->next)!=NULL){
		MixPDF* mp=noIter->mp;
		mixIdx=noIter->mixIdx;
		MuAcc* ma=(MuAcc *) GetHook(mp->mean);
		VaAcc* va = (VaAcc *) GetHook(mp->cov.var);
		NATAcc* accIter=ma->acc_k;
		while((accIter=accIter->next)!=NULL){/*look up the noise model given its identifier*/
			 NoiseModel *nmIter=nlistHead;
			 while((nmIter=nmIter->next)!=NULL){
				 if(!strcmp(nmIter->id,accIter->noise_k->id)){
					 accIter->noise_k=nmIter;
					 break;
				 }
			 }
		 }
		 for(iter=0;iter<Iters;iter++){
			 for(mode=0;mode<3;mode++){
				 for(i=1;i<=numCepCoef;i++){/*Load current "clean" models*/
					 for(j=0;j<3;j++){
						 x_mu_up[j][i]=x_mu_bk[mixIdx][numCepCoef*j+i];
						 x_cv_up[j][i]=x_cv_bk[mixIdx][numCepCoef*j+i];
					 }
				 }
				 /*Calculate the derivative for NAT update, and the obj value before update*/
				 origobj=CalcMiniAuxAndStat(&x, ma->acc_k, dct, idct, x_mu_up, x_cv_up, mp, mode<2?"m":"v");
				 switch(mode){
				 case 0:/*static mean*/
					 upMu(&x, ma->GtPG[0], iGtPG, ma->GtPacc1[0], x_mu_up[0]);break;
				 case 1:/*dynamic mean*/
					 upMu(&x, ma->GtPG[1], iGtPG, ma->GtPacc1[1], x_mu_up[1]);
					 upMu(&x, ma->GtPG[2], iGtPG, ma->GtPacc1[2], x_mu_up[2]);break;
				 case 2:/*variance*/
					 for(j=0;j<3;j++){
						 upCv(&x, va->hessian[j], ihessian, va->gradient[j], x_cv_up[j]);
					 }break;
				 }
				 lrate=glrate;
				 for(k=0;k<bkOffLoop;k++){
					 for(i=1;i<=numCepCoef;i++){/*Reload the current clean model*/
						 for(j=0;j<3;j++){
							 x_mu_work[j][i]=x_mu_bk[mixIdx][numCepCoef*j+i];
							 x_cv_work[j][i]=x_cv_bk[mixIdx][numCepCoef*j+i];
						 }
					 }
					 /*perform the backoff update, given current model and update model*/
					 switch(mode){
					 case 0:
						 backoffUpdate(x_mu_work[0],x_mu_up[0],lrate);break;
					 case 1:
						 backoffUpdate(x_mu_work[1],x_mu_up[1],lrate);
						 backoffUpdate(x_mu_work[2],x_mu_up[2],lrate);break;
					 case 2:
						 for(j=0;j<3;j++){
							 backoffUpdate(x_cv_work[j],x_cv_up[j],lrate);
						 }break;
					 }
					 /*check the new obj value given the new acoustic model*/
					 curobj=CalcMiniAuxAndStat(&x, ma->acc_k, dct, idct, x_mu_work, x_cv_work, mp, "");
					 if(curobj>origobj+1){
						 updateCnt++;
						 totalSubK+=ma->acc_subK;
						 totalFulK+=ma->acc_K;
						 switch(mode){
						 case 0:
							 for(i=1;i<=numCepCoef;i++){
								 x_mu_bk[mixIdx][i]=mp->mean[i]=x_mu_work[0][i];
							 }
							 printf("Update Gaussian[%d,%d,%d] static mean:origobj=%e, curobj=%e, lrate=%e, subK=%d, fullK=%d\n",mixIdx,iter,mode,origobj,curobj,lrate,ma->acc_subK,ma->acc_K);
							 break;
						 case 1:
							 for(i=1;i<=numCepCoef;i++){
								 for(j=1;j<3;j++){
									 x_mu_bk[mixIdx][numCepCoef*j+i]=mp->mean[numCepCoef*j+i]=x_mu_work[j][i];
								 }
							 }
							 printf("Update Gaussian[%d,%d,%d] dynamic mean:origobj=%e, curobj=%e, lrate=%e, subK=%d, fullK=%d\n",mixIdx,iter,mode,origobj,curobj,lrate,ma->acc_subK,ma->acc_K);
							 break;
						 case 2:
							 for(i=1;i<=numCepCoef;i++){
								 for(j=0;j<3;j++){
									 x_cv_bk[mixIdx][numCepCoef*j+i]=mp->cov.var[numCepCoef*j+i]=x_cv_work[j][i];
								 }
							 }
							 printf("Update Gaussian[%d,%d,%d] variance:origobj=%e, curobj=%e, lrate=%e, subK=%d, fullK=%d\n",mixIdx,iter,mode,origobj,curobj,lrate,ma->acc_subK,ma->acc_K);
							 break;
						 }
						 FixDiagGConst(mp);
						 origobj=curobj;
						 break;
					 }
					 lrate*=(k==0?0.5:lrate);
				 }
			 }
		 }
	}
	noHead->totalSubK=totalSubK;
	noHead->totalFulK=totalFulK;
	noHead->updateCnt=updateCnt;
	/*DeleteHeap(&x);this may cause problem, as it will edit a global variable*/
	return NULL;
}

void NATUpdateAcousticModel(HMMSet *hset, NoiseModel* nlistHead){
	 HMMScanState hss;
	 /*go through acc of each mix*/
	 int i;
	 int mixIdx=0;
	 float totalSubK=0, totalFullK=0, upCount=0;
	 MixPDF* mp;
	 ForceDiagC(hset);
	 int nt;
	 pthread_t *pths=(pthread_t*)New(&gstack,sizeof(pthread_t)*numThreads);
	 NATObj* nos=(NATObj*)New(&gstack,sizeof(NATObj)*numThreads);
	 for(nt=0;nt<numThreads;nt++){
		 nos[nt]=(NATObj)New(&gstack,sizeof(struct _NATObj));
		 nos[nt]->next=NULL;
		 nos[nt]->mixIdx=nt;
	 }
	 NewHMMScan(hset,&hss);
	 nt=0;
	 while(GoNextMix(&hss,FALSE)){
		 mp=hss.mp;
		 mixIdx=mp->mIdx;
		 MuAcc* ma=(MuAcc *) GetHook(mp->mean);
		 if(parMode<0&&mixIdx>0){/*just in case of p<0, default -1, make sure they are clean*/
			 for(i=1;i<=numCepCoef*3;i++){
				 mp->mean[i]=x_mu_bk[mixIdx][i];
				 mp->cov.var[i]=x_cv_bk[mixIdx][i];
			 }
			 FixDiagGConst(mp);
		 }
		 if(mixIdx<=0||ma->acc_subK==0)/*critical bug found*/
			 continue;
		 NATObj no=(NATObj)New(&gstack,sizeof(struct _NATObj));
		 no->mixIdx=mixIdx;
		 no->mp=mp;
		 NATObj temp=nos[nt%numThreads]->next;
		 nos[nt%numThreads]->next=no;
		 no->next=temp;
		 nt++;
	 }
	 EndHMMScan(&hss);
	 for(nt=0;nt<numThreads;nt++){
		 pthread_create(&pths[nt],NULL,NATThread,nos[nt]);
	 }
	 for(nt=0;nt<numThreads;nt++){
		 pthread_join(pths[nt], NULL /* void ** return value could go here */);
	 }
	 for(nt=0;nt<numThreads;nt++){
		 totalSubK+=nos[nt]->totalSubK;
		 totalFullK+=nos[nt]->totalFulK;
		 upCount+=nos[nt]->updateCnt;
	 }
	 printf("#######################################################################################################\n");
	 printf("#########Summary: Updates=%d, Average subK=%.2f, Average fullK=%.2f, Approximate rate=%.4f##########\n",(int)upCount, totalSubK/upCount, totalFullK/upCount, totalSubK/totalFullK);
	 printf("#######################################################################################################\n");

}
/* ----------------------------------------------------------- */
/*                      END:  HERest.c                         */
/* ----------------------------------------------------------- */
