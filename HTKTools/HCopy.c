/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright: Microsoft Corporation                    */
/*          1995-2000 Redmond, Washington USA                  */
/*                    http://www.microsoft.com                 */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*      File: HCopy.c: Copy one Speech File to another         */
/* ----------------------------------------------------------- */

char *hcopy_version = "!HVER!HCopy:   3.4.1 [CUED 12/03/09]";
char *hcopy_vc_id = "$Id: HCopy.c,v 1.1.1.1 2006/10/11 09:54:59 jal58 Exp $";

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HWave.h"
#include "HVQ.h"
#include "HAudio.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HAdapt.h"
/* -------------------------- Trace Flags & Vars ------------------------ */

#define T_TOP     001           /* basic progress reporting */
#define T_KINDS   002           /* report file formats and parm kinds */
#define T_SEGMENT 004           /* output segment label calculations */
#define T_MEM     010           /* debug memory usage */

static int  trace  = 0;         /* Trace level */
typedef struct _TrList *TrPtr;  /* simple linked list for trace info */
typedef struct _TrList {      
   char *str;                   /* output string */
   TrPtr next;                  /* pointer to next in list */
} TrL;
static TrL trList;              /* 1st element in trace linked list */
static TrPtr trStr = &trList;   /* ptr to it */

static int traceWidth = 70;     /* print this many chars before wrapping ln */

static ConfParam *cParm[MAXGLOBS];
static int nParm = 0;            /* total num params */

/* ---------------------- Global Variables ----------------------- */

FileFormat srcFF     = UNDEFF;   /* I/O configuration options */
FileFormat tgtFF     = UNDEFF;
FileFormat srcLabFF  = UNDEFF;
FileFormat tgtLabFF  = UNDEFF;
ParmKind srcPK       = ANON;
ParmKind tgtPK       = ANON;
HTime srcSampRate    = 0.0;
HTime tgtSampRate    = 0.0;
Boolean saveAsVQ = FALSE;

static int swidth0 = 1;

static HTime st=0.0;            /* start of samples to copy */
static HTime en=0.0;            /* end of samples to copy */
static HTime xMargin=0.0;       /* margin to include around extracted labs */
static Boolean stenSet=FALSE;   /* set if either st or en set */
static int labstidx=0;          /* label start index (if set) */
static int labenidx=0;          /* label end index (if set) */
static int curstidx=0;          /* label start index (if set) */
static int curenidx=0;          /* label end index (if set) */
static int labRep=1;            /* repetition of named label */
static int auxLab = 0;          /* auxiliary label to use (0==primary) */
static Boolean chopF = FALSE;   /* set if we should truncate files/trans */

static LabId labName = NULL;    /* name of label to extract (if set) */
static Boolean useMLF=FALSE;    /* set if we are saving to an mlf */
static Boolean labF=FALSE;      /* set if we should  process label files too */
static char *labDir = NULL;     /* label file directory */
static char *outLabDir = NULL;  /* output label dir */
static char *labExt = "lab";    /* label file extension */
static char outLabFile[256];

static int xpwin=1;
static Boolean onlyStatic=FALSE;
static MemHeap hmmStack;   /*For Storage of all dynamic structures created...*/
HMMSet hset;             /* Set of HMMs to be re-estimated */
HMMSet synrhset;
char hmmListFn[MAXFNAMELEN];
char synrListFn[MAXFNAMELEN];
Boolean useGmmSynr=FALSE;
Boolean useTiedStatePos=FALSE;
TreeNode tspSynr=NULL;
HLink* phnSynr=NULL;
Boolean saveCompressed=FALSE;
/*Boolean saveLogPos=FALSE;	/*this is for output features*/

static int numStrs=0;
static int strWidth[5];

static Wave wv;                 /* main waveform; cat all input to this */
static ParmBuf pb;              /* main parmBuf; cat input, xform wv to this */
static Transcription *trans=NULL;/* main labels; cat all input to this */
static Transcription *tr;       /* current transcription */
static char labFile[255];       /* current source of trans */
static HTime off = 0.0;         /* length of files appended so far */

static char * posDir = NULL;	 /* directory to store posterior files */
static char * posExt = "pos";	 /* extension of posterior files */
static XFInfo xfInfo;
static HDFInfo hdfInfo;
static FileFormat dff=UNDEFF;       /* data file format */
static Boolean fMPE=FALSE;
static Boolean cmn=FALSE;
static Boolean cvn=FALSE;
static Boolean forceLog=FALSE;		/*this is for input features*/

/* ---------------- Memory Management ------------------------- */

#define STACKSIZE 100000        /* assume ~100K wave files */
static MemHeap iStack;          /* input stack */
static MemHeap iStack2;
static MemHeap oStack;          /* output stack */
static MemHeap cStack;          /* chop stack */
static MemHeap lStack;          /* label i/o  stack */
static MemHeap tStack;          /* trace list  stack */


static InputXForm* xf=NULL;
static Boolean g2ixform=FALSE;
/*static Vector gmean=NULL;*/
/* ---------------- Process Command Line ------------------------- */

#define MAXTIME 1E13            /* maximum HTime (1E6 secs) for GetChkdFlt */

void ReportUsage(void)
{
   printf("\nUSAGE: HCopy [options] src [ + src ...] tgt ...\n\n");
   printf(" Option                                       Default\n\n");
   printf(" -a i     Use level i labels                  1\n");
   printf(" -e t     End copy at time t                  EOF\n");
   printf(" -i mlf   Save labels to mlf s                null\n");
   printf(" -l dir   Output target label files to dir    current\n");
   printf(" -m t     Set margin of t around x/n segs     0\n");
   printf(" -n i [j] Extract i'th [to j'th] label        off\n");
   printf(" -s t     Start copy at time t                0\n");
   printf(" -t n     Set trace line width to n           70\n");
   printf(" -x s [n] Extract [n'th occ of] label  s      off\n");
   printf(" -E i     Feature expansion to window size i  1\n");
   printf(" -o       Only extract static features        FALSE\n");
   printf(" -W file  MMF for posterior synthesiser       null\n");
   printf(" -p file  phone list for MMF                  null\n");
   printf(" -u       save uncompressed (pos)features     FALSE\n");
   printf(" -g       save log posterior features         FALSE\n");
   printf(" -J dir ext  input xform directory			 null\n");
   printf(" -Y dir ext  posterior feature directory      null\n");
   PrintStdOpts("FGILPOX");
}

/* SetConfParms: set conf parms relevant to this tool */
void SetConfParms(void)
{
   int i;
   Boolean b;
   char buf[MAXSTRLEN];

   nParm = GetConfig("HCOPY", TRUE, cParm, MAXGLOBS);
   if (nParm>0){
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfBool(cParm,nParm,"SAVEASVQ",&b)) saveAsVQ = b;
      if (GetConfInt(cParm,nParm,"NSTREAMS",&i)) swidth0 = i;
      if (GetConfStr(cParm,nParm,"SOURCEFORMAT",buf))
         srcFF = Str2Format(buf);
      if (GetConfStr(cParm,nParm,"TARGETFORMAT",buf))
         tgtFF = Str2Format(buf);
      if (GetConfStr(cParm,nParm,"SOURCEKIND",buf))
         srcPK = Str2ParmKind(buf);
      if (GetConfStr(cParm,nParm,"TARGETKIND",buf)) {
         tgtPK = Str2ParmKind(buf);
         if (tgtPK&HASNULLE) 
            HError(1019, "SetConfParms: incompatible TARGETKIND=%s for coding", buf);
      }
   }
}

/* FixOptions: Check and set config options */
void FixOptions(void)
{
   if (stenSet && (labstidx>0 || labName != NULL))
      HError(1019,"FixOptions: Specify -s/-e or -x but not both");
   if (labstidx>0 && labName != NULL)
      HError(1019,"FixOptions: Specify label index or name but not both");
   if (srcFF == UNDEFF) srcFF = HTK;
   if (tgtFF == UNDEFF) tgtFF = HTK;
   if (tgtPK == ANON) tgtPK = srcPK;
}

/*void loadGmean(char* hmmfn){
        HLink hmmLink;
        MLink macroLink;
        LabId  hmmId  = NULL;
        char base[MAXSTRLEN];
        char path[MAXSTRLEN];
        char ext[MAXSTRLEN];

        /* Load HMM defs *
        if(MakeOneHMM(&hset,BaseOf(hmmfn,base))<SUCCESS)
          HError(2028,"Initialise: MakeOneHMM failed");
        if(LoadHMMSet(&hset,PathOf(hmmfn,path),ExtnOf(hmmfn,ext))<SUCCESS)
          HError(2028,"Initialise: LoadHMMSet failed");

        /* Get a pointer to the physical HMM *
        hmmId = GetLabId(base,FALSE);
        macroLink = FindMacroName(&hset,'h',hmmId);
        if (macroLink==NULL)
          HError(2020,"Initialise: cannot find hmm %s in hset",hmmfn);
        hmmLink = (HLink)macroLink->structure;
        gmean = (((hmmLink->svec+2)->info->pdf+1)->spdf.cpdf+1)->mpdf->mean;
}*/


int main(int argc, char *argv[])
{
   char *s;                     /* next file to process */
   char clean[256];
   char noisy[256];

   char *datafn;
   char posfn[MAXSTRLEN];
   void OpenSpeechFile(char *s);
   HTime MergeParmFile(char **srcs);
   void AppendSpeechFile(char *s);
   void PutTargetFile(char *s);
   HTime ExNoiseSpeechFile(char *noisy, char *clean);

   if(InitShell(argc,argv,hcopy_version,hcopy_vc_id)<SUCCESS)
      HError(1000,"HCopy: InitShell failed");
   InitMem();   InitLabel();
   InitMath();  InitSigP();
   InitWave();  InitAudio();
   InitVQ();    InitModel();
   InitAdapt(&xfInfo);
   InitUtil();
   if(InitParm()<SUCCESS)  
      HError(1000,"HCopy: InitParm failed");

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(0);

   SetConfParms();
   /* initial trace string is null */
   trList.str = NULL;

   CreateHeap(&iStack, "InBuf",   MSTAK, 1, 0.0, STACKSIZE, LONG_MAX);
   CreateHeap(&iStack2, "InBuf",   MSTAK, 1, 0.0, STACKSIZE, LONG_MAX);
   CreateHeap(&oStack, "OutBuf",  MSTAK, 1, 0.0, STACKSIZE, LONG_MAX);
   CreateHeap(&cStack, "ChopBuf", MSTAK, 1, 0.0, STACKSIZE, LONG_MAX);
   CreateHeap(&lStack, "LabBuf",  MSTAK, 1, 0.0, 10000, LONG_MAX);
   CreateHeap(&tStack, "Trace",   MSTAK, 1, 0.0, 100, 200);
   CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   CreateHMMSet(&hset,&hmmStack,TRUE);
   CreateHMMSet(&synrhset,&hmmStack,TRUE);
   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      /*if (strlen(s)!=1)
         HError(1019,"HCopy: Bad switch %s; must be single letter",s);*/
      switch(s[0]){
      case 'a':
         if (NextArg() != INTARG)
            HError(1019,"HCopy: Auxiliary label index expected");
         auxLab = GetChkedInt(1,100000,s) - 1;
         break;
      case 'e':              /* end time in seconds, max 10e5 secs */
         en = GetChkedFlt(-MAXTIME,MAXTIME,s);
         stenSet = TRUE; chopF = TRUE;
         break;
      case 'i':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Output MLF name expected");
         strcpy(outLabFile,GetStrArg());
         if(SaveToMasterfile(outLabFile)<SUCCESS)
            HError(1014,"HCopy: Cannot write to MLF");
         useMLF = TRUE; labF = TRUE; break;
      case 'l':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Target label file directory expected");
         outLabDir = GetStrArg();
         labF = TRUE; break;
      case 'm':
         xMargin = GetChkedFlt(-MAXTIME,MAXTIME,s);
         chopF = TRUE; break;
      case 'n':
         if (NextArg() != INTARG)
            HError(1019,"HCopy: Label index expected");
         labstidx= GetChkedInt(-100000,100000,s);
         if (NextArg() == INTARG)
            labenidx = GetChkedInt(-100000,100000,s);
         chopF = TRUE; break;          
      case 's':      /* start time in seconds */
         st = GetChkedFlt(0,MAXTIME,s);
         stenSet = TRUE; chopF = TRUE; break;
      case 't':
    	  if (!strcmp(s,"tiedStatePos")){
    		  useTiedStatePos=TRUE;
    		  break;
    	  }
         if (NextArg() != INTARG)
            HError(1019,"HCopy: Trace line width expected");
         traceWidth= GetChkedInt(10,100000,s); break;
      case 'x':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Label name expected");
         labName = GetLabId(GetStrArg(),TRUE);
         if (NextArg() == INTARG)
            labRep = GetChkedInt(1,100000,s);
         chopF = TRUE; labF = TRUE; break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Source file format expected");
         if((srcFF = Str2Format(GetStrArg())) == ALIEN)
            HError(-1089,"HCopy: Warning ALIEN src file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Source label File format expected");
         if((srcLabFF = Str2Format(GetStrArg())) == ALIEN)
            HError(-1089,"HCopy: Warning ALIEN Label output file format set");
         labF= TRUE; break;
      case 'I':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: MLF file name expected");
         LoadMasterFile(GetStrArg());
         labF = TRUE; break;
      case 'L':
         if (NextArg()!=STRINGARG)
            HError(1019,"HCopy: Label file directory expected");
         labDir = GetStrArg();
         labF = TRUE; break;
      case 'P':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Label File format expected");
         if((tgtLabFF = Str2Format(GetStrArg())) == ALIEN)
            HError(-1089,"HCopy: Warning ALIEN Label file format set");
         labF = TRUE; break;
      case 'O':
         if (NextArg() != STRINGARG)
            HError(1019,"HCopy: Target file format expected");
         if((tgtFF = Str2Format(GetStrArg())) == ALIEN)
            HError(-1089,"HCopy: Warning ALIEN target file format set");
         break;
      case 'T':
         trace = GetChkedInt(0,16,s); break;
      case 'X':
         if (NextArg()!=STRINGARG)
            HError(1019,"HCopy: Label file extension expected");
         labExt = GetStrArg();
         labF = TRUE; break;     
      case 'E':
    	  if (NextArg()!=INTARG)
    		HError(1019,"HCopy: Feature expansion window size expected");
    	  xpwin=GetIntArg();break;
      case 'o':
    	  if (!strcmp(s,"os"))
    		  onlyStatic=TRUE;
    	  else
    		  HError(1019,"HCopy: Unknown switch %s",s);
    	  break;
      case 'M':
    	  if (NextArg()!=INTARG){
    		  HError(1019, "HCopy: number of feature streams expected");
    	  }
    	  numStrs=strWidth[0]=GetIntArg();
    	  int i;
    	  for(i=1;i<=numStrs;i++){
        	  if (NextArg()!=INTARG){
        		  HError(1019, "HCopy: stream[%d] width expected",i);
        	  }
        	  strWidth[i]=GetIntArg();
    	  }
    	  break;
      case 'W':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCopy: -W: Posterior synthesiser model file expected");
    	 AddMMF(&synrhset,GetStrArg());
    	 useGmmSynr=TRUE;
    	 break;
      case 'p':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCopy: -p: Posterior synthesiser list file expected");
    	 strcpy(synrListFn,GetStrArg());
    	 break;
      case 'Z':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCopy: -Z: MMF(fMPE) definition file expected");
    	 AddMMF(&hset,GetStrArg());
    	 break;
      case 'z':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HCopy: -z: model(fMPE) phone list file expected");
    	 strcpy(hmmListFn,GetStrArg());
    	 break;
      case 'f':
    	  if(!strcmp(s,"fl")){
    		  forceLog=TRUE;
    	  }
    	  if(!strcmp(s,"fMPE")){
    		  fMPE=TRUE;
    	  }
    	  break;
      case 'j':
    	  if(NextArg()!=STRINGARG){
    		  HError(2319,"HCopy: input transform file path expected");
    	  }
    	  char *fname=GetStrArg();
    	  char macroname[32];
    	  NameOf(fname, macroname);
    	  xf=LoadInputXForm(NULL, macroname, fname);
    	  break;
      case 'g':/*force global xform to input xform*/
          /*if(!strcmp(s,"gm")){
             loadGmean(GetStrArg());
             break;
          }*/
    	  xf = (InputXForm *)New(hset.hmem,sizeof(InputXForm));
    	  xfInfo.useInXForm = TRUE;
    	  g2ixform=TRUE;
    	  break;
      case 'h':
    	  if (NextArg()!=STRINGARG)
    		  HError(1,"Speaker name pattern expected");
    	  xfInfo.outSpkrPat = GetStrArg();
    	  break;
      case 'J':
         if (NextArg()!=STRINGARG)
            HError(2319,"HCopy: input transform directory expected");
         AddInXFormDir(&hset,GetStrArg());
         if (NextArg()==STRINGARG) {
            if (xfInfo.inXFormExt == NULL)
               xfInfo.inXFormExt = GetStrArg();
            else
               HError(2319,"HCopy: only one input transform extension may be specified");
         }
         if (NextArg() != SWITCHARG)
        	 HError(2319,"HCopy: cannot have -J as the last option");
         break;
      case 'Y':
    	 if (NextArg()!=STRINGARG)
    		HError(2319,"HCopy: -Y: Posterior directory expected");
    	 posDir = GetStrArg();
    	 if (NextArg()==STRINGARG)
    	     posExt = GetStrArg();
    	 break;
      case 'c':
    	  if(!strcmp(s,"crc")){
    		  saveCompressed=TRUE;
    		  break;
    	  }
    	  if(NextArg()==STRINGARG){
    		  char* ops=GetStrArg();
    		  cmn=strchr(ops,'m')!=NULL?TRUE:FALSE;
    		  cvn=strchr(ops,'v')!=NULL?TRUE:FALSE;
    		  break;
    	  }
    	  break;
      default:
         HError(1019,"HCopy: Unknown switch %s",s);break;
      }
   }
   if (NumArgs() == 1)
      HError(1019,"HCopy: Target file or + operator expected");
   FixOptions();
   if(useGmmSynr){/*load gmmcluster */
	   if(useTiedStatePos)
		   tspSynr=BuildClusteredGMMTree(&synrhset, synrListFn);
	   else
		   phnSynr=LoadPosSynth(&synrhset,synrListFn);
   }
   if(fMPE||g2ixform){
	   if(MakeHMMSet( &hset, hmmListFn, NULL )<SUCCESS) /* make the containers for each HMM */
	      HError(2321,"Initialise: MakeHMMSet failed");
	   if(LoadHMMSet( &hset,NULL,NULL)<SUCCESS) /* load the parameters into the containers from the specified macro model file */
	      HError(2321,"Initialise: LoadHMMSet failed");
	   SetParmHMMSet(&hset);
	   InitHDFInfo(&hdfInfo, &hset);
   }
   if(g2ixform){
	   if (xfInfo.inSpkrPat == NULL) xfInfo.inSpkrPat = xfInfo.outSpkrPat;
	   if (xfInfo.paSpkrPat == NULL) xfInfo.paSpkrPat = xfInfo.outSpkrPat;
   }
   char **srcs=NULL;
   int i;
   if(numStrs>1){
	   srcs=(char**)New(&gstack,sizeof(char*)*(numStrs+1));
	   for(i=1;i<=numStrs;i++){
		   srcs[i]=(char*)New(&gstack,sizeof(char)*256);
	   }
   }
   while (NumArgs()>0) { /* process group S1 + S2 + ... TGT */
      off = 0.0;
      if (NextArg()!=STRINGARG)
         HError(1019,"HCopy: Source file name expected");    
      s = GetStrArg();
      datafn=s;
      if(g2ixform){
    	  UpdateSpkrStats(&hset,&xfInfo, datafn);
    	  xf->xform=xfInfo.inXForm->xformSet->numXForms==1?xfInfo.inXForm->xformSet->xforms[1]:NULL;
      }
      strcpy(noisy,s);
      if(posDir!=NULL){
    	  MakeFN(datafn,posDir,posExt,posfn);
    	  printf("Reading posterior from file: %s\n", posfn); fflush(stdout);
      }
      if(phnSynr!=NULL){
    	  printf("Synthesise posterior for file: %s\n", datafn); fflush(stdout);
      }
      if(fMPE){
    	  if(phnSynr&&((xfInfo.inXForm&&xfInfo.inXForm->xformSet->numXForms>1)||hset.xf)){
    		  HError(1,"Main: Current posterior synthesiser only support one linear transformation, **non input xform version**");
    	  }
    	  if(useTiedStatePos){
    		  LoadHDF4fMPE2(&hdfInfo, dff, datafn, tspSynr, hset.TVWR);
    	  }
    	  else{
    		  LoadHDF4fMPE(&hdfInfo, dff, posfn, datafn, &synrhset, phnSynr, hset.TVWR, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL);
    	  }
    	  /*LoadHDF4fMPE(&hdfInfo, dff, posfn, datafn, &synrhset, phnSynr, FALSE, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL);*/
      }

      if(numStrs>1){
    	  i=0;
    	  strcpy(srcs[++i],s);
          if (NextArg()!=STRINGARG)
        	  HError(1019,"HCopy: Target file or + operator expected");
          s = GetStrArg();


          while (strcmp(s,"+") == 0) {     /* Append + S2 + S3 ... */
             if (NextArg()!=STRINGARG)
                HError(1019,"HCopy: Append file name expected");
             s = GetStrArg();
             strcpy(srcs[++i],s);
             if (NextArg()!=STRINGARG)
                HError(1019,"HCopy: Target file or + operator expected");
             s = GetStrArg();
          }
          assert(i==numStrs);
          MergeParmFile(srcs);
      }
      else{
          OpenSpeechFile(s);               /* Load initial file  S1 */
          /*skip the target file*/
          if (NextArg()!=STRINGARG)
        	  HError(1019,"HCopy: Target file or + operator expected");
          s = GetStrArg();

          while (strcmp(s,"+") == 0) {     /* Append + S2 + S3 ... */
             if (NextArg()!=STRINGARG)
                HError(1019,"HCopy: Append file name expected");
             s = GetStrArg();
             strcpy(clean,s);
             AppendSpeechFile(s);
             if (NextArg()!=STRINGARG)
                HError(1019,"HCopy: Target file or + operator expected");
             s = GetStrArg();
          }
      }
      PutTargetFile(s);
      if(trace & T_MEM) PrintAllHeapStats();
      if(trans != NULL){
         trans = NULL;
         ResetHeap(&lStack);
      }
      ResetHeap(&iStack);
      ResetHeap(&oStack);
      if(chopF) ResetHeap(&cStack);
   }
   if(useMLF) CloseMLFSaveFile();
   if (NumArgs() != 0) HError(-1019,"HCopy: Unused args ignored");
   Exit(0);
   return (0);          /* never reached -- make compiler happy */
}

/* ----------------- Trace linked list handling ------------------------ */

/* AppendTrace: insert a string to trStr for basic tracing */
void AppendTrace(char *str)
{
   TrPtr tmp = trStr;

   /* Seek to end of list */
   while (tmp->str != NULL) tmp = tmp->next;
   tmp->str =  CopyString(&tStack, str);
   tmp->next = (TrPtr)New(&tStack,sizeof(trList));
   tmp->next->str = NULL;
   tmp->next->next = NULL;
}

/* PrintTrace: Print trace linked list */
void PrintTrace(void)
{
   int linelen = 0;
   TrPtr tmp = trStr;

   /* print all entries in list */
   while (tmp->next != NULL){
      printf("%s ",tmp->str);
      linelen += strlen(tmp->str) + 1;
      if (linelen > traceWidth && tmp->next->next!=NULL){
         printf("\n    ");  /* wrap line where appropriate */
         linelen = 0;
      }
      tmp = tmp->next;
   }
   if(linelen > 0) printf("\n");
}

/* ------------------- Utility Routines ------------------------ */

/* ClampStEn: set/clamp  st/en times */
void ClampStEn(HTime length, HTime *st, HTime *en)
{  
   *st -= xMargin;
   if (*st < 0) *st = 0;
   
   if( *en > 0.0 ){             /* Absolute time */
      *en += xMargin;
   }
   else if( *en < 0.0 ){        /* Relative to end */
      *en = length + *en + xMargin;
      if (*en >= length) *en = length;
      if (*en < *st) *en = *st;
      if (*st > *en) *st = *en;
   }
   else                         /* default to eof */
      *en = length - xMargin;

   /* Now clamp */
   if (*en >= length) *en = length;
   if (*en < *st) *en = *st;
   if (*st > *en) *st = *en;
}

/* ----------------- Label Manipulation ------------------------ */

/* FixLabIdxs: -ve idxs count from end, so set +ve and check */
void FixLabIdxs(int nlabs)
{
   if (labstidx<0) curstidx = nlabs + 1 + labstidx;
   else curstidx = labstidx;
   if (labenidx<0) curenidx = nlabs + 1 + labenidx;
   else curenidx = labenidx;
   if (curstidx < 0 || curstidx > nlabs)
      HError(1030,"FixLabIdxs: label start index [%d] out of range",curstidx);
   if (curenidx < curstidx || curstidx > nlabs)
      HError(1030,"FixLabIdxs: label end index  [%d] out of range",curenidx);
}

/* SetLabSeg: Set st and en for label (sequence) */
void SetLabSeg(Transcription *tr)
{
   LabList *ll = tr->head;    /* use first lab list */
   LLink p,q;
   
   if (tr->numLists > 1)
      HError(-1031,"SetLabSeg: label lists 2 to %d will be ignored",
             tr->numLists);
   if (labName != NULL) {  /* extract labName */
      if (auxLab==0) {
         p = GetCase(ll,labName,labRep);
         st = p->start; en = p->end;
      } else {
         p = GetAuxCase(ll,labName,labRep,auxLab);
         st = p->start; en = AuxLabEndTime(p,auxLab);
      }        
   } else {                /* extract labstidx to labenidx */
      if (auxLab==0){
         FixLabIdxs(CountLabs(ll));
         p = GetLabN(ll,curstidx);
         q = GetLabN(ll,curenidx);
         st = p->start; en = q->end;
      }else{
         FixLabIdxs(CountAuxLabs(ll,auxLab));
         p = GetAuxLabN(ll,curstidx,auxLab);
         q = GetAuxLabN(ll,curenidx,auxLab);
         st = p->start; en = AuxLabEndTime(q,auxLab);
      }
   }  
   if(trace & T_SEGMENT)
      printf("Extracting %8.0f to %8.0f\n",st,en);
}

/* LoadTransLabs: Load transcription from file */
Transcription *LoadTransLabs(char *src)
{
   Transcription *t;
   
   MakeFN(src,labDir,labExt,labFile);
   if(trace & T_SEGMENT)
      printf("Loading label file %s\n",labFile);
   t = LOpen(&lStack,labFile,srcLabFF);
   if(chopF && ! stenSet) SetLabSeg(t);
   return t;
}

/* SaveLabs: save trans t to label file corresponding to tgt */
void SaveLabs(char *tgt, Transcription *t)
{
   MakeFN(tgt,outLabDir,labExt,labFile);
   if(trace & T_SEGMENT)
      printf("Saving label file %s\n",labFile);
   if (CountLabs(trans->head) == 0) 
      HError(-1031,"SaveLabs: No labels in transcription %s", labFile);
   if(LSave(labFile,t,tgtLabFF)<SUCCESS)
      HError(1014,"SaveLabs: Could not save label file %s", labFile);
}

/* AppendLabs: append label file corresponding to src to trans,
   len is time length of this file; accumulate to find offset for
   concatenated files */
void AppendLabs(Transcription *t, HTime len)
{
   LabList *ll,*transll;
   LLink p,q;
   int maxAux;

   if(trace & T_SEGMENT)
      printf("Adding labels, len: %.0f off: %.0f\n",len,off);
   if(trans == NULL) 
      trans = CopyTranscription(&lStack, t);
   else
      for (ll = t->head,transll = trans->head; ll != NULL;
           ll = ll->next,transll = transll->next){
         if (transll == NULL)
            HError(1031,"AppendLabs: lablist has no target to append to");
         maxAux = ll->maxAuxLab;
         if (maxAux > transll->maxAuxLab){
            HError(-1031,"AppendLabs: truncating num aux labs from %d down to %d",
                   maxAux, transll->maxAuxLab);
            maxAux = transll->maxAuxLab;
         }
         for (p=ll->head->succ; p->succ!= NULL; p=p->succ){
            q = AddLabel(&lStack,transll,
                         p->labid,p->start + off,p->end + off,p->score);
            if (maxAux>0)
               AddAuxLab(q,maxAux,p->auxLab,p->auxScore);
         }
      }
   /* accumulate length of this file for total offset */
   off += len;
}

/* ChopLabs: Chop trans around stime to etime. end = 0 means to end of file */
void ChopLabs(Transcription *t, HTime start, HTime end)
{
   LabList *ll;
   LLink p;
   HTime del = tgtSampRate;
   
   /* assume ClampStEn has already been called, so start and end are OK */
   if(trace & T_SEGMENT)
      printf("ChopLab: extracting labels %.0f to %.0f\n",start,end);
   for (ll = t->head; ll != NULL; ll = ll->next){
      for (p=ll->head->succ; p->succ!= NULL; p=p->succ){
         if((p->start < start-del) || ((end > 0.0 ) && (p->end > end+del))) {
            DeleteLabel(p);
         }
         else {
            p->start -= start;
            p->end -= start;
         }
      }
   }
}

/* ----------------------- Wave File Handling ------------------------ */

/* ChopWave: return wave chopped to st and end. end = 0 means all */
Wave ChopWave(Wave srcW, HTime start, HTime end, HTime sampRate)
{
   Wave tgtW;
   HTime length;                /* HTime length of file */
   long stSamp, endSamp, nSamps;
   short *data;
   
   data = GetWaveDirect(srcW,&nSamps);
   length = nSamps * sampRate;
   if(start >= length)
      HError(1030,"ChopWave: Source too short to get data from %.0f",start); 
   ClampStEn(length,&start,&end);
   if(trace & T_SEGMENT)
      printf("ChopWave: Extracting data %.0f to %.0f\n",start,end);
   stSamp = (long) (start/sampRate);
   endSamp = (long) (end/sampRate);
   nSamps = endSamp - stSamp;
   if(nSamps <= 0)
      HError(1030,"ChopWave: Truncation options result in zero-length file"); 
   tgtW = OpenWaveOutput(&cStack,&sampRate,nSamps);
   PutWaveSample(tgtW,nSamps,data + stSamp);
   CloseWaveInput(srcW);
   if(chopF && labF) ChopLabs(tr,start,end);
   return(tgtW);
}

/* IsWave: check config parms to see if target is a waveform */
Boolean IsWave(char *srcFile)
{
   FILE *f;
   long nSamp,sampP, hdrS;
   short sampS,kind;
   Boolean isPipe,bSwap,isWave;
   
   isWave = tgtPK == WAVEFORM;
   if (tgtPK == ANON){
      if ((srcFF == HTK || srcFF == ESIG) && srcFile != NULL){
         if ((f=FOpen(srcFile,WaveFilter,&isPipe)) == NULL)
            HError(1011,"IsWave: cannot open File %s",srcFile);
         switch (srcFF) {
         case HTK:
            if (!ReadHTKHeader(f,&nSamp,&sampP,&sampS,&kind,&bSwap))
               HError(1013, "IsWave: cannot read HTK Header in File %s",
                      srcFile);
            break;
         case ESIG:
            if (!ReadEsignalHeader(f, &nSamp, &sampP, &sampS,
                                   &kind, &bSwap, &hdrS, isPipe))
               HError(1013, "IsWave: cannot read Esignal Header in File %s",
                      srcFile);             
            break;
         }
         isWave = kind == WAVEFORM;
         FClose(f,isPipe);
      } else
         isWave = TRUE;
   }
   return isWave;
}

/* OpenWaveFile: open source wave file and extract portion if indicated */
HTime OpenWaveFile(char *src)
{
   Wave w, cw;
   long nSamps;
   short *data;

   if((w = OpenWaveInput(&iStack,src,srcFF,0,0,&srcSampRate))==NULL)
      HError(1013,"OpenWaveFile: OpenWaveInput failed");
   srcPK = WAVEFORM;
   tgtSampRate = srcSampRate;
   cw = (chopF)?ChopWave(w,st,en,srcSampRate) : w;
   data = GetWaveDirect(cw,&nSamps);
   wv = OpenWaveOutput(&oStack, &srcSampRate, nSamps);
   PutWaveSample(wv,nSamps,data);
   CloseWaveInput(cw);
   return(nSamps*srcSampRate);
}

/* AppendWave: append the src file to global wave wv */
HTime AppendWave(char *src)
{
   Wave w, cw;
   HTime period=0.0;
   long nSamps;
   short *data;

   if((w = OpenWaveInput(&iStack,src, srcFF, 0, 0, &period))==NULL)
      HError(1013,"AppendWave: OpenWaveInput failed");
   if(trace & T_KINDS )
      printf("Appending file %s format: %s [WAVEFORM]\n",src,
             Format2Str(WaveFormat(w)));   
   if(period != srcSampRate)
      HError(1032,"AppendWave: Input file %s has inconsistent sampling rate",src);
   cw = (chopF)? ChopWave(w,st,en,srcSampRate) : w;
   data = GetWaveDirect(cw,&nSamps);
   PutWaveSample(wv,nSamps,data);
   CloseWaveInput(cw);
   return(nSamps*period);
}

/* ----------------------- Parm File Handling ------------------------ */

/* ChopParm: return parm chopped to st and end. end = 0 means all */
ParmBuf ChopParm(ParmBuf b, HTime start, HTime end, HTime sampRate)
{  
   int stObs, endObs, nObs, i;
   HTime length;
   short swidth[SMAX];
   Boolean eSep;
   ParmBuf cb;
   Observation o;
   BufferInfo info;

   length =  ObsInBuffer(b) * sampRate;
   ClampStEn(length,&start,&end);
   if(start >= length)
      HError(1030,"ChopParm: Src file too short to get data from %.0f",start);
   if(trace & T_SEGMENT)
      printf("ChopParm: Extracting segment %.0f to %.0f\n",start,end);
   stObs = (int) (start/sampRate);
   endObs = (int) (end/sampRate);
   nObs = endObs -stObs;
   if(nObs <= 0)
      HError(1030,"ChopParm: Truncation options result in zero-length file");
   GetBufferInfo(b,&info);
   ZeroStreamWidths(swidth0,swidth);
   SetStreamWidths(tgtPK,info.tgtVecSize,swidth,&eSep);
   o = MakeObservation(&cStack, swidth, info.tgtPK, saveAsVQ, eSep);
   if (saveAsVQ){
      if (info.tgtPK&HASNULLE){
         info.tgtPK=DISCRETE+HASNULLE;
      }else{
         info.tgtPK=DISCRETE;
      }
   }
   cb =  EmptyBuffer(&cStack, nObs, o, info);
   for (i=stObs; i < endObs; i++){
      ReadAsTable(b, i, &o);
      AddToBuffer(cb, o);
   }
   CloseBuffer(b);
   if(chopF && labF) ChopLabs(tr,start,end);
   return(cb);
}

/* AppendParm: append the src file to current Buffer pb. Return appended len */
HTime AppendParm(char *src)
{  
   int i;
   char bf1[MAXSTRLEN]; 
   char bf2[MAXSTRLEN]; 
   short swidth[SMAX];
   Boolean eSep;
   ParmBuf b, cb;
   Observation o;
   BufferInfo info;

   if((b =  OpenBuffer(&iStack,src,0,srcFF,TRI_UNDEF,TRI_UNDEF))==NULL)
      HError(1050,"AppendParm: Config parameters invalid");
   GetBufferInfo(b,&info);
   if(trace & T_KINDS ){
      printf("Appending file %s format: %s [%s]->[%s]\n",src,
             Format2Str(info.srcFF), ParmKind2Str(info.srcPK,bf1),
             ParmKind2Str(info.tgtPK,bf2));
   }
   if  (tgtSampRate != info.tgtSampRate)
      HError(1032,"AppendParm: Input file %s has inconsistent sample rate",src);
   if ( BaseParmKind(tgtPK) != BaseParmKind(info.tgtPK))
      HError(1032,"AppendParm: Input file %s has inconsistent tgt format",src);
   cb = (chopF)?ChopParm(b,st,en,info.tgtSampRate) : b;
   ZeroStreamWidths(swidth0,swidth);
   SetStreamWidths(info.tgtPK,info.tgtVecSize,swidth,&eSep);
   o = MakeObservation(&iStack, swidth, info.tgtPK, saveAsVQ, eSep);
   for (i=0; i < ObsInBuffer(cb); i++){
      ReadAsTable(cb, i, &o);
      AddToBuffer(pb, o);
   }
   CloseBuffer(cb);
   return(i*info.tgtSampRate);
}

/* OpenParmFile: open source parm file and return length */
HTime OpenParmFile(char *src)
{
   int i,j,k,halfwin,numObs,idx;
   ParmBuf b, cb;
   short swidth[SMAX],swidthx[SMAX];
   Boolean eSep;
   Observation o,ox;
   BufferInfo info,infox;
   int tgtvsize,obsvsize,lobsvsize;
   LogDouble tempsum,x,temp;
   Boolean offsetObs;
   Vector tempObs;
   Vector xpObs;
   if((b =  OpenBuffer(&iStack,src,0,srcFF,TRI_UNDEF,TRI_UNDEF))==NULL)
      HError(1050,"OpenParmFile: Config parameters invalid");
   GetBufferInfo(b,&info);
   srcSampRate = info.srcSampRate;
   tgtSampRate = info.tgtSampRate;
   srcPK = info.srcPK;
   cb = chopF?ChopParm(b,st,en,info.tgtSampRate):b;
   ZeroStreamWidths(swidth0,swidth);
   SetStreamWidths(info.tgtPK,info.tgtVecSize,swidth,&eSep);
   o = MakeObservation(&oStack, swidth, info.tgtPK, saveAsVQ, eSep);
   if (saveAsVQ){
      if (info.tgtPK&HASNULLE){
         info.tgtPK=DISCRETE+HASNULLE;
      }else{
         info.tgtPK=DISCRETE;
      }
   }
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
   xpObs=CreateVector(&gstack,lobsvsize);
   if(xpwin>1&&!useGmmSynr){/*feature expansion*/
	   GetBufferInfo(b,&infox);
	   tgtvsize=lobsvsize;
	   infox.tgtVecSize=tgtvsize;
	   infox.tgtPK=tgtPK;
	   infox.saveWithCRC=infox.saveCompressed=saveCompressed;/*force compress*/
	   ZeroStreamWidths(swidth0,swidthx);
	   SetStreamWidths(infox.tgtPK,infox.tgtVecSize,swidthx,&eSep);
	   ox =  MakeObservation(&oStack, swidthx, infox.tgtPK, saveAsVQ, eSep);
	   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb), ox, infox);
	   halfwin=xpwin/2;
	   offsetObs=(xpwin-1)%2;/*e.g. win=9->no offset, win=8->offset*/
	   tempObs=CreateVector(&gstack,info.tgtVecSize);
	   numObs=ObsInBuffer(cb);
	   for(i=0; i < numObs; i++){
		   if(offsetObs){
			   ReadAsTable(cb, i, &o);
			   CopyVector(o.fv[swidth0],tempObs);
		   }
		   for(j=-halfwin;j<=halfwin;j++){
			   if(offsetObs&&j==0){/*no need to read again the current frame*/
				   continue;
			   }
			   k=i+j;
			   k=k<0?0:(k>=numObs?numObs-1:k);
			   ReadAsTable(cb, k, &o);
			   for(idx=1;idx<=obsvsize;idx++){/*tgtvsize refers to the part of features we care about*/
				   if(offsetObs){
					   xpObs[idx+obsvsize*((j>0?j-1:j)+halfwin)]=o.fv[swidth0][idx];/*-tempObs[idx];/*offset the current frame*/
				   }
				   else{
					   xpObs[idx+obsvsize*(j+halfwin)]=o.fv[swidth0][idx];
				   }
			   }
		   }
		   CopyVector(xpObs,ox.fv[swidth0]);
	       AddToBuffer(pb, ox);
	   }
   }
   else if(useGmmSynr){/*synthesise posterior using gmm*/
	   GetBufferInfo(b,&infox);
	   tgtvsize=synrhset.numPhyHMM;/*switch to posterior dimension*/
	   infox.tgtVecSize=tgtvsize;/*number of phones*/
	   infox.tgtPK=tgtPK;
	   /*compression may cause problem for rare occurred posterior*/
	   /*check out function CalcCompress from HParm, the way to calculate A[nx]*/
	   infox.saveWithCRC=infox.saveCompressed=saveCompressed;
	   ZeroStreamWidths(swidth0,swidthx);
	   SetStreamWidths(infox.tgtPK,infox.tgtVecSize,swidthx,&eSep);
	   ox =  MakeObservation(&oStack, swidthx, infox.tgtPK, saveAsVQ, eSep);
	   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb), ox, infox);
	   halfwin=xpwin/2;
	   offsetObs=(xpwin-1)%2;/*e.g. win=9->no offset, win=8->offset*/
	   tempObs=CreateVector(&gstack,info.tgtVecSize);
	   numObs=ObsInBuffer(cb);
	   for(i=0; i < numObs; i++){
		   if(offsetObs){
			   ReadAsTable(cb, i, &o);
			   CopyVector(o.fv[swidth0],tempObs);
		   }
		   for(j=-halfwin;j<=halfwin;j++){
			   if(offsetObs&&j==0){/*no need to read again the current frame*/
				   continue;
			   }
			   k=i+j;
			   k=k<0?0:(k>=numObs?numObs-1:k);
			   ReadAsTable(cb, k, &o);
			   for(idx=1;idx<=obsvsize;idx++){/*tgtvsize refers to the part of features we care about*/
				   if(offsetObs){
					   xpObs[idx+obsvsize*((j>0?j-1:j)+halfwin)]=o.fv[swidth0][idx];/*-tempObs[idx];/*offset the current frame*/
				   }
				   else{
					   xpObs[idx+obsvsize*(j+halfwin)]=o.fv[swidth0][idx];
				   }
			   }
		   }
		   tempsum=LZERO;
		   for(j=1;j<=tgtvsize;j++){
			   x=StateOutprob(&synrhset,(phnSynr[j]->svec+2)->info->pdf+1,xpObs,NULL);
			   ox.fv[swidth0][j]=x;
			   tempsum=LAdd(tempsum,x);
		   }
		   float maxLike=-1E10;
		   for(j=1;j<=tgtvsize;j++){
			   temp=ox.fv[swidth0][j]-tempsum;
			   ox.fv[swidth0][j]=forceLog?temp:exp(temp);
			   if(ox.fv[swidth0][j]>maxLike){
				   maxLike=ox.fv[swidth0][j];
			   }
		   }
	       AddToBuffer(pb, ox);
	   }
	   FreeVector(&gstack,tempObs);
   }
   else if(xf){/*Input feature transformation*/
	   GetBufferInfo(b,&infox);
	   Matrix mat=xf->xform?xf->xform->xform[1]:NULL;
	   if(mat==NULL){
		   infox.tgtVecSize=tgtvsize=obsvsize;
	   }
	   else{
		   assert(obsvsize==NumCols(mat));
	   	   infox.tgtVecSize=tgtvsize=NumRows(mat);
	   }
	   infox.tgtPK=tgtPK;/*force MFCC*/
	   /*infox.tgtPK=info.tgtPK;
	   /*compression may cause problem for rare occurred posterior*/
	   /*check out function CalcCompress from HParm, the way to calculate A[nx]*/
	   infox.saveWithCRC=infox.saveCompressed=saveCompressed;
	   ZeroStreamWidths(swidth0,swidthx);
	   SetStreamWidths(infox.tgtPK,infox.tgtVecSize,swidthx,&eSep);
	   ox =  MakeObservation(&oStack, swidthx, infox.tgtPK, saveAsVQ, eSep);
	   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb), ox, infox);

	   tempObs=CreateVector(&gstack,info.tgtVecSize);
	   numObs=ObsInBuffer(cb);
	   for(i=0; i < numObs; i++){
		   ReadAsTable(cb, i, &o);
		   for(j=1;j<=tgtvsize;j++){
			   ox.fv[swidth0][j]=forceLog?(o.fv[swidth0][j]>0?log(o.fv[swidth0][j]):LSMALLP):o.fv[swidth0][j];
		   }
		   if(xf->xform){
			   ApplyXForm2Vector(xf->xform, ox.fv[swidth0]);
		   }
		   /*for(j=1;j<=tgtvsize;j++){
			   if(mat==NULL){
				   ox.fv[swidth0][j]=o.fv[swidth0][j]-(gmean?gmean[j]:0);
			   }
			   else{
				   ox.fv[swidth0][j]=0;
				   for(k=1;k<=obsvsize;k++){
					   temp=forceLog?(o.fv[swidth0][k]>0?log(o.fv[swidth0][k]):LZERO):o.fv[swidth0][k];
					   ox.fv[swidth0][j]+=(temp-(gmean?gmean[k]:0))*mat[j][k];
				   }
			   }
		   }*/
		   AddToBuffer(pb, ox);
	   }
	   FreeVector(&gstack,tempObs);
   }
   else if(cmn||cvn){/*cepstral mean variance normalization*/
	   info.saveWithCRC=info.saveCompressed=saveCompressed;
	   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb), o, info);
	   Vector mean=CreateVector(&gstack, info.tgtVecSize);
	   Vector var=CreateVector(&gstack, info.tgtVecSize);
	   ZeroVector(mean);
	   ZeroVector(var);
	   int T=ObsInBuffer(cb);
	   for(i=0; i < T; i++){
		   ReadAsTable(cb, i, &o);
		   for(j=1;j<=info.tgtVecSize;j++){
			   mean[j]+=o.fv[swidth0][j];
			   var[j]+=o.fv[swidth0][j]*o.fv[swidth0][j];
		   }
	   }
	   for(j=1;j<=info.tgtVecSize;j++){
		   mean[j]/=T;
		   var[j]=var[j]/T-mean[j]*mean[j];
	   }
	   for(i=0; i < T; i++){
		   ReadAsTable(cb, i, &o);
		   for(j=1;j<=info.tgtVecSize;j++){
			   if(cmn)/*only static cepstral parameters are subtracted*/
				   o.fv[swidth0][j]-=mean[j];
			   if(cvn)/*both static and dynamic parameters are scaled*/
				   o.fv[swidth0][j]/=(var[j]>0?sqrt(var[j]):1);/*in case of zero variance*/
		   }
		   AddToBuffer(pb, o);
	   }
   }
   else{
	   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb), o, info);
	   for(i=0; i < ObsInBuffer(cb); i++){
		   ReadAsTable(cb, i, &o);
		   if(fMPE){/*fMPE feature transformation*/
			   CopyVector(ApplyFprojection(&hdfInfo, i+1),o.fv[swidth0]);
		   }
		   if(forceLog){
			   for(j=1;j<=info.tgtVecSize;j++){
				   o.fv[swidth0][j]=o.fv[swidth0][j]>0?log(o.fv[swidth0][j]):LZERO;
			   }
		   }
		   AddToBuffer(pb, o);
	   }
   }
   FreeVector(&gstack,xpObs);
   CloseBuffer(cb);
   if( info.nSamples > 0 )
      return(info.nSamples*srcSampRate);
   else
      return(ObsInBuffer(pb)*info.tgtSampRate);
}


HTime MergeParmFile(char **srcs)
{
   int i,j,numObs,s;
   ParmBuf b[5], cb[5];
   short swidth[SMAX],swidthx[SMAX];
   char pk[MAXSTRLEN];
   Boolean eSep;
   Observation o[5],ox;
   BufferInfo info[5],infox;
   int tgtvsize=0;
   swidthx[0]=numStrs;
   Vector memMark=CreateVector(&iStack,1);
   for(s=1;s<=numStrs;s++){
	   sprintf(pk,"HPARM%d",s);
	   SetNewConfig(pk);
	   if((b[s] =  OpenBuffer(&iStack,srcs[s],0,srcFF,TRI_UNDEF,TRI_UNDEF))==NULL)
	      HError(1050,"OpenParmFile: Config parameters invalid");
	   GetBufferInfo(b[s],&info[s]);
	   srcSampRate = info[s].srcSampRate;
	   tgtSampRate = info[s].tgtSampRate;
	   srcPK = info[s].srcPK;
	   cb[s] = chopF?ChopParm(b[s],st,en,info[s].tgtSampRate):b[s];
	   ZeroStreamWidths(swidth0,swidth);
	   SetStreamWidths(info[s].tgtPK,info[s].tgtVecSize,swidth,&eSep);
	   o[s] = MakeObservation(&oStack, swidth, info[s].tgtPK, saveAsVQ, eSep);
	   GetBufferInfo(b[s],&infox);
	   swidthx[s]=strWidth[s];
	   tgtvsize+=strWidth[s];
   }
   infox.tgtVecSize=tgtvsize;
   infox.tgtPK=tgtPK;/*force MFCC*/
   /*compression may cause problem for rare occurred posterior*/
   /*check out function CalcCompress from HParm, the way to calculate A[nx]*/
   infox.saveWithCRC=infox.saveCompressed=saveCompressed;
   ox =  MakeObservation(&oStack, swidthx, infox.tgtPK, saveAsVQ, eSep);
   pb =  EmptyBuffer(&oStack, ObsInBuffer(cb[1]), ox, infox);
   numObs=ObsInBuffer(cb[1]);
   for(i=0; i < numObs; i++){
	   for(s=1;s<=numStrs;s++){
		   ReadAsTable(cb[s], i, &o[s]);
		   for(j=1;j<=swidthx[s];j++){
			   ox.fv[s][j]=o[s].fv[swidth0][j];
		   }
	   }
	   AddToBuffer(pb, ox);
   }
   for(s=1;s<=numStrs;s++)
	   CloseBufferWithoutClean(cb[s]);
   FreeVector(&iStack,memMark);
   if( info[1].nSamples > 0 )
      return(info[1].nSamples*srcSampRate);
   else
      return(ObsInBuffer(pb)*info[1].tgtSampRate);
}


/* --------------------- Speech File Handling ---------------------- */
void OpenSpaf(char *spafn){/*Particularly designed for PCA*/
	short swidthx[SMAX];
	Boolean eSep;
	Observation ox;
	BufferInfo infox;
	int ivSize,ovSize;
	int i,j,k,dnn_row,dnn_col;
	char ch;
	char dnn_buf[32];
	float value;
	FILE* dnn_f;
	if((dnn_f=fopen(spafn,"r"))==NULL){
		HError(1,"Cannot open file: %s",spafn);
	}
	fscanf(dnn_f,"%d %d\n",&dnn_row,&dnn_col);/*row for time, col for dim*/
	ivSize=dnn_col;
	Matrix mat=xf->xform?xf->xform->xform[1]:NULL;
	Vector bias=xf->xform?xf->xform->bias:NULL;
	if(mat==NULL){
		infox.tgtVecSize=ovSize=ivSize;
	}
	else{
		infox.tgtVecSize=ovSize=NumRows(mat);
		assert(ivSize==NumCols(mat));
	}
	infox.tgtSampRate=100000;
	infox.tgtPK=tgtPK;/*force MFCC*/
	/*compression may cause problem for rare occurred posterior*/
	/*check out function CalcCompress from HParm, the way to calculate A[nx]*/
	infox.saveWithCRC=infox.saveCompressed=saveCompressed;
	ZeroStreamWidths(swidth0,swidthx);
	SetStreamWidths(infox.tgtPK,infox.tgtVecSize,swidthx,&eSep);
	Vector tempV=CreateVector(&oStack,ivSize);
	ox =  MakeObservation(&oStack, swidthx, infox.tgtPK, saveAsVQ, eSep);
	pb =  EmptyBuffer(&oStack, dnn_row, ox, infox);

	for(i=1;i<=dnn_row;i++){
		for(j=1;j<=ivSize;j++){/*set all small posterior to be floor rather than zero to avoid LZERO*/
			tempV[j]=forceLog?LSMALLP:SMALLP;
		}
		j=0;
		while((ch=getc(dnn_f))!=EOF){
			if(ch==' '||ch=='\n'){
				dnn_buf[j++]='\0';
				if(j>1){
					sscanf(dnn_buf,"%d:%e",&j,&value);
					tempV[j]=forceLog?(value>0?log(value):LSMALLP):value;
				}
				j=0;
				if(ch=='\n')
					break;
			}
			else{
				dnn_buf[j++]=ch;
			}
		}
		for(j=1;j<=ovSize;j++){
			if(mat==NULL){
				ox.fv[swidth0][j]=tempV[j];
			}
			else{
				ox.fv[swidth0][j]=0;
				for(k=1;k<=ivSize;k++){
					ox.fv[swidth0][j]+=tempV[k]*mat[j][k];
				}
				ox.fv[swidth0][j]+=(bias?bias[j]:0);/*check the definition of bias*/
			}
		}
		AddToBuffer(pb, ox);
	}
	fclose(dnn_f);
}

/* OpenSpeechFile: open waveform or parm file */
void OpenSpeechFile(char *s)
{
   HTime len;
   char ext[12];
   char buf[MAXSTRLEN];
   ExtnOf(s,ext);
   if(!strcmp(ext,"spa")){
	   OpenSpaf(s);
   }
   else{
	   if (labF) tr = LoadTransLabs(s);
	   if(IsWave(s))
		   len = OpenWaveFile(s);
	   else
		   len = OpenParmFile(s);
	   if(labF) AppendLabs(tr,len);
   }
   if (trace & T_TOP) AppendTrace(s);
   if (tgtPK == ANON) tgtPK = srcPK;      
   if(trace & T_KINDS){
      printf("Source file format: %s [%s]\n",
             Format2Str(srcFF), ParmKind2Str(srcPK,buf));
      printf("Target file format: %s [%s]\n",
             Format2Str(tgtFF), ParmKind2Str(tgtPK,buf));
      printf("Source rate: %.0f Target rate: %.0f \n",
             srcSampRate,tgtSampRate);
   }
}

/* AppendSpeechFile: open waveform or parm file */
void AppendSpeechFile(char *s)
{
   HTime len;
   
   if (labF) tr = LoadTransLabs(s);
   if(tgtPK == WAVEFORM)
      len = AppendWave(s);
   else
      len = AppendParm(s);
   if(labF){
      AppendLabs(tr,len);
   }
   if (trace & T_TOP) { 
      AppendTrace("+"); AppendTrace(s);
   }
}

/* PutTargetFile: close and store waveform or parm file */
void PutTargetFile(char *s)
{
	if(tgtPK == WAVEFORM) {
		if(CloseWaveOutput(wv,tgtFF,s)<SUCCESS)
			HError(1014,"PutTargetFile: Could not save waveform file %s", s);
	}
	else {
		if(SaveBuffer(pb,s,tgtFF)<SUCCESS)
			HError(1014,"PutTargetFile: Could not save parm file %s", s );
		CloseBuffer(pb);
	}
   if (trace & T_TOP){
      AppendTrace("->"); AppendTrace(s);
      PrintTrace();     
      ResetHeap(&tStack);
      trList.str = NULL;
   }
   if(trans != NULL)
      SaveLabs(s,trans);
}

/* ----------------------------------------------------------- */
/*                      END:  HCopy.c                          */
/* ----------------------------------------------------------- */
