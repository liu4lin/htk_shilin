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
/* ----------------------------------------------------------- */
/*         Copyright:                                          */
/*                                                             */
/*         2002-2004  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*        File: HMMIRest.c: Discriminatave (MMI/MPE)  training */
/*     Using Frame Discrimination.                             */
/* ----------------------------------------------------------- */

#define EXITSTATUS 0 /*2 for gprof.*/

char *hmmirest_version = "!HVER!HMMIRest:   3.4.1 [CUED 12/03/09]";
char *hmmirest_vc_id = "$Id: HMMIRest.c,v 1.1.1.1 2006/10/11 09:55:01 jal58 Exp $";

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
#include "HFB.h"
#include "HDict.h"
#include "HNet.h"
#include "HLM.h"
#include "HLat.h"
#include "HArc.h"
#include "HFBLat.h"
#include "HExactMPE.h"


#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define ABS(a) ((a)>0?(a):-(a))
#define FINITE(x) (!isnan(x) && x<1.0e+30 && x>-1.0e+30)


/* Trace Flags */
#define T_TOP   0001    /* Top level tracing */
#define T_TIM   0002    /* Output timings */

/* possible values of updateMode */
#define UPMODE_DUMP 1
#define UPMODE_UPDATE 2
#define UPMODE_BOTH 3

#define PJMMODE_UNIFORM 11
#define PJMMODE_GSMPL   12
#define PJMMODE_MSMPL	13


static char *denLatDir [MAXLATS];  /*MMI lattices.*/
int nDenLats = 0;
static char *numLatDir [MAXLATS];  /*Numerator-alignment lattices.*/
int nNumLats = 0;
static char *latExt    = "lat";


static char *latFileMask = NULL;
static char *LatMask_Numerator   = NULL;
static char *LatMask_Denominator = NULL;

static char numLatSubDirPat[MAXSTRLEN] = "\0";  /* path mask of numerator lattices */
static char denLatSubDirPat[MAXSTRLEN] = "\0";  /* path mask of denominator lattices */

/* Global Settings */

static char * hmmDir = NULL;     /* directory to look for hmm def files */
static char * hmmExt = NULL;     /* hmm def file extension */
static char * newDir = NULL;     /* directory to store new hmm def files */
static char * newExt = NULL;     /* extension of new reestimated hmm files */
static char * posDir[5];	 /* directory to store posterior files */
static int numposdirs=0;
static char * posExt = "pos";	 /* extension of posterior files */
static char * statFN=0;            /* stats file, if any */
static char *dictFN=0;             /* not needed at the moment. */
static float minVar  = 0.0;      /* minimum variance (diagonal only). 
                                    This gives a constant floor to all dimensions' vars. */
      
static int updateMode = UPMODE_UPDATE; /* dump summed accs, update models or do both? */


static float mixWeightFloor=MINMIX*2; /* Floor for mixture weights */
/* There *was* an initialisation of 0.0.  However, this interacts very badly
   with the discriminative technique, which really relies on them being
   taken back up (? I wrote this earlier, dont know if it is still true) 
   I have changed it further from 1.01 to 2, which I had usually put in
   the HTE file.  This enables less stuff in the HTE without me having to
   do experiments to determine the best value.  
*/

static int nSnt      = 0;        /* num sentences from current speaker */

static UPDSet uFlags = 0;   /* update flags */
static UPDSet uFlagsAccs = 0;   /* used in storing accs. */

/*static UPDSet uFlags = UPMEANS|UPVARS|UPTRANS|UPMIXES|UPREGS;   /* update flags */
/*static UPDSet uFlagsAccs = UPMEANS|UPVARS|UPTRANS|UPMIXES|UPREGS;   /* used in storing accs. */
static UPDSet uFlagsMLE = 0; /*which we only update with MLE, ignoring the MMI parameters.*/


static int parMode   = -1;       /* enable one of the parallel modes */
/* i.e.  0 for reestimation, 1,2,3... for accumulating .acc files. */
 

static Boolean stats = FALSE;    /* enable statistics reports */
static char * mmfFn  = NULL;     /* output MMF file, if any */
static long int trace     = 1;        /* Trace level */

static Boolean saveBinary = FALSE;  /* save output HMMs in binary  */

static Boolean ML_MODE = FALSE;     /* when only one set of accs are supplied. */
static Boolean MPE = FALSE;         /* when we are doing MPE/MWE. */
static Boolean MMIPrior = FALSE;    /* use MMI prior as I-smoothing prior */
static float MMITauI = 0.0;        /* I-smoothing tau for MMI prior in MPE training */

static float E = 2.0;               /* constant used in BW updatel */
static float DFACTOR = 2.0;         /* not a config, a constant used in estimating
                                       constant D based on a multiple of the min D for vars > 0 */
static float CWeights = 1.0;
static float CTrans = 1.0;


static Boolean MPEStoreML=FALSE;   /*  Set TRUE if we need to accumulate ML stats while doing MPE. */
static Boolean THREEACCS = FALSE;  /*  Set TRUE if 3 sets of accs need to be stored (for MPE, or possibly for MMI-MAP) */
static int NumAccs;                /*  Set in Initialise() to  1 or 2 or 3. */

static float ISmoothTau = 0.0;        /* I-smoothing: a h-crit-like thing.  Set to 100 for MMI or 50 for MPE, or 25 for MWE. */

/* If unset and ISmoothTau is set, try setting these to 1/10 the occupancy for ISmoothTau. */
static float ISmoothTauTrans = 0.0;
static float ISmoothTauWeights = 0.0;
static float ISmoothTauRWeights = 0.0;
static Boolean ISmoothTauTransSet=FALSE;
static Boolean ISmoothTauWeightsSet=FALSE;
static Boolean ISmoothTauRWeightsSet=FALSE;

/* R.E. priors: */
static float PriorTau = 0.0;         /* tau value [e.g. 10,25] for use in discriminative MAP with -Hprior option. */
static float PriorK = 0.0;           /* e.g. use 1 if disc. trained prior model is better than ML trained new model.  */

static float PriorTauWeights = 0.0;
static float PriorTauTrans = 0.0;


static Boolean twoDataFiles = FALSE; /* For training using two data files, probably never needed & code not tested. */

static float MinOcc = 10;              /* Minimum numerator (ML) occupancy for a Gaussian to be updated */
static float MinOccTrans = 10;         /* Minimum numerator (ML) occupancy for a transition row */
static float MinOccWeights = 10;       /* Minimum numerator (ML) occupancy for a set of weights. */
/* MinOccTrans & MinOccWeights are really only applicable in the way
   intended if ISmoothTauTrans>0 and ISmoothTauTrans>0.  To get around this,
   the function UpdateWeightsOrTrans applies an absolute minimum occupancy
   of 0.0.  This avoids crashes from there being no occupancy. */
         

static float hcrit = 1.0; /* Scale on denominator (MMI) part. Not really useful.  */

static float varFloorPercent = 0;
static float varSmooth = 0;

static Boolean useLLF = FALSE;          /* use directory based LLF files instead of individual lattices */
/* Global non-config variables */



static FileFormat dff=UNDEFF;       /* data file format */
static ConfParam *cParm[MAXGLOBS];   /* configuration parameters */
static int nParm = 0;                /* total num params */

static int maxM = 0;             /* max mixtures in any model */
static int vSize;                /* input vector size */
static int S;                    /* number of data streams */   /*! Equals 1 or error![?] */
static HSetKind hsKind;          /* kind of the HMM system */ /*!Must be PLAINHS | SHAREDHS (|TIEDHS?) */

static int L;                        /* number of logical HMM's */
static int P;                        /* number of physical HMM's */
static HMMSet hset;                  /* Set of HMMs to be re-estimated */
static HMMSet hset_prior;            /* Usually uninitialised, except for MMI-MAP/MPE-MAP */
static HMMSet hset_fullc;
Boolean hset_prior_initialised = FALSE;
static char *hset_prior_dir = NULL;

static FBLatInfo  fbInfo;            /* Structure for discriminative forward-backward. */

static int totalConst=0,nonFlooredConst=0; /*TODO: print.*/


static long totalT=0;                             /* total number of frames in training data */
static LogDouble totalPr1=0,totalPr2=0,totalPr3=0;              /* total log prob upto current utterance, totalPr3 for MMI den in MPE with MMI Prior*/


static Vector vFloor[SMAX];          /* variance floor - default is all zero */

static MemHeap accStack;           /* accumulated statistics */
static MemHeap hmmStack;           /* HMM defs and related structures */
static MemHeap transStack;         /* Transcriptions... comment in original HERest says transformations, but it's wrong.*/

static MemHeap latStack;           /* Lattices. */


Vocab vocab;

/* MPE: */
static int TotalNWords=0;
static double TotalCorr=0;

/* information about transforms */
static XFInfo xfInfo;

/* static prior */
static Boolean STATICPRIOR = FALSE;

/*------------------------fMPE-------------------------*/

static int pnorm=2;
static float RE = 2.0;               /* constant used to control learning rate in fMPE */
static Boolean fMPE=FALSE;
static int targetDim=39;
static int inputDimM=120;
static int inputDimO=39;
static int inputDimP=120;
static int octxwin=0;
static int pctxwin=1;
static int numRegions=0;
static char* hldaname=NULL;
static float Tau=100;

/*----------------------TVWR------------------------*/
static short posize=0;
static short poswidth[SMAX];
static char *dnn_label_list=NULL;

void UpdateWeight(int s, StreamElem *ste);
/* ------------------ Process Command Line -------------------------- */
   
/* SetConfParms: set conf parms relevant to HMMIRest  */

void SetConfParms(void)
{
   int i;
   Boolean b;
   double f;
   char buf[MAXSTRLEN];
  
   nParm = GetConfig("HMMIREST", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfFlt(cParm,nParm,"VARFLOORPERCENTILE",&f)) varFloorPercent = f;
      if (GetConfFlt(cParm,nParm,"VARSMOOTH",&f)) varSmooth = f;
      if (GetConfFlt(cParm,nParm,"C",&f)){ CWeights = CTrans = f; }
      if (GetConfFlt(cParm,nParm,"CW",&f)) CWeights = f;
      if (GetConfFlt(cParm,nParm,"CT",&f)) CTrans = f;
      if (GetConfFlt(cParm,nParm,"MINOCC",&f)){ MinOcc = f; }
      if (GetConfFlt(cParm,nParm,"MINOCCTRANS",&f)){ MinOccTrans = f; }
      if (GetConfFlt(cParm,nParm,"MINOCCWEIGHTS",&f)){ MinOccWeights = f; }
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&b)) saveBinary = b;
      if (GetConfFlt(cParm,nParm,"E",&f)) E = f;
      /*------------fMPE-------------*/
      if (GetConfFlt(cParm,nParm,"RE",&f)) RE = f;
      if (GetConfBool(cParm,nParm,"FMPE",&b)){
    	  MPE=fMPE=b;
    	  THREEACCS=MPE&&MPEStoreML;
      }
      if (GetConfStr(cParm,nParm,"TRANSKIND",buf)){
    	  if (!strcmp (buf, "FPROJ"))  fMPE = TRUE;
    	  MPE=fMPE;
    	  THREEACCS=MPE&&MPEStoreML;
      }
      if (GetConfInt(cParm,nParm,"TARGETDIM",&i)){
    	  targetDim=i;
      }
      if (GetConfInt(cParm,nParm,"INPUTDIMO",&i)){
    	  inputDimO=i;
      }
      if (GetConfInt(cParm,nParm,"INPUTDIMP",&i)){
    	  inputDimP=i;
      }
      if (GetConfInt(cParm,nParm,"OCTXWIN",&i)){
    	  octxwin=i;
      }
      if (GetConfInt(cParm,nParm,"PCTXWIN",&i)){
    	  pctxwin=i;
      }
      if (GetConfInt(cParm,nParm,"NUMREGIONS",&i)){
    	  numRegions=i;
      }
      if (GetConfStr(cParm,nParm,"HLDAINIT",buf)){
    	  hldaname= (char*)malloc(strlen(buf)+1);
    	  strcpy(hldaname, buf);
      }
      /*--------------------------*/
      if (GetConfFlt(cParm,nParm,"DFACTOROCC",&f)) E = f; /*Back-compat. */
      if (GetConfFlt(cParm,nParm,"HCRIT",&f)) hcrit = f;
      if (GetConfBool(cParm,nParm,"MPE",&b)){  MPE = b; THREEACCS=MPE&&MPEStoreML; }
      if (GetConfBool(cParm,nParm,"MWE",&b)){  MPE = b; THREEACCS=MPE&&MPEStoreML; } /* "MWE" has identical effects here but differs in HFBLat.c */
      if (GetConfBool(cParm,nParm,"MEE",&b)){  MPE = b; THREEACCS=MPE&&MPEStoreML; } /* Back-compat. */
      if (GetConfBool(cParm,nParm,"MLE",&b)){  ML_MODE=TRUE; THREEACCS=FALSE; 
      uFlagsMLE =  UPMEANS|UPVARS|UPTRANS|UPMIXES; }
      if (GetConfBool(cParm,nParm,"MMIPRIOR",&b)){  MMIPrior = b;}
      if (GetConfFlt(cParm,nParm,"MMITAUI",&f)){ MMITauI = f;}
      if (GetConfFlt(cParm,nParm,"ISMOOTHTAU",&f)){ ISmoothTau = f; MPEStoreML=TRUE; THREEACCS=MPE&&MPEStoreML; }
      if (GetConfFlt(cParm,nParm,"ICRITOCC",  &f)){ ISmoothTau = f; MPEStoreML=TRUE; THREEACCS=MPE&&MPEStoreML; } /*back-compat. */

      if (GetConfFlt(cParm,nParm,"ISMOOTHTAUT",&f)){ ISmoothTauTrans = f; ISmoothTauTransSet=TRUE; }
      if (GetConfFlt(cParm,nParm,"ISMOOTHTAUW",&f)){ ISmoothTauWeights = f; ISmoothTauWeightsSet=TRUE; }
      if (GetConfFlt(cParm,nParm,"ISMOOTHTAUR",&f)){ ISmoothTauRWeights = f; ISmoothTauRWeightsSet=TRUE; MPEStoreML=TRUE; THREEACCS=MPE&&MPEStoreML;}


      if (GetConfFlt(cParm,nParm,"PRIORTAU",&f)){ PriorTau = f; } 
      if (GetConfFlt(cParm,nParm,"PRIORTAUW",&f)){ PriorTauWeights = f; } 
      if (GetConfFlt(cParm,nParm,"PRIORTAUT",&f)){ PriorTauTrans = f; } 
      if (GetConfFlt(cParm,nParm,"PRIORK",&f)){ PriorK = f; }
      if (GetConfBool(cParm,nParm,"STATICPRIOR",&b)) {
         STATICPRIOR=b;
         if (STATICPRIOR) PriorK = 1.0; else PriorK = 0.0;
      }

      if (GetConfFlt(cParm,nParm,"MIXWEIGHTFLOOR",&f)){ mixWeightFloor = MINMIX * f; }


      if (GetConfStr(cParm,nParm,"LATFILEMASK",buf)) {
         latFileMask = (char*)malloc(strlen(buf)+1); 
         strcpy(latFileMask, buf);
      }

      if (GetConfStr(cParm,nParm,"LATMASKNUM",buf)) {
         LatMask_Numerator = (char*)malloc(strlen(buf)+1); 
         strcpy(LatMask_Numerator,buf);
      }
      if (GetConfStr(cParm,nParm,"LATMASKDEN",buf)) {
         LatMask_Denominator = (char*)malloc(strlen(buf)+1); 
         strcpy(LatMask_Denominator,buf);
      }
      if (GetConfStr(cParm,nParm,"INXFORMMASK",buf)) {
         xfInfo.inSpkrPat = (char*)malloc(strlen(buf)+1);
         strcpy(xfInfo.inSpkrPat,buf);
      }
      if (GetConfStr(cParm,nParm,"PAXFORMMASK",buf)) {
         xfInfo.paSpkrPat = (char*)malloc(strlen(buf)+1);
         strcpy(xfInfo.paSpkrPat,buf);
      }

      if (GetConfBool(cParm,nParm,"USELLF",&b))  useLLF = b;

      if (GetConfStr(cParm,nParm,"UPDATEMODE",buf)) {
         if (!strcmp (buf, "DUMP")) updateMode = UPMODE_DUMP;
         else if (!strcmp (buf, "UPDATE")) updateMode = UPMODE_UPDATE;
         else if (!strcmp (buf, "BOTH")) updateMode = UPMODE_BOTH;
         else HError(2319, "Unknown UPDATEMODE specified (must be DUMP, UPDATE or BOTH)");
      }
      /*--------------------TVWR-----------------------*/
      if (GetConfInt(cParm,nParm,"POSIZE",&i)){
    	  posize=i;
    	  poswidth[1]=i;
    	  poswidth[0]=1;
      }
      if (GetConfInt(cParm,nParm,"POSIZE1",&i)){
    	  posize+=i;
    	  poswidth[1]=i;
    	  poswidth[0]=1;
      }
      if (GetConfInt(cParm,nParm,"POSIZE2",&i)){
    	  posize+=i;
    	  poswidth[2]=i;
    	  poswidth[0]=2;
      }
      if (GetConfInt(cParm,nParm,"POSIZE3",&i)){
    	  posize+=i;
    	  poswidth[3]=i;
    	  poswidth[0]=3;
      }
      if (GetConfInt(cParm,nParm,"POSIZE4",&i)){
    	  posize+=i;
    	  poswidth[4]=i;
    	  poswidth[0]=4;
      }
   }

   if (MPE && uFlagsMLE) HError(1, "Can't combine MPE with ML update of some parameters (code could be simply added).");
   if (MMIPrior && !THREEACCS) HError(999, "MMI Prior must be used in MPE update (THREEACCS).");
   if ((STATICPRIOR && PriorK==0.0) || (!STATICPRIOR && PriorK==1.0)) HError(999, "Specify either PRIORK or STATICPRIOR (PRIORK overwrites value given by STATICPRIOR).");
   /*   if(ISmoothTau && !ISmoothTauTransSet){ ISmoothTauTrans = 10; printf("Smoothing transitions with tau=%f since ISMOOTHTAUT not set\n",ISmoothTauTrans); }
        if(ISmoothTau && !ISmoothTauWeightsSet){ ISmoothTauTrans = 10; printf("Smoothing weights with tau=%f since ISMOOTHTAUW not set\n",ISmoothTauWeights); } */
}

void ReportUsage(void)
{
   printf("\nUSAGE: HMMIRest [options] hmmList dataFiles...\n\n");
   printf(" Option                                   Default\n\n");
   printf(" -a      Use an input linear transform        off\n");
   printf(" -d s    dir to find hmm definitions       current\n");
   printf(" -D f    dictionary file.                  none   \n");
   printf(" -g      MLE updates only.                   \n");
   printf(" -h s    set output speaker name pattern   *.%%%%%%\n");
   printf("         to s, optionally set input and parent patterns\n");
   printf(" -l N    set max sentences (useful for debug) all\n"); 
   printf(" -m N    set min examples needed per model   3\n");
   printf(" -o s    extension for new hmm files        as src\n");
   printf(" -p N    set parallel mode to N             off\n");
   printf(" -q dir   Directory for numerator lats.        [needed. May use >1 -q option]\n");
   printf(" -qp s   Subdirectory pattern string for numerator lats.    none\n");
   printf(" -r dir   Directory for denominator lats.       [needed. May use >1 -r option]\n");
   printf(" -rp s   Subdirectory pattern string for denominator lats.  none\n");
   printf(" -s s    print statistics to file s         off\n"); 
   printf(" -three  (in re-est phase) expect 3 accs, for MPE  off\n"); 
   printf(" -twodatafiles  two data files for single-pass retrain  off\n"); 
   printf(" -u tmvwr update t)rans m)eans v)ars w)ghts r)regs...\n");
   printf(" -umle tmvw do ML updates for specified parameters.\n");
   printf(" -v f    set minimum variance to f          0.0\n");
   printf(" -w f    set mix weight floor to f*MINMIX   0.0\n");
   printf(" -x s    extension for hmm files            none\n");
   /* printf(" -y f    dictionary file.                   none\n");  not needed now. */
   /* printf(" -z   combine all accs into one acc, HDR0.acc.{1,2,3}  off \n"); */
   printf(" -Q      Lattice file extension              lat \n");
   printf(" -Y f    posterior file list, htk format required for posterior file\n");
   printf(" -Z s    posterior file extension            pos \n");
   PrintStdOpts("BEFHIMSTJX"); /*K,G,L removed*/
   printf(" Note: doesn't work if means and variances are shared independently.\n");
   printf("\n\n");
}

void SetuFlags(UPDSet *uFlags)
{
   char *s;
   
   s=GetStrArg();
   *uFlags=0;        
   while (*s != '\0')
      switch (*s++) {
      case 't': (*uFlags)+=UPTRANS; break;
      case 'm': (*uFlags)+=UPMEANS; break;
      case 'v': (*uFlags)+=UPVARS;  break;
      case 'w': (*uFlags)+=UPMIXES; break;
      case 'a': (*uFlags)+=UPXFORM; break;
      case 'p': (*uFlags)+=UPFPROJ; break;
      case 's': (*uFlags)+=UPSEMIT;break;
      case 'r': (*uFlags)+=UPRWGTS; break;
      }
}


void PrintCriteria(){
   printf("\nMMI criterion per frame is: %f (%f - %f)\n", (totalPr1-totalPr2)/totalT,totalPr1/totalT,totalPr2/totalT);
   if(MPE) printf("\nMPE/MWE criterion is: %f ( %f / %d )\n", TotalCorr/TotalNWords, TotalCorr, TotalNWords);
   if(!MPE || MPEStoreML) printf("\nML criterion per frame is: %f (%f/%d)\n", totalPr1/totalT,totalPr1, totalT);
}


int main(int argc, char *argv[]) 
{
   char datafn1[MAXSTRLEN], *datafn, *datafn2, *s,  latfn[MAXSTRLEN], datafn_lat[MAXFNAMELEN], posfn[MAXSTRLEN];
   Lattice *denLats[MAXLATS], *numLats[MAXLATS]; int latn;
   char *statfn=NULL;
   int maxSnt=0;
   void Initialise(char *hmmListFn);
   void UpdateModels(void);
   void UpdateFproj(void);
   void StatReport(void);
   void LoadMPEMLEStat(char *statfn);

   char newFn[MAXSTRLEN];
   FILE *f;
   Source src;
   float x; int i;
   HMMSet synrhset;
   char synrListFn[MAXFNAMELEN];
   Boolean useGmmSynr=FALSE;
   Boolean useTiedStatePos=FALSE;
   HLink* phnSynr=NULL;
   TreeNode tspSynr=NULL;

   int pos_stream;
   posDir[1]=NULL;
   
   setbuf(stdout, NULL); /*unbuffered output.*/
   InitShell(argc,argv,hmmirest_version,hmmirest_vc_id);
   InitMem();
   InitMath();
   InitSigP();
   InitAudio();
   InitWave();
   InitVQ();
   InitModel();
   if(InitParm()<SUCCESS)    HError(1/*was 2300*/,"HMMIRest: InitParm failed");
   InitLabel();
   InitTrain();
   InitUtil();
   InitFBLat();
   InitExactMPE();
   InitArc();
   InitDict();
   InitLat();
   InitNet();
   InitAdapt(&xfInfo); 

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(EXITSTATUS);
   SetConfParms(); /*!Set up config variables.*/
   CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   CreateHMMSet(&hset,&hmmStack,TRUE);
   CreateHMMSet(&synrhset,&hmmStack,TRUE);
   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      switch(s[0]){
      case 'd':
        if(!strcmp(s, "d")){
          if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: HMM definition directory expected"); 
          hmmDir = GetStrArg(); 
        }
        else if (!strcmp(s, "dll")){
  		  dnn_label_list=GetStrArg();
  		  break;
        }
        else if(!strcmp(s, "dprior")){
          if(!hset_prior_initialised){
            hset_prior_initialised=TRUE;
            CreateHMMSet(&hset_prior,&hmmStack,TRUE);
          }
          hset_prior_dir = GetStrArg();
        } else HError(1, "Unknown option %s",s);
        break;   

      case 'g': ML_MODE=TRUE; THREEACCS=FALSE;/*This is the option used during re-estimation when we are only using one set of accs.*/
         uFlagsMLE =  UPMEANS|UPVARS|UPTRANS|UPMIXES; /*TODO, check if necessary. */
         break; 
      case 'l':
         maxSnt = GetChkedInt(0,1000,s); break;
      case 'o':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: HMM file extension expected");
         newExt = GetStrArg(); break;
      case 'p':
         parMode = GetChkedInt(0,500,s);    break; 
      case 'q':
         if (!strcmp(s, "q"))
            numLatDir[nNumLats++] = GetStrArg();
         else if (!strcmp(s, "qp")){
            strcpy(numLatSubDirPat, GetStrArg());
            if (strchr(numLatSubDirPat,'%')==NULL)
               HError(2319,"HMMIRest: Numerator path mask invalid");
         }
         else HError(1, "Unknown option %s",s);
         break;
      case 'r':
         if (!strcmp(s, "r"))
            denLatDir[nDenLats++] = GetStrArg();
         else if (!strcmp(s, "rp")){
            strcpy(denLatSubDirPat,GetStrArg());
            if (strchr(denLatSubDirPat,'%')==NULL)
               HError(2319,"HMMIRest: Denominator path mask invalid");
         }
         else HError(1, "Unknown option %s",s);
         break; 
      case 's':
         stats = TRUE;
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: Stats file name expected");
         statFN = GetStrArg(); break;
      case 't':
    	  if (!strcmp(s,"tiedStatePos")){
    		  useTiedStatePos=TRUE;
    	  }
    	  else if (!strcmp(s, "three"))
            THREEACCS=TRUE;
         else if (!strcmp(s, "twodatafiles"))
            twoDataFiles = TRUE;
         else HError(1, "Unknown option %s",s);
         break;
      case 'u':
         if(!strcmp(s, "u"))
            SetuFlags(&uFlags);
         else if(!strcmp(s, "umle")){
            SetuFlags(&uFlagsMLE);
            if((uFlagsMLE&UPMEANS) && !(uFlagsMLE&UPVARS)){
               HError(-1, "Updating means with MLE but note that the formula used to get the value of the smoothing.. ");
               HError(-1, "constant D still uses discriminative mean estimate.  This is a minor issue in practice since most");
               HError(-1, "of the constants D are set to the denominator occupancy times E.  But be aware..." );
            }
         }
         else HError(1, "Unrecognised -u option %s", s);
         break;
      case 'v':
         minVar = GetChkedFlt(0.0,10.0,s); break;
      case 'w':
         if(!strcmp(s, "w"))
            mixWeightFloor = MINMIX * GetChkedFlt(0.0,10000.0,s); 
         else HError(1, "Unrecognised -w option %s", s);
         break;
      case 'x':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: HMM file extension expected");
         hmmExt = GetStrArg(); break;
      case 'y':  /*This is never used-- back-compat. */
         if (NextArg() != STRINGARG)
            HError(1, "HMMIRest: -y expects filename.");
         dictFN = GetStrArg(); 
         break;
      case 'z':
         if (NextArg() != STRINGARG)
            HError(2319,"HMMIRest: output TMF file expected");
         xfInfo.xformTMF = GetStrArg(); 
         break;
      case 'B':
         saveBinary=TRUE;
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(2319,"HMMIRest: Data File format expected");
         if((dff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2389,"HMMIRest: Warning ALIEN Data file format set");
         break;
      case 'H':
         if (NextArg() != STRINGARG)
            HError(2319,"HMMIRest: HMM macro file name expected");
         if(!strcmp(s,"H")){
           char *x = GetStrArg();
            AddMMF(&hset,x);
         } else if(!strcmp(s,"Hprior")){
           if(!hset_prior_initialised){
             hset_prior_initialised=TRUE;
             CreateHMMSet(&hset_prior,&hmmStack,TRUE);
           }
           AddMMF(&hset_prior, GetStrArg());
         } else if(!strcmp(s,"Hfullc")){
        	 CreateHMMSet(&hset_fullc,&hmmStack,TRUE);
        	 AddMMF(&hset_fullc, GetStrArg());
         }
         else HError(1, "Unknown option %s", s);
         break;            
      case 'I':
         if (NextArg() != STRINGARG)
            HError(2319,"HMMIRest: MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'M':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: Output macro file directory expected");
         newDir = GetStrArg();
         break;

      case 'Y':
    	  if (NextArg()!=STRINGARG)
    		  HError(2319,"HERest: -Y: Posterior directory expected");
    	  posDir[++numposdirs] = GetStrArg();
    	  if (NextArg()==STRINGARG)
    		  posExt = GetStrArg();
    	 break;

      case 'Q':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: -Q: Lattice extension expected");
         latExt = GetStrArg();
         break;
      case 'T':
         trace = GetChkedInt(0,0100000,s);
         break;
      /* additional options for transform support */
      case 'a':
        xfInfo.useInXForm = TRUE; break;
      case 'h':
        if (NextArg()!=STRINGARG)
          HError(1,"Speaker name pattern expected");
        xfInfo.outSpkrPat = GetStrArg();
        if (NextArg() != SWITCHARG)
          HError(2319,"HERest: cannot have -h as the last option");
        break;
      case 'E':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: parent transform directory expected");
         xfInfo.usePaXForm = TRUE;
         xfInfo.paXFormDir = GetStrArg();
         if (NextArg()==STRINGARG)
           xfInfo.paXFormExt = GetStrArg();
         if (NextArg() != SWITCHARG)
           HError(2319,"HMMIRest: cannot have -E as the last option");
         break;
      case 'J':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: input transform directory expected");
         AddInXFormDir(&hset,GetStrArg());
         if (NextArg()==STRINGARG)
           xfInfo.inXFormExt = GetStrArg();
         if (NextArg() != SWITCHARG)
           HError(2319,"HMMIRest: cannot have -J as the last option");
         break;
      case 'K':
         if (NextArg()!=STRINGARG)
            HError(2319,"HMMIRest: output transform directory expected");
         xfInfo.outXFormDir = GetStrArg();
         if (NextArg()==STRINGARG)
           xfInfo.outXFormExt = GetStrArg();
         if (NextArg() != SWITCHARG)
           HError(2319,"HMMIRest: cannot have -K as the last option");
         break;
      case 'L':
    	  if (NextArg()!=STRINGARG)
    		  HError(2319,"HMMIRest: MPE&MLE statistics script expected");
    	  statfn=GetStrArg();
    	  break;
      case 'W':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HMMIRest: -W: Posterior synthesiser model file expected");
    	 AddMMF(&synrhset,GetStrArg());
    	 useGmmSynr=TRUE;
    	 numposdirs=1;
    	 break;
      case 'e':
    	 if(NextArg()!=STRINGARG)
    		 HError(2319,"HMMIRest: -e: Posterior synthesiser list file expected");
    	 strcpy(synrListFn,GetStrArg());
    	 break;
      break;
      default:
         HError(2319,"HMMIRest: Unknown switch %s",s);break;
      }
   } /*matches while(NextArg() == SWITCHARG)*/


   InitVocab(&vocab); /* The actual dict is not needed, only the structure; this relates to HNet and reading lattices. */

   if (NextArg() != STRINGARG)
      HError(2319,"HMMIRest: file name of hmm list expected");
   if(useGmmSynr){
	   if(useTiedStatePos)
		   tspSynr=BuildClusteredGMMTree(&synrhset, synrListFn);
	   else
		   phnSynr=LoadPosSynth(&synrhset,synrListFn);
   }
   if(parMode>0)
	   assert(numposdirs==poswidth[0]);
   Initialise(GetStrArg()); /*GetStrArg() will return the hmmList.*/
   if(hset.fMPE&&(uFlags&UPFPROJ)&&parMode!=0){
	   if(statfn==NULL)
		   HError(2319,"HMMIRest: MPE&MLE statistics for fMPE regression parameters update expected");
	   printf("fMPE: load MPE&MLE statistics for fMPE regression parameters update\n");
	   LoadMPEMLEStat(statfn);
   }

   do { 
      char *accfn;
      if(parMode == 0){  /*The default is -1.  0 means gather together the parallel files.*/

         /*check no of files:*/
         if(THREEACCS && MMIPrior){   
            if((NumArgs() % 4) != 0)
               HError(1, "HMMIRest: nAccs should divide into 4 for MPE with MMI Prior");
         } else if(THREEACCS){
            if((NumArgs() % 3) != 0)
               HError(1, "HMMIRest: nAccs should divide into 3 for MEE&&MPEStoreML");
         } else if(uFlags&UPFPROJ){
        	 if((NumArgs() % 1) != 0)
        	   HError(1, "HMMIRest: nAccs should divide into 1 for fMPE");
         } else if(uFlags&UPSEMIT){
        	 if((NumArgs() % 1) != 0)
        	   HError(1, "HMMIRest: nAccs should divide into 1 for SEMIT");
         }
         else if (!ML_MODE && (NumArgs() % 2) != 0)/*need an even number because .1 and .2 files.*/
            HError(2319,"HMMIRest: Odd num of training files for single pass training");
       
         accfn = GetStrArg(); 

         /* Load first accumulator file (.1) */

         /* Warn if it is not a .acc.1 file */
         if(uFlags&UPFPROJ){
        	 if(!strstr(accfn,".acc.4")) HError(-1, "Expecting a *.acc.4 file, got %s", accfn);
        	 src=LoadAccsParallel(&hset, accfn,  uFlagsAccs, 0);   /*src is the still-open file.*/
         }
         else{
             if(!strstr(accfn,".acc.1")) HError(-1, "Expecting a *.acc.1 file, got %s", accfn);
             src=LoadAccsParallel(&hset, accfn,  uFlagsAccs, 0);   /*src is the still-open file.*/
         }

         ReadFloat(&src,&x,1,TRUE);    /*x is the average log prob (MPE->avg correctness) */
         ReadInt(&src,&i,1,TRUE);      /* i the number of timeperiods (MPE->num correct words).*/


         if(!MPE){ totalPr1 += x*i;  totalT += i; }
         else {  TotalCorr += x*i; TotalNWords += i; }

         CloseSource( &src );
         /* Load second accumulator file (.2) */

         if(!ML_MODE&&!(uFlags&UPFPROJ)&&!(uFlags&UPSEMIT)){                   /* ..Then load MMI acc. */
            accfn = GetStrArg();
            if(!strstr(accfn,".acc.2")) HError(-1, "Expecting a *.acc.2 file, got %s", accfn);
            src=LoadAccsParallel(&hset, accfn,  uFlagsAccs, 1);

            ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
            ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
            CloseSource( &src );

            /*pr2 contains the MMI stats.*/
            totalPr2 += x*i;  /*totalT += i;*/
         }

         /* Load third accumulator file (.3) */

         if(THREEACCS){             /* (MPE case). */
            accfn = GetStrArg();
            if(!strstr(accfn, ".3")) HError(1, "Error, expecting a HDR?.acc.3 file, got %s", accfn);
            src=LoadAccsParallel(&hset, accfn, uFlagsAccs, 2/*third position(1st==0) of hset*/); /*src is the still-open file.*/

            ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
            ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
            CloseSource( &src );

            totalPr1 += x*i;  totalT += i;
         }

         /* Load fourth accumulator file (.4) */

         if(MMIPrior){             /* (MPE case). */
            accfn = GetStrArg();
            if(!strstr(accfn, ".4")) HError(1, "Error, expecting a HDR?.acc.4 file, got %s", accfn);
            src=LoadAccsParallel(&hset, accfn, uFlagsAccs, 3/*fourth position(1st==0) of hset*/); /*src is the still-open file.*/

            ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
            ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
            CloseSource( &src );

            /*pr3 contains the MMI stats.*/
            totalPr3 += x*i;
         }
      } else {
         /*parMode not zero -> load data files & align..*/
         Boolean isPipe;
       
         if(NextArg() != STRINGARG)
            HError(2319,"HERest: data file name expected");
       
         if ( maxSnt != 0  && nSnt>maxSnt ) GetStrArg(); /*Pass over file. */
         else {   /* apply F-B.  */

            if (twoDataFiles){
               if ((NumArgs() % 2) != 0)
                  HError(2319,"HERest: Odd num of training files for single pass training");
               strcpy(datafn1,GetStrArg());
               datafn = datafn1;   datafn2 = GetStrArg();
            } else {
	       datafn = GetStrArg();
	       datafn2 = NULL;
	    }
	 
            if (UpdateSpkrStats(&hset,&xfInfo, datafn)) nSnt=0 ;
            fbInfo.inXForm = xfInfo.inXForm;
            fbInfo.paXForm = xfInfo.paXForm;

            /* derive lattice base file name from segment name using LATFILEMASK 
               this can be used to discard extra info (various cluster IDs, etc) */
            if (latFileMask) {
               if (!MaskMatch (latFileMask, datafn_lat, datafn))
                  HError(2319,"HERest: LATFILEMASK %s has no match with segemnt %s", latFileMask, datafn);
            }
            else
               strcpy (datafn_lat, datafn);

            if(nDenLats > 0){ /* Load denominator (recognition) lattices. */
               char buf1[1024],buf2[1024],buf3[1024];
               for(latn = 0; latn<nDenLats;latn++){
                  if ( denLatSubDirPat[0] ){
                     if ( !MaskMatch( denLatSubDirPat , buf1 , datafn_lat ) )
                        HError(2319,"HERest: mask %s has no match with segemnt %s" , denLatSubDirPat , datafn_lat );
                     MakeFN(buf1,denLatDir[latn],NULL,buf2);
                  }
                  else
                     strcpy(buf2,denLatDir[latn]);
                  if ( LatMask_Denominator != NULL ){
                     if ( !MaskMatch( LatMask_Denominator , buf1 , datafn_lat ) )
                        HError(2319,"HERest: mask %s has no match with segemnt %s" , LatMask_Denominator , datafn_lat );
                     MakeFN(buf1,buf2,NULL,buf3);
                     strcpy (buf2, buf3);
                  }
                  
                  if (useLLF) { 
                     denLats[latn] = GetLattice(datafn_lat,buf2, latExt,
                                                &latStack, &vocab, FALSE/*shortArc*/, TRUE/*add2Dict*/);
                  }
                  else {
                     MakeFN(datafn_lat,buf2,latExt,latfn);
                     f = FOpen(latfn, NetFilter, &isPipe);
                     if(!f) HError(1, "Couldn't open file %s\n", latfn);
                     printf("Reading lattice from file: %s\n", latfn); fflush(stdout);
                     denLats[latn] = ReadLattice(f, &latStack, &vocab, FALSE/*shortArc*/, TRUE/*add2Dict*/);
                     FClose(f, isPipe);
                  }
               }
            }

            if(nNumLats > 0){  /* Load numerator (correct transcription) lattices. */
               char buf1[1024],buf2[1024],buf3[1024];
               for(latn=0;latn<nNumLats;latn++){
                  if ( numLatSubDirPat[0] ){
                     if ( !MaskMatch( numLatSubDirPat , buf1 , datafn_lat ) )
                        HError(2319,"HERest: mask %s has no match with segemnt %s" , numLatSubDirPat , datafn_lat );
                     MakeFN(buf1,numLatDir[latn],NULL,buf2);
                  }
                  else
                     strcpy(buf2,numLatDir[latn]);
                  if ( LatMask_Numerator != NULL ){
                     if ( !MaskMatch( LatMask_Numerator , buf1 , datafn_lat ) )
                        HError(2319,"HERest: mask %s has no match with segemnt %s" , LatMask_Numerator , datafn_lat );
                     MakeFN(buf1,buf2,NULL,buf3);
                     strcpy (buf2, buf3);
                  }

                  if (useLLF) {
                     numLats[latn] = GetLattice(datafn_lat,buf2, latExt,
                                                &latStack, &vocab, FALSE/*shortArc*/, TRUE/*add2Dict*/);
                  }
                  else {
                     MakeFN(datafn_lat,buf2,latExt,latfn);
                     f = FOpen(latfn, NetFilter, &isPipe);
                     if(!f)  HError(1, "Couldn't open file %s\n", latfn);
                     printf("Reading lattice from file: %s\n", latfn); fflush(stdout);
                     numLats[latn] = ReadLattice(f, &latStack, &vocab, FALSE/*shortArc*/, TRUE/*add2Dict*/);
                     FClose(f, isPipe);
                  }
               }
            }
            if(hset.TVWR){
            	for(pos_stream=1;pos_stream<=numposdirs;pos_stream++){
            		if(posDir[pos_stream]!=NULL){
            			MakeFN(datafn,posDir[pos_stream],posExt,posfn);
            			printf("Reading posterior from file: %s\n", posfn); fflush(stdout);
            		}
            		if(phnSynr!=NULL){
            			printf("Synthesize posterior for file: %s\n", datafn); fflush(stdout);
            		}
            		LoadHDF4TVWR(&(fbInfo.hdfInfo), dff, posDir[pos_stream], posExt, datafn, &synrhset, phnSynr, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL, pos_stream);
            	}
            }
            if(hset.fMPE){
            	if(useTiedStatePos){
            		LoadHDF4fMPE2(&(fbInfo.hdfInfo), dff, datafn, tspSynr, hset.TVWR);
            	}
            	else{
            		LoadHDF4fMPE(&(fbInfo.hdfInfo), dff, posfn, datafn, &synrhset, phnSynr, hset.TVWR, xfInfo.inXForm?xfInfo.inXForm->xformSet->xforms[1]:NULL);
            	}
            }
            { /*apply F-B*/
               Boolean DoCorrectSentence,DoRecogLattice;
               int CorrIndex,RecogIndex1, RecogIndex2;
               DoCorrectSentence = !MPE || (MPE&&MPEStoreML);
               DoRecogLattice = !ML_MODE;

               CorrIndex = MPE&&!ML_MODE ? 2 : 0;   /* If MPE then the correct transcription ("mle" acc, for I smoothing) goes in position 2, if MMI then in 0. */
               RecogIndex1 = MPE ? 0 : 1;  /* If MPE then the first of the indices of the recognition lattice is the "num" acc (0).
                                              If MMI then it is the "den" acc (1). */
               RecogIndex2 = MPE ? 1 : 999;  /* If MPE then the second of the indices when aligning the recognition lat is the "den" acc,
                                                where arcs with negative differentials go.  If not MPE then it's a don't-care, and never read. */

               if(DoCorrectSentence && !nNumLats)  HError(-1, "No correct-transcription lattices specified so , use -q option.");
               if(DoRecogLattice && !nDenLats) HError(1, "No recognition lattices specified, use -r option.");
           

               if(DoCorrectSentence){
                  int i; 
                  for(i=0;i<nNumLats;i++)
                	  FBLatAddLattice(&fbInfo, numLats[i]);
                  FBLatFirstPass(&fbInfo, dff, datafn, datafn2, NULL/*MPE-related*/);
                  if(fbInfo.failedPass)
                	  printf("File %s is ignored\n",datafn);
                  FBLatSecondPass(&fbInfo, CorrIndex, 999/*dont-care*/);
                  totalT += fbInfo.T;
                  totalPr1 += fbInfo.pr;
               }

               if(DoRecogLattice){
                  int i,j; 
                  for(i=0;i<nDenLats;i++)
                	  FBLatAddLattice(&fbInfo, denLats[i]);
                  for(i=0;i<nNumLats;i++){
                     Boolean UseLat = TRUE; 
                     for(j=0;j<nDenLats;j++)
                    	 if (LatInLat(numLats[i],denLats[j]))
                    		 UseLat=FALSE; /*  Don't add redundant num lattices. */
                     if(UseLat){
                    	 if(trace&T_TOP) printf("[+num]");
                    	 FBLatAddLattice(&fbInfo, numLats[i]);
                     }
                  }
                  if(MMIPrior) SetDoingFourthAcc(TRUE,3);
                  FBLatFirstPass(&fbInfo, dff, datafn, datafn2, MPE ? numLats[0] : NULL);
                  if(fbInfo.failedPass)
                	  printf("File %s is ignored\n",datafn);
                  /* MPE only uses one of the num lats, if there are multiple ones (unlikely anyway) */
                  FBLatSecondPass(&fbInfo, RecogIndex1, RecogIndex2);
                  if(MMIPrior){   
                     SetDoingFourthAcc(FALSE,999);
                     totalPr3 += fbInfo.pr;
                  }
             
                  if(!DoCorrectSentence) totalT += fbInfo.T;
                  totalPr2 += fbInfo.pr;
                  if(MPE){  TotalNWords += fbInfo.MPEFileLength; TotalCorr += fbInfo.AvgCorr; }
               }

               nSnt++;
               if(hset.TVWR){
            	   ClearRwgtCache(&hset);
               }
               ResetHeap(&transStack);
               ResetHeap(&latStack);
            }
         }
      } /*[parMode]*/
   } while (NumArgs()>0);
   
   
   if (parMode>0 || (parMode==0 && (updateMode&UPMODE_DUMP))){
	   if(hset.fMPE&&(uFlags&UPFPROJ)){
		   MakeFN("HDR$.acc.4",newDir,NULL,newFn);
		   f=DumpAccsParallel(&hset,newFn,parMode, uFlagsAccs, 0);
		   x = TotalCorr/TotalNWords;
		   WriteFloat(f, &x, 1, TRUE);
		   WriteInt(f, &TotalNWords, 1, TRUE);
		   fclose(f);
	   }
	   else{
		  MakeFN("HDR$.acc.1",newDir,NULL,newFn);
		  f=DumpAccsParallel(&hset,newFn,parMode, uFlagsAccs, 0);

		  if(MPE){ /* .1 acc contains the MPE crit at the end. */
			 x = TotalCorr/TotalNWords;
			 WriteFloat(f, &x, 1, TRUE);
			 WriteInt(f, &TotalNWords, 1, TRUE);
		  }  else {
			 x = totalPr1/totalT;
			 WriteFloat(f,&x,1,TRUE);
			 WriteInt(f,&totalT,1,TRUE);
		  }
		  fclose( f );
		  if(!ML_MODE){ /* the dump .2 accs. */
			 MakeFN("HDR$.acc.2",newDir,NULL,newFn);
			 f=DumpAccsParallel(&hset,newFn,parMode, uFlagsAccs, 1);
			 x = totalPr2/totalT; /*MMI den prob in either MMI or MPE case.*/
			 WriteFloat(f,&x,1,TRUE);
			 WriteInt(f,&totalT,1,TRUE);
			 fclose( f );
		  }

		  if(THREEACCS){
			 MakeFN("HDR$.acc.3",newDir,NULL,newFn);
			 f=DumpAccsParallel(&hset,newFn,parMode, uFlagsAccs,2/*3rd position (numbered 2) on hset1(during alignment);*/);

			 x = totalPr1/totalT; /* This is where the MLE prob is stored. */
			 WriteFloat(f, &x, 1, TRUE);
			 WriteInt(f, &totalT, 1, TRUE);
			 fclose( f );
		  }

		  if(MMIPrior){   /* 4th acc for MMI den in MPE with MMI prior */
			 MakeFN("HDR$.acc.4",newDir,NULL,newFn);
			 f=DumpAccsParallel(&hset,newFn,parMode, uFlagsAccs, 3/*4th position (numbered 3) on hset1(during alignment);*/);

			 x = totalPr3/totalT; /* This is where the MMI den prob is stored. */
			 WriteFloat(f, &x, 1, TRUE);
			 WriteInt(f, &totalT, 1, TRUE);
			 fclose( f );
		  }
	   }
      if(trace&T_TOP) PrintCriteria();

   } 
   if (parMode <= 0){      /*parMode <= 0, so do re-estimation.*/
 	  if (hset.TVWR){
 		  printf("TVWR tied parameter?=%s\n",hset.tiedTVWR?"yes":"no");
 	  }
      if (stats) StatReport();
      if (uFlags&UPXFORM) {/* ensure final speaker correctly handled */
         UpdateSpkrStats(&hset,&xfInfo, NULL);
         if (trace&T_TOP) {
          printf("Reestimation complete - average log prob per frame = %e (%d frames)\n",
                 totalPr1/totalT, totalT);
        }
      } else
      if (updateMode&UPMODE_UPDATE)
        UpdateModels(); 
   }
   Exit(EXITSTATUS);
   return (0);          /* keep compiler happy */
}

/* -------------------------- Initialisation ----------------------- */

void Initialise(char *hmmListFn)
{  
   char buf[256];
   int s,regwsize=0;
   CreateHeap(&transStack,   "transStore",    MSTAK, 1, 0.5, 1000,  10000);
   CreateHeap(&accStack,   "accStore",    MSTAK, 1, 1.0, 50000,  500000);
   CreateHeap(&latStack,"latStore", MSTAK, 1, 1.0, 50000, 500000);


   /* Load HMMs and init HMMSet related global variables */
   MakeHMMSet( &hset, hmmListFn,dnn_label_list );/*read hmmlist and build logical list*/
   regwsize=posize*hset.temp_ctx_size;
   if(regwsize){
	   hset.TVWR=TRUE;
	   hset.regwSize=regwsize;
	   hset.pwidth[0]=poswidth[0];
	   for(s=1;s<=poswidth[0];s++){
		   hset.pwidth[s]=poswidth[s];
	   }
   }
   LoadHMMSet( &hset,hmmDir,hmmExt);
   SetParmHMMSet(&hset);
   if(hset_prior_initialised){
     MakeHMMSet( &hset_prior, hmmListFn,NULL );
     LoadHMMSet( &hset_prior, hmmDir,hmmExt);
   }
   if(MMIPrior) NumAccs=4;
   else if(ML_MODE) NumAccs=1;
   else if(THREEACCS/*MPE||MPEStoreML*/) NumAccs=3;
   else NumAccs=2;

   if(fMPE){
	   hset.fMPE=fMPE;
	   if(uFlags&UPFPROJ){
		   if(parMode==0){
			   THREEACCS=FALSE;
			   NumAccs=1;
		   }
		   else{
			   THREEACCS=TRUE;/*the third one only has the valid trans&weight but invalid mean&variance acc*/
			   NumAccs=3;
		   }
	   }
	   if(hset.fproj==NULL){
		   InitFprojection(&hset, numRegions, targetDim, inputDimO, inputDimP, octxwin, pctxwin, "FPROJ", hldaname);
		   CreateFprojMacro(&hset,hset.fproj, "FPROJ");
	   }
   }
   if(uFlags&UPSEMIT){
	   THREEACCS=FALSE;
	   uFlags = uFlags|UPMEANS|UPVARS;
	   NumAccs=1;
   }
   {
      uFlagsAccs =  uFlags|(uFlags&UPMEANS||uFlags&UPVARS ? UPMEANS|UPVARS : 0);  
      /*That modification to uFlags means: if either mean or var is updated, accumulate both.*/
      AttachAccsParallel(&hset, &accStack, uFlagsAccs, NumAccs);
      ZeroAccsParallel(&hset, uFlagsAccs, NumAccs); 
   }
   

   P = hset.numPhyHMM;
   L = hset.numLogHMM;
   vSize = hset.vecSize;
   S = hset.swidth[0];
   maxM = MaxMixInSet(&hset);
   hsKind = hset.hsKind;

   if(S > 1)
      HError(-1, "HMMIRest: Code is intended to support multiple streams but code has not been debugged.  Be warned!");
   /*Check that hset is the right kind.*/
   if(!(hsKind == PLAINHS || hsKind == SHAREDHS || hsKind == TIEDHS))
      HError(1, "HMMIRest: hset kind not PLAIN or SHARED or TIED.");
   if((!(hset.ckind == NULLC)) && (!(hset.ckind == DIAGC)))
      HError(1, "HMMIRest: cov kind not DIAGC.");
   
   /* Additional code for the adaptation updates */
   /*!Deleted*/

   /*Initialise those modules.*/
   InitialiseFBInfo(&fbInfo, &hset,   uFlags|(uFlags&UPMEANS||uFlags&UPVARS ? UPMEANS|UPVARS : 0), twoDataFiles);
   /*That modification to uFlags means: if either mean or var is updated, accumulate both.*/
 

   /* Set the variance floor */
   SetVFloor( &hset, vFloor, minVar); /*This sets the array minVar up... but it only works if
                                        minimum variances have been set within the hmmset (e.g,
                                        by HCompV.*/

   if (trace&T_TOP) {
      printf("HMMIRest  Updating: ");
      if (uFlags&UPTRANS) printf("Transitions "); 
      if (uFlags&UPMEANS) printf("Means "); 
      if (uFlags&UPVARS)  printf("Variances "); 
      if (uFlags&UPMIXES && maxM>1)  printf("MixWeights ");
      if (uFlags&UPRWGTS && maxM>1)  printf("RegWeights ");
      if (uFlags&UPXFORM)  {
        if (xfInfo.inSpkrPat == NULL) xfInfo.inSpkrPat = xfInfo.outSpkrPat;
        if (xfInfo.paSpkrPat == NULL) xfInfo.paSpkrPat = xfInfo.outSpkrPat;
        if (uFlags != UPXFORM)
          HError(999,"Can only update linear transforms OR model parameters!");
        xfInfo.useOutXForm = TRUE;
        /* This initialises things - temporary hack - THINK!! */
        CreateAdaptXForm(&hset, "tmp"); 
      }
      if (uFlags&UPFPROJ)  printf("Regressions parameters ");
      printf("\n ");
      if (parMode>=0) printf("Parallel-Mode[%d] ",parMode);
      printf(" - covkind is %s; ",CovKind2Str(hset.ckind,buf));
      printf(" system is ");
      switch (hsKind){
      case PLAINHS:  printf("PLAIN\n");  break;
      case SHAREDHS: printf("SHARED\n"); break;
      case TIEDHS:   printf("TIED\n");   break;
      case DISCRETEHS: printf("DISCRETE\n"); break;
      }

      printf("%d Logical/%d Physical Models Loaded, VecSize=%d\n",L,P,vSize);
      if (hset.numFiles>0)
         printf("%d MMF input files\n",hset.numFiles);
      if (mmfFn != NULL)
         printf("Output to MMF file:  %s\n",mmfFn); 

      if (uFlags&UPFPROJ){
    	  printf("fMPE: regression parameters update iteration: Steepest Gradient ascent method\n");
      }
      fflush(stdout);
   }

   if (parMode != 0){
      ConvDiagC(&hset,TRUE); /*!Converts DIAGC to INVDIAGC: changing over the CovKind in each HMM and inverting all the elements.*/
      ConvLogWt(&hset);
   }   
}
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
void LoadMPEMLEStat(char* fn){
	FILE* file;
	Source src;
	int numacc=0;
	char accfn[256];   /* buffer for current script arg */
	float x;
	int i;
	UPDSet uFlagsAccs = UPMEANS|UPVARS|UPTRANS|UPMIXES;   /* used in storing accs. */
	HMMScanState hss;
	int k;
	/* 3rd set of occs only required in MPE case where ML update has been specified for certain parms.. */
	float occ1,occ2,occ3,sqAcc1,sqAcc2,Acc1,Acc2, ivar, u,Snum,Sden;
	int vSize = hset.swidth[1];
	Vector mean;
	Covariance cov;
	VaAcc *va1,*va2,*va3;
	MuAcc *ma1,*ma2,*ma3;
	MixPDF *mp;
	Vector initchk=CreateVector(&gstack, vSize);
	if ((file = fopen(fn,"r")) == NULL){  /* Don't care if text/binary */
		HRError(5010,"LoadStatfMPE: Cannot open script file %s",fn);
	}
	while (ScriptWord(file,accfn) != NULL){
		++numacc;
	}
	rewind(file);
    /*check no of files:*/
    if((numacc % 3) != 0)
    	HError(1, "HMMIRest: nAccs should divide into 3 for MEE&&MPEStoreML");
    totalPr1=0;
    totalT=0;
    totalPr2=0;
    TotalCorr=0;
    TotalNWords=0;

	do{
        /* Load numerator accumulator file (.1) */
		ScriptWord(file,accfn);numacc--;
        if(!strstr(accfn,".acc.1")) HError(1, "Expecting a *.acc.1 file, got %s", accfn);
        src=LoadAccsParallel(&hset, accfn,  uFlagsAccs, 0);   /*src is the still-open file.*/
        ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
        ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
        TotalCorr += x*i; TotalNWords += i;
        CloseSource( &src );

        /* Load denominator accumulator file (.2) */
        ScriptWord(file,accfn);numacc--;
        if(!strstr(accfn,".acc.2")) HError(1, "Expecting a *.acc.2 file, got %s", accfn);
        src=LoadAccsParallel(&hset, accfn,  uFlagsAccs, 1);
        ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
        ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
        totalPr2 += x*i;  /*totalT += i;*/
        CloseSource( &src );

        /* Load third accumulator file (.3) */
        ScriptWord(file,accfn);numacc--;
        if(!strstr(accfn, ".acc.3")) HError(1, "expecting a HDR?.acc.3 file, got %s", accfn);
        src=LoadAccsParallel(&hset, accfn, uFlagsAccs, 2/*third position(1st==0) of hset*/); /*src is the still-open file.*/
        ReadFloat(&src,&x,1,TRUE); /*!x must be the average log prob*/
        ReadInt(&src,&i,1,TRUE); /*!and i the number of timeperiods.*/
        totalPr1 += x*i;  totalT += i;
        CloseSource( &src );
	}while(numacc>0);
	fclose(file);
	if(trace&T_TOP) PrintCriteria();
    totalPr1=0;
    totalT=0;
    totalPr2=0;
    TotalCorr=0;
    TotalNWords=0;
	/*RestoreAccsParallel(&hset, 0);  /*Removes mu offsets, 0 is index of numerator */
	/*RestoreAccsParallel(&hset, 1);/*Removes va offsets, 1 is index of denominator */
	/*RestoreAccsParallel(&hset, 2);*/

	NewHMMScan(&hset, &hss);
	while (GoNextMix(&hss,FALSE)) {
		mp=hss.mp;
		cov = mp->cov; /*type Covariance. */
		va1 = GetHook(cov.var);
		va2 = va1+1;
		va3 = va1+2;
		mean = mp->mean;
		ma1 = GetHook(mean);/*numerator acc*/
		ma2 = ma1+1;/*denominator acc*/
		ma3 = ma1+2;

		occ1 = va1->occ;
		occ2 = va2->occ;
		occ3 = va3->occ;
		for (k=1; k<=vSize; k++){ /*For each vector component.*/
			sqAcc1 = va1->cov.var[k];
			sqAcc2 = va2->cov.var[k];
			Acc1   = ma1->mu[k];
			Acc2   = ma2->mu[k];
			u	   = mean[k];
			ivar = mp->ckind == INVDIAGC?cov.var[k]:1/cov.var[k];
			/*ma3->udif[k]=ivar*(Acc1-Acc2-u*(occ1-occ2));
			Snum=(sqAcc1-2*Acc1*u+occ1*u*u);
			Sden=(sqAcc2-2*Acc2*u+occ2*u*u);
			va3->vdif[k]=0.5*(Snum*ivar*ivar-occ1*ivar)-0.5*(Sden*ivar*ivar-occ2*ivar);*/
			ma3->udif[k]=ivar*(Acc1-Acc2);
			va3->vdif[k]=0.5*(ivar*ivar*(sqAcc1-sqAcc2)-ivar*(occ1-occ2));
			/*ma3->udif[k]*=(1/occ3);
			va3->vdif[k]*=(2/occ3);*/
		}
	}
	EndHMMScan(&hss);

	ZeroVector(initchk);
	NewHMMScan(&hset, &hss);
	while (GoNextMix(&hss,FALSE)) {
		mp=hss.mp;
		cov = mp->cov; /*type Covariance. */
		va1 = GetHook(cov.var);
		va2 = va1+1;
		va3 = va1+2;
		mean = mp->mean;
		ma1 = GetHook(mean);/*numerator acc*/
		ma2 = ma1+1;/*denominator acc*/
		ma3 = ma1+2;

		occ3 = va3->occ;
		if(occ3==0)continue;
		for (k=1; k<=vSize; k++){ /*For each vector component.*/
			initchk[k]+=(ma3->udif[k]*occ3+2*va3->vdif[k]*ma3->mu[k])/occ3;
		}
	}
	EndHMMScan(&hss);
	printf("init check total indirect differential:\n");
	for(k=1;k<=vSize;k++){
		printf("ind[%d]=%e\n",k,initchk[k]);
	}
	FreeVector(&gstack,initchk);

}

/* ------------------- Statistics Reporting  -------------------- */

/* PrintStats: for given hmm */
void PrintStats(FILE *f, int n, HLink hmm, int numEgs)
{
   WtAcc *wa;
   char buf[MAXSTRLEN];
   StateInfo *si;
   int i,N;
    
   N = hmm->numStates;
   ReWriteString(HMMPhysName(&hset,hmm),buf,DBL_QUOTE);
   fprintf(f,"%4d %14s %4d ",n,buf,numEgs);
   for (i=2;i<N;i++) {
      si = hmm->svec[i].info;
      wa = (WtAcc *)((si->pdf+1)->hook);
      fprintf(f," %10f",wa->occ);
   }
   fprintf(f,"\n");
}

/* StatReport: print statistics report */
void StatReport(void) /*This is used by other programs so I have had to change it back to how it was.*/
{
   HMMScanState hss;
   HLink hmm;
   FILE *f;
   int px;

   if ((f = fopen(statFN,"w")) == NULL){
      HError(2311,"StatReport: Unable to open stats file %s",statFN);
      return;
   }
   NewHMMScan(&hset,&hss);
   px=1;
   do {
      hmm = hss.hmm;
      PrintStats(f,px,hmm,(int)hmm->hook);
      px++;
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);
   fclose(f);
}

/* --------------------------- Model Update --------------------- */

static int nFloorVar = 0;     /* # of floored variance comps */
static int nFloorVarMix = 0;  /* # of mix comps with floored vars */
static long int nVar = 0;     /* # total of vars */
static int nMix = 0;  /*  total # of gaussians. */

static int nWeight = 0;  /*  total # of weights. */
static int nFloorWeight = 0;  /* # of floored weights. */


/* FloorMixes: apply floor to given mix set */
void FloorMixes(MixtureElem *mixes, int M, float floor)
{
   float sum,fsum,scale;
   MixtureElem *me;
   int m;
   
   sum = fsum = 0.0;
   for (m=1,me=mixes; m<=M; m++,me++) {
      if (me->weight>floor)
         sum += me->weight;
      else {
         fsum += floor; me->weight = floor; nFloorWeight++;
      }
      nWeight++;
   }
   if (fsum>1.0) HError(2327,"FloorMixes: Floor sum too large");
   if (fsum == 0.0) return; /*all >= floor*/
   if (sum == 0.0) HError(2328,"FloorMixes: No mixture weights above floor");
   /*all <= floor*/

   scale = (1.0-fsum)/sum; /*fsum is the sum of floor * (no of vals <= floor), sum is sum of >= floor.*/

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
         nFloorWeight++;
      }
      nWeight++;
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



Boolean SolveQuadratic(double a, double b, double c, double *ans1, double *ans2){
   double temp;

   if(a==0){ /* bx + c = 0 --> x = -c/b */
      (*ans1)=(*ans2) = -c/b;
      return TRUE;
   }

   temp  = b*b - (4*a*c);
   if(temp<0){
      if(b != 0 && fabs(temp/(b*b)) < 0.0001){
         *ans1=*ans2=(-b/(2*a)); /*Probably due to arithmetic error*/
         return TRUE;
      }
      else
         return FALSE;
   }
   temp = sqrt(temp);

   *ans1 = (temp-b)/(2*a);
   *ans2 = (-temp-b)/(2*a);
   return TRUE;
}


void UpdateWeightsOrTrans(int M, float *acc1, float *acc2, float *mixes, float *oldMixes, float C);



/* UpdateTrans: use acc values to calc new estimate for transP */
void UpdateTran(int px, HLink hmm)
{
   int i,j,N=hmm->numStates;
   TrAcc *ta1, *ta2, *ta3;
   Vector NewWghts = CreateVector(&gstack, N);
   Vector OldWghts = CreateVector(&gstack, N);

   
   ta1 = GetHook(hmm->transP);
   ta2 = (ML_MODE?NULL:ta1+1);
   ta3 = (THREEACCS?ta1+2:NULL); /*non-NULL in MPE case, where it is the ML accs. */

   if (ta1==NULL) return;   /* already done */
   N = hmm->numStates;
   for (i=1;i<N;i++) {
      if ((ta3 ? ta3->occ[i]:ta1->occ[i]) > MinOccTrans) {   /* more than MinOccTrans ML stats */
         for(j=1;j<=N;j++){
            NewWghts[j] = OldWghts[j] = (hmm->transP[i][j]>MINEARG?exp(hmm->transP[i][j]):0.0);
         }
         if(uFlagsMLE&UPTRANS && ta2) for(j=1;j<=N;j++) ta2->tran[i][j] = 0.0;
         UpdateWeightsOrTrans(N, ta1->tran[i], ta2?ta2->tran[i]:NULL, NewWghts, OldWghts, CTrans);
         for(j=1;j<=N;j++) 
           if(NewWghts[i] == 0 && OldWghts[i] != 0)
             HError(-1, "Transitions going to zero: advise setting e.g. ISMOOTHTAUT = 10 ");
         for(j=1;j<=N;j++)
            hmm->transP[i][j]=(NewWghts[j]>0.0?log(NewWghts[j]):LZERO);
      }
   }
   
   SetHook(hmm->transP,NULL);
   Dispose(&gstack, NewWghts);
}
  




float GiveDimMixD(MixPDF *mp, int k, int priortype){
   /*Gives the constant D, used in the EBW updates.*/

   float occ1, occ2, D=0;
   double temp1, temp2;
   VaAcc *va1=NULL, *va2=NULL;
   MuAcc *ma1=NULL, *ma2=NULL;
   float var_acc1, var_acc2, mu_acc1, mu_acc2, old_mu, old_sigmasq, a,b,c;

   if (priortype == 1){  /* MMI Prior */
      ma1=((MuAcc*)GetHook(mp->mean)) + 2;
      ma2=ma1+1;
      va1=((VaAcc *)GetHook(mp->cov.var)) + 2;
      va2=va1+1;
   } else if (priortype == 0){  /* standard MPE */
      ma1=GetHook(mp->mean);
      ma2=ma1+1;
      va1=GetHook(mp->cov.var);
      va2=va1+1;
   } else
      HError(999,"GiveMixD: Wrong prior type!");
                         
   occ1=ma1->occ; occ2=ma2->occ;
   if(occ2==0 || (uFlagsMLE&UPVARS)) return 0;

   var_acc1=va1->cov.var[k]; var_acc2=va2->cov.var[k];
   mu_acc1=ma1->mu[k]; mu_acc2=ma2->mu[k];
   old_mu=mp->mean[k]; old_sigmasq=mp->cov.var[k];
  
   /*Equation is:
     sigma_new = (var_acc1 - var_acc2 + D(old_sigmasq + old_mu^2)) / (occ1-occ2+D)    -    mu_new^2
     where mu_new = (mu_acc1 - mu_acc2 + D.old_mu) / (occ1-occ2+D)
    
     so, rearranging, we have, where sigma_new = 0,
     (var_acc1 - var_acc2 + D(old_sigmasq+old_mu^2)) * (occ1-occ2+D) = (mu_acc1-mu_acc2+D.old_mu)*(mu_acc1-mu_acc2+D.old_mu)
     =>   (var_acc1-var_acc2)(occ1-occ2) + (var_acc1-var_acc2)(D) + D(old_sigmasq+old_mu^2)(occ1-occ2) + D^2(old_sigmasq+old_mu^2)
     = (mu_acc1-mu_acc2)^2 + 2D(mu_acc1-mu_acc2)(old_mu) + D^2(old_mu^2)
    
     ==>
     D^2 (old_sigmasq+old_mu^2-old_mu^2) + 
     D [ ((old_sigmasq+old_mu^2)*(occ1-occ2))  + var_acc1 - var_acc2 - 2old_mu(mu_acc1-mu_acc2)]
     + (var_acc1-var_acc2)(occ1-occ2) - (mu_acc1 - mu_acc2)^2
     = 0
   */
  
   a = old_sigmasq;
   b = (old_sigmasq+old_mu*old_mu)*(occ1-occ2) + var_acc1 - var_acc2 - 2*old_mu*(mu_acc1-mu_acc2);
   c = (var_acc1-var_acc2)*(occ1-occ2) -  (mu_acc1 - mu_acc2)*(mu_acc1 - mu_acc2);
  
   if(SolveQuadratic(a,b,c,&temp1,&temp2)) /*If there's a solution...*/
      D = MAX(temp1,temp2) * DFACTOR;

  
   totalConst++;
   if(D > occ2*E) nonFlooredConst++;
  
   D=MAX(D,occ2*E);

   if(! FINITE(D))
      HError(1, "NaN in GiveMixD");
  
   return D;
}

float GiveMixD(MixPDF *mp, int stream, int priortype){
   /*Gives the constant D, used in the EBW updates.*/
   /*Returns FALSE if no use of variance/means.*/

   float  D=0;
   int k, vSize = hset.swidth[stream];

   if(ML_MODE) return 0;

   for(k=1;k<=vSize;k++){
      D = MAX(D, GiveDimMixD(mp,k,priortype));
   }
   return D;
}


Boolean UpdateGauss(int stream, MixPDF *mp){
   float D;
   int k;
   /* 3rd set of occs only required in MPE case where ML update has been specified for certain parms.. */
   float occ1,occ2,occ3,sqAcc1,sqAcc2,sqAcc3,Acc1,Acc2,Acc3, oldVar,newVar,newMean; 
   int vSize = hset.swidth[stream];
   Vector mean;
   Covariance cov;
   Vector minV = vFloor[stream]; /*!The vFloor for this stream [vFloor a global var.]*/
   VaAcc *va1,*va2,*va3;
   MuAcc *ma1,*ma2,*ma3;
   Boolean mixFloored=FALSE;

   cov = mp->cov; /*type Covariance. */
   va1 = GetHook(cov.var);
   va2 = (ML_MODE?NULL:va1+1);
   va3 = (THREEACCS?va1+2:NULL);
   mean = mp->mean;
   ma1 = GetHook(mean);/*numerator acc*/
   ma2 = (ML_MODE?NULL:ma1+1);/*denominator acc*/
   ma3 = (THREEACCS?ma1+2:NULL);

   occ1 = va1->occ;
   occ2 = (va2?va2->occ:0.0);
   occ3 = (va3?va3->occ:0.0);

   /* already checked that there is enough occupancy. */    


   if(fabs(va1->occ - ma1->occ) > 0.1) HError(1, "Mean and variance occupation counts differ.  Means and variances must not be separately tied.");


   /*if((va3 ? va3->occ : va1->occ) < MinOcc) return FALSE; /* if ML occ < MinOcc dont update. */
   if((va3&&ISmoothTau ? va3->occ : va1->occ) < MinOcc) return FALSE; /* if ML occ < MinOcc dont update. */

   if(mp->ckind != DIAGC)
      HError(999,"UpdateGauss: unknown ckind %d",mp->ckind);


   D= GiveMixD(mp,stream,0); /*Use constants on a mix level.*/


   for (k=1; k<=vSize-hset.projSize; k++){ /*For each vector component.*/
      sqAcc1 = va1->cov.var[k]; sqAcc2 = va2?va2->cov.var[k]:0.0; sqAcc3 = va3?va3->cov.var[k]:0.0; oldVar = cov.var[k];
      Acc1 = ma1->mu[k], Acc2 = ma2?ma2->mu[k]:0.0; Acc3 = ma3 ? ma3->mu[k]:0.0;
    

	  newMean =
		 uFlags&UPMEANS  ?
		 /*ML:*/        (uFlagsMLE&UPMEANS ? (ma3 ? Acc3/occ3 : Acc1/occ1) :
						 /*MMI/MPE:*/   (Acc1 - Acc2 + D*mean[k])/(occ1-occ2+D)) :
		 /*old: */   mean[k];
	  newVar =/*if MMI/MPE, newVar=(((sqAcc1-sqAcc2)+D*(cov.var[k]+mean[k]*mean[k]))/((occ1-occ2)+D))-newMean*newMean*/
		 uFlags&UPVARS ?
		 /*ML:*/            (uFlagsMLE&UPVARS ? (ma3 ? (sqAcc3 - 2*Acc3*newMean + occ3*newMean*newMean)/occ3 : (sqAcc1 - 2*Acc1*newMean + occ1*newMean*newMean)/occ1) :
							 /*MMI/MPE:*/       ((sqAcc1 - 2*Acc1*newMean + occ1*newMean*newMean) -  (sqAcc2 - 2*Acc2*newMean + occ2*newMean*newMean)
												 +  D*(cov.var[k] + (mean[k]-newMean)*(mean[k]-newMean))) / (occ1-occ2+D)) :
		 /*old:*/          cov.var[k];

      if (newVar<minV[k]) {
         newVar = minV[k]; nFloorVar++; mixFloored = TRUE;
      }
      nVar++;
      mean[k] = newMean;
      cov.var[k] = newVar;
   }
   if (mixFloored == TRUE) nFloorVarMix++;
   nMix ++;
   return TRUE;
}

float min(float a, float b){
	return a>b?b:a;
}
/*static Boolean isSil=FALSE;*/
/**
 * A solver for \sum_i^N frac{num_i}{den_i-lambda}=1
 */
Boolean NewtonMethod(DVector num, DVector den, double *xlambda){
	double lambda0,lambda,f_lambda,f_lambda_dif,temp,min_deni,min_dmni,lowerUpperBound;
	const double f_goal=0.0001;
	const double minor_shift=0.0001;
	int i;
	double num_occ;
	Boolean found;
	int iter;
	int N=DVectorSize(num);
	min_deni=1E100;
	num_occ=0;
	min_dmni=1E100;
	for(i=1;i<=N;i++){
		num_occ+=num[i];
		if((den[i]<min_deni)&&(num[i]>1E-100))
			min_deni=den[i];
		temp=den[i]-num[i];
		if((temp<min_dmni)&&(num[i]>1E-100))
			min_dmni=temp;
	}
	min_deni=min_deni-minor_shift;
	lowerUpperBound=min(min_dmni,min_deni);
	lambda0=-num_occ;
	lambda=min(lambda0,lowerUpperBound);
	f_lambda=-1;
	f_lambda_dif=0;
	for(i=1;i<=N;i++){
		temp=den[i]-lambda;
		f_lambda+=((num[i]<=1E-100)?0:num[i]/temp);
		f_lambda_dif+=((num[i]<=1E-100)?0:num[i]/(temp*temp));
	}
	iter=0;
	while(fabs(f_lambda)>f_goal&&++iter<=100){/*Search root of f(lambda)=\sum_i{num_i/(den_i*N-lambda)}-1*/
		if(f_lambda_dif==0){
			break;
		}
		lambda=lambda-f_lambda/f_lambda_dif;/*lambda=lambda-f(lambda)/f'(lambda)*/
		lambda=min(lambda,lowerUpperBound);
		f_lambda=-1;
		f_lambda_dif=0;
		for(i=1;i<=N;i++){
			temp=den[i]-lambda;
			f_lambda+=((num[i]<=1E-100)?0:num[i]/temp);/*f(lambda)=\sum_i{num_i/(den_i*N-lambda)}-1*/
			f_lambda_dif+=((num[i]<=1E-100)?0:num[i]/(temp*temp));/*f'(lambda)=\sum_i{num_i/(den_i*N-lambda)^2}*/
		}
		/*if(isSil){
			printf("f=%f, df=%f, lambda=%f, upper=%f\n",f_lambda,f_lambda_dif,lambda,min_deni);
		}*/
	}
	if(fabs(f_lambda)<=f_goal&&lambda<=lowerUpperBound){/*bug found from fabs(f_lambda<=f_goal)*/
		found=TRUE;
	}
	else{
		found=FALSE;
	}
	*xlambda=lambda;
	return found;
}
/**
 * Note that the numerator is I-Smoothed, for the denominator, if pre-cancel of static weight is done, let cur_swgt be that static weight.
 * If in tied system, no pre-cancelling is needed so that keep cur_swgt always ONE.
 * f(\lambda)=\sum_i {n_i/(d_i*N-\lambda)} -1
 */

void optimizerwjmi(HMMSet *hset, RWtAcc* num, RWtAcc* den, RWtAcc* ml, RegWgt* rwgt){
	double lambda;
	float wjmi;
	int i,s, _s, N,offset=0;
	Boolean found;
	Vector memMark=CreateVector(&gstack,1);
	for(s=1;s<=hset->pwidth[0];s++){
		for(_s=1;_s<=hset->temp_ctx_size;_s++){
			N=hset->rkind>=HIGH_SPA_REG?hset->confuTbl[rwgt->spIdx][_s][0]:hset->pwidth[s];
			DVector temp_num=CreateDVector(&gstack,N);
			DVector temp_den=CreateDVector(&gstack,N);
			for(i=1;i<=N;i++){
				/*if(hset.tiedTVWR)*/
				temp_num[i]=num->hjm[offset+i]+(ISmoothTauRWeightsSet&&ml->occ>0?(ISmoothTauRWeights*ml->hjm[offset+i]/ml->occ):0);
				/*else
					temp_num[i]=num->hjm[i]+(ISmoothTauWeightsSet&&mlstateocc>0?(ISmoothTauWeights*ml->hjm[i]/mlstateocc):0);*/
				wjmi=rwgt->regparm[offset+i];
				temp_den[i]=wjmi>MINMIX?den->hjm[offset+i]/wjmi:0;
			}
			found=NewtonMethod(temp_num, temp_den, &lambda);
			if(found){/*MPE updating*/
				for(i=1;i<=N;i++){
					rwgt->regparm[offset+i]=temp_num[i]/(temp_den[i]-lambda);
				}
			}
			offset+=N;
		}
	}
	FreeVector(&gstack,memMark);
}



/**
 * Remember to do pre-ISmoothing before application
 */
void optimizercjm(StreamElem *ste){
	double lambda;
	Boolean found;
	int m,M=ste->nMix;
	MixtureElem *me;
	WtAcc *wa1,*wa2,*wa3;
	wa1 = (WtAcc*)ste->hook;
	wa2 = wa1+1;
	wa3 = wa2+1;
	DVector temp_num=CreateDVector(&gstack,M);
	DVector temp_den=CreateDVector(&gstack,M);
	for(m=1;m<=M;m++){
		me=ste->spdf.cpdf+m;
		temp_num[m]=wa1->c[m]+((ISmoothTauWeightsSet&&wa3->occ)?ISmoothTauWeights*wa3->c[m]/wa3->occ:0);
		temp_den[m]=me->sweight>MINMIX?wa2->c[m]/me->sweight:0;
	}
	found=NewtonMethod(temp_num, temp_den, &lambda);
	if(found){
		for(m=1;m<=M;m++){
			me=ste->spdf.cpdf+m;
			me->weight=me->sweight=temp_num[m]/(temp_den[m]-lambda);
		}
	}
	FreeDVector(&gstack,temp_num);
}


void UpdateRWeight(HMMSet *hset, int px, HLink hmm)
{
   int s,m,M=0,N,S,j;
   WtAcc *wa1,*wa2,*wa3;
   RWtAcc *rwa1,*rwa2,*rwa3;
   StateElem *se;
   StreamElem *ste;
   MixtureElem *me;
   N = hmm->numStates;
   se = hmm->svec+2;
   S = hset->swidth[0];
   for (j=2; j<N; j++,se++){
      ste = se->info->pdf+1;
      for (s=1;s<=S; s++,ste++){
         wa1 = (WtAcc *)ste->hook;
         wa2 = wa1+1;
         wa3 = wa1+2;
         M=ste->nMix;
         if(wa1&&wa3->occ>MinOccWeights)/*update c_jm*/
        	 optimizercjm(ste);
         for(m=1;m<=M;m++){
        	 me=ste->spdf.cpdf+m;
        	 rwa1 = (RWtAcc*)me->rwgt->hook;
        	 rwa2 = rwa1+1;
        	 rwa3 = rwa1+2;
        	 if(rwa1&&rwa3->occ>MinOccWeights)/*update w_jmi*/
        		 optimizerwjmi(hset, rwa1, rwa2, rwa3, me->rwgt);
        	 me->rwgt->hook=NULL;
         }
         ste->hook=NULL;/*this is very important to make it right*/
      }
   }
}

void UpdateTiedTVWR(HMMSet *hset, int px, HLink hmm)
{
	int m,M=0,N,j,s;
	TVWRElem *ze;
	RWtAcc *rwa1,*rwa2,*rwa3;
	RegWgt **rwgts;
	N = hmm->numStates;
	ze = hmm->zvec+2;
	/*char *name=HMMPhysName(hset,hmm);*/
	for (j=2; j<N; j++,ze++){
		for(s=1;s<=hset->swidth[0];s++){
			rwgts=ze->info->tstr[s].rwgts;
			M=ze->info->tstr[s].nMix;
			for(m=1;m<=M;m++){
				rwa1=(RWtAcc *) (rwgts[m]->hook);
				rwa2 = rwa1+1;
				rwa3 = rwa1+2;
				if(rwa1&&rwa3->occ>MinOccWeights)/*update w_jmi*/
					optimizerwjmi(hset, rwa1, rwa2, rwa3, rwgts[m]);
				rwgts[m]->hook=NULL;
			}
		}
	}
}


void UpdateTrans(void){
   HMMScanState hss;
   HLink hmm;
   int px,n;
   void UpdateTran(int px, HLink hmm);

   NewHMMScan(&hset,&hss);
   px=1;
   do{
      hmm=hss.hmm;
      /*n = (int)hmm->hook; /*The number of training egs seen*/
      /* n is NO LONGER USED. */
      UpdateTran(px,hmm);
      px++;
   }while (GoNextHMM(&hss));
   EndHMMScan(&hss);
}

void UpdateWeights(void){
   HMMScanState hss;

   NewHMMScan(&hset,&hss);
   while(GoNextStream(&hss,FALSE)){
	   UpdateWeight(hss.s, hss.ste);
   }
   EndHMMScan(&hss);
}



void UpdateRWeights(void){
   HMMScanState hss;
   HLink hmm;
   int px,n;

   NewHMMScan(&hset,&hss);
   px=1;
   do{
      hmm=hss.hmm;
      n = (int)hmm->hook; /*The number of training egs seen*/
      /* n is NO LONGER USED. */
      if(hset.tiedTVWR){
    	  UpdateTiedTVWR(&hset, px, hmm);
      }
      else{
    	  UpdateRWeight(&hset, px, hmm);
      }
      px++;
   }while (GoNextHMM(&hss));
   EndHMMScan(&hss);
   /*printf("TVWR update stats:MPE update count=%d, ML update count=%d, skip update count=%d, cur average obj=%f, prv average obj=%f\n",mpe_up_count,ml_up_count,no_up_count,curAvgObj/mpe_up_count,prvAvgObj/mpe_up_count);*/
}

static void FixHMMForICrit();


static void FixWeightsForICrit(float Tau, Boolean THREEACCS){
   HMMScanState hss;
   NewHMMScan(&hset,&hss); 
   while(GoNextStream(&hss,FALSE)){
      WtAcc *wa_src, *wa_dst; int m,M;
      M = hss.M; wa_dst = (WtAcc*)hss.ste->hook; wa_src = (THREEACCS ? wa_dst+2 : wa_dst); /* THREEACCS should be true for the forseeable use of this. */

      for(m=1;m<=M;m++) wa_dst->c[m] += Tau * (wa_src->occ ? wa_src->c[m]/wa_src->occ : 1/M);
      if(!wa_src->occ) HError(-1, "wa_src->occ zero, in FixWeightsForICrit.");
      wa_dst->occ += Tau;
   }
   EndHMMScan(&hss); 
}

static void FixTransForICrit(float Tau, Boolean THREEACCS){
   HMMScanState hss;
   NewHMMScan(&hset,&hss); 
   do{
      TrAcc *ta_src, *ta_dst; int m,M;
      ta_dst = (TrAcc*) GetHook(hss.hmm->transP); ta_src = (THREEACCS ? ta_dst+2 : ta_dst); /* THREEACCS should be true for the forseeable use of this. */
      M=hss.hmm->numStates;
      for(m=1;m<M;m++){
         int n;
         if(ta_src->occ[m] != 0){
            for(n=1;n<=M;n++){
               ta_dst->tran[m][n] += Tau /*not * M!*/  *  ta_src->tran[m][n]/ta_src->occ[m];
            }
            ta_dst->occ[m] += Tau;
         }
      }
   } while(GoNextHMM(&hss));
   EndHMMScan(&hss); 
}

/* Calclulate MMI acc and save it in the orignial ML acc position, i.e., 3rd acc */
static void GetMMIAccMix(int stream, MixPDF *mp)
{
   int i,k,vSize;
   float occ1,occ2,D,s,mmimean,mmivar;
   Vector minV = vFloor[stream]; /*!The vFloor for this stream [vFloor a global var.]*/
   MuAcc *ma1, *ma2;
   VaAcc *va1, *va2;

   vSize = hset.swidth[stream];

   /* num for MMI is the ML acc */
   ma1 = ((MuAcc*)GetHook(mp->mean)) + 2; 
   va1 = ((VaAcc*)GetHook(mp->cov.var)) + 2;
   occ1 = ma1->occ;
   /* den for MMI is the 4th acc */
   ma2 = ((MuAcc*)GetHook(mp->mean)) + 3; 
   va2 = ((VaAcc*)GetHook(mp->cov.var)) + 3;
   occ2 = ma2->occ;

   if (occ1 >= MinOcc){   /* if ML occ< MinOcc dont do I-smoothing, keep ML in ma1 */

      /* I-smoothing for MMI acc: count smoothing, save in ML acc */
      if(occ1>0.001 && MMITauI > 0){
         s = MMITauI/occ1;
         for(i=1;i<=vSize;i++){
            ma1->mu[i] += s*ma1->mu[i];
            va1->cov.var[i] += s*va1->cov.var[i]; 
         }
         occ1 += MMITauI;
      }
      
      /* calculate D */
      D= GiveMixD(mp,stream,1); /*Use constants on a mix level. 1 means use 3 and 4 acc to calculate D*/
      
      /* calulate MMI statistics and save in 3rd acc */ 
      for(k=1; k<=vSize; k++){ /*For each vector component.*/
         /* MMI statistics */
         mmimean = (ma1->mu[k]-ma2->mu[k]+D*mp->mean[k])/(occ1-occ2+D);
         mmivar = (va1->cov.var[k]-va2->cov.var[k]+D*(mp->cov.var[k]+mp->mean[k]*mp->mean[k]))/(occ1-occ2+D); /* - mmimean*mmimean; omitted as in below it has to be added again */
         if (mmivar<minV[k] + mmimean*mmimean)  mmivar = minV[k] + mmimean*mmimean;      
         /* save to 3rd acc, the occupancy 100 does not matter as it will be divided in following I-smoothing procedure*/
         ma1->mu[k] = 100*mmimean;
         va1->cov.var[k] = 100*mmivar;
      }
      ma1->occ = 100;
      va1->occ = 100;
   } 
}


static void _FixHMMForICrit(float Tau, Boolean THREEACCS){
   /* Normally both Mu and Var occupancy will be identical. */
   HMMScanState hss;
   NewHMMScan(&hset,&hss); 

   while(GoNextMix(&hss,FALSE)){
      int i;
      MixPDF *mp;
      MuAcc *ma_src, *ma_dst; VaAcc *va_src, *va_dst;
      float s;
      mp = hss.mp;

      if(MMIPrior) GetMMIAccMix(hss.s, mp);

      /* _src refers to the place the MLE accs are stored [ hset,index 0 for MMI, hset,index 1 for MEE ] 
         _dst refers to the place we add the MLE accs to, i.e just the same place for MMI, or the 'real' num accs for MEE; in
         both cases these are located in hset, acc index 1.
      */
      /* MLE accs  */
      ma_src = ((MuAcc*)GetHook(mp->mean)) + (THREEACCS ? 2 : 0); /*Not used if OldParmsForICrit==TRUE. */
      va_src = ((VaAcc*)GetHook(mp->cov.var)) + (THREEACCS ? 2 : 0);
      /* Num accs [may be different]*/
      ma_dst = GetHook(mp->mean);
      va_dst = GetHook(mp->cov.var);
    
      {
         if(ma_src->occ>0.001){
            s= Tau/ma_src->occ;
	
            for(i=1;i<=VectorSize(mp->mean);i++){
               float srcmu =ma_src->mu[i], srcvar = va_src->cov.var[i];
               ma_dst->mu[i] += s*srcmu;
               va_dst->cov.var[i] += s*srcvar;

            }
            ma_dst->occ += s*ma_src->occ;
            va_dst->occ += s*va_src->occ;
         }
      }
   }
   EndHMMScan(&hss);
}


void AddPriorsFromPriorHMM(int dst_index, float Tau, float K, Boolean IsMMI, float ISmoothTau){
   /* Normally both Mu and Var occupancy will be identical. */
   HMMScanState hss, hss_prior;
   NewHMMScan(&hset,&hss); NewHMMScan(&hset_prior,&hss_prior); 

   while(GoNextMix(&hss,FALSE) && GoNextMix(&hss_prior,FALSE)){
      int i;
      MixPDF *mp;
      MuAcc *ma_dst; VaAcc  *va_dst;
      mp = hss.mp;

      /* _src refers to the place the MLE accs are stored [ hset,index 0 for MMI, hset,index 1 for MEE ] 
         _dst refers to the place we add the MLE accs to, i.e just the same place for MMI, or the 'real' num accs for MEE; in
         both cases these are located in hset, acc index 1.
      */

      ma_dst = ((MuAcc*)GetHook(mp->mean)) + dst_index; 
      va_dst = ((VaAcc*)GetHook(mp->cov.var)) + dst_index; 
      {
        
        if(hss_prior.mp->ckind!=DIAGC) HError(1, "Wrong ckind in prior HMMSet.");
        for(i=1;i<=VectorSize(mp->mean);i++){
          float srcmu = hss_prior.mp->mean[i], 
            srcvar = srcmu*srcmu + hss_prior.mp->cov.var[i];         
          float priormu,priorvar;

          priormu  = K*srcmu   +  (1-K)  *  ( Tau * srcmu + ma_dst->mu[i] ) / ( Tau + ma_dst->occ );
          priorvar = K*srcvar   +  (1-K)  *  ( Tau * srcvar + va_dst->cov.var[i] ) / ( Tau + va_dst->occ );

          if(!IsMMI){
             ma_dst->mu[i] = 100 * priormu; /* doesnt actually matter, 100 is the occupancy I randomly choose, should be >MinOcc though. */
             va_dst->cov.var[i] = 100 * priorvar; /* doesnt actually matter, 100 is the occupancy I randomly choose, should be >MinOcc though. */
          } else {
             ma_dst->mu[i] +=  ISmoothTau *  priormu; /* here ma_dst is the numerator *and* MLE occs... */
             va_dst->cov.var[i] += ISmoothTau * priorvar;
          }
        }
        if(!IsMMI){
           ma_dst->occ  = 100; /* This is an arbitrary value which does not affect results,*not* a tau value. */
           va_dst->occ  = 100;
        }else{
          ma_dst->occ += ISmoothTau; va_dst->occ += ISmoothTau;
        }
      }
   }
   EndHMMScan(&hss);    EndHMMScan(&hss_prior);
}

static void SmoothWeightsFromPriorHMM(int index, float Tau){
   HMMScanState hss,hss_prior;
   NewHMMScan(&hset,&hss);  NewHMMScan(&hset_prior,&hss_prior);
   while(GoNextStream(&hss,FALSE) && GoNextStream(&hss_prior,FALSE)){
      WtAcc *wa_dst; int m,M; 
      M = hss.M; 
      wa_dst = ((WtAcc*)hss.ste->hook) + index;
      for(m=1;m<=M;m++) wa_dst->c[m] += Tau * hss_prior.ste->spdf.cpdf[m].weight;
      wa_dst->occ += Tau;
   }
   EndHMMScan(&hss); EndHMMScan(&hss_prior); 
}


static void SmoothTransFromPriorHMM(int index, float Tau){
   HMMScanState hss,hss_prior;
   NewHMMScan(&hset,&hss);    NewHMMScan(&hset_prior,&hss_prior); 
   do{
      TrAcc *ta_dst; int m,M;
      ta_dst = ((TrAcc*)GetHook(hss.hmm->transP)) + index; 
      M=hss.hmm->numStates;
      for(m=1;m<M;m++){
         int n;
         
         for(n=1;n<M;n++)   ta_dst->tran[m][n] += Tau * hss_prior.hmm->transP[m][n];
         ta_dst->occ[m] += Tau;
      }
   } while(GoNextHMM(&hss) && GoNextHMM(&hss_prior));
   EndHMMScan(&hss); EndHMMScan(&hss_prior);
}




static void FixHMMForICrit(){
   Boolean ISmoothingDone=FALSE;

   if(PriorTau>0||PriorK>0||PriorK>0||PriorTauTrans>0) {
     if(!hset_prior_initialised)  HError(-1, "Config indicates that you intend to use a prior model (-Hprior), but none supplied.");
   } else {
     if(hset_prior_initialised)  HError(1, "Config indicates that you are not making use of the prior model (-Hprior), which has been supplied.");
   }

   if(hset_prior_initialised){  /* Using a prior HMM set */
     if(hset.ckind == FULLC) HError(1, "Prior HMM set not supported with FULLC.");
     if(PriorTauWeights) SmoothWeightsFromPriorHMM(THREEACCS ? 2 : 0, PriorTauWeights); 
     if(PriorTauTrans) SmoothTransFromPriorHMM(THREEACCS ? 2 : 0, PriorTauTrans); 
     if (THREEACCS){ /* if MPE... */
       AddPriorsFromPriorHMM(2 /* to index-pos 2 */, PriorTau, PriorK, FALSE,0);          /* The ML accs are in position 2 for MPE. */
     } else if (ML_MODE){ /* MLE */
       AddPriorsFromPriorHMM(0 /* to index-pos 0 */, PriorTau, PriorK, FALSE,0);          /* The ML accs are in position 0 for MMI. */
     } else { /* MMI --> do I-smoothing and MAP in one go. */
       AddPriorsFromPriorHMM(0 /* to index-pos 0 */, PriorTau, PriorK, TRUE, ISmoothTau); /* Combines I-smoothing with estimating the "center of the prior" in a MAP fashion. */
       ISmoothingDone=TRUE;
     }

   }

   if(ISmoothTau>0 && !ISmoothingDone && hset.ckind != FULLC) _FixHMMForICrit(ISmoothTau, THREEACCS);
  
   if(ISmoothTauWeights>0)
      FixWeightsForICrit(ISmoothTauWeights, THREEACCS);

   if(ISmoothTauTrans>0)
      FixTransForICrit(ISmoothTauTrans, THREEACCS);
}




void UpdateWeightsOrTrans(int M, float *acc1, float *acc2, float *mixes, float *oldMixes, float C){ 
   int iter=0;
   int m;
   float objective = 0, last_objective=0, last_last_objective=0;
   float fmax,csum;  

   if(C == 0) return; /*Don't update!  This is really c=infty.*/

   { 
      float sum=0; for(m=1;m<=M;m++) sum += acc1[m]; 
      if(sum < 1) return;
   }/*the original objective function (6.18) doesn't have a closed-form solution to optimise it, then iterative method is used*/
   while(1){

      iter++;
    
      /*calc. objective. */
      last_last_objective = last_objective; last_objective = objective;  objective=0;
      for(m=1;m<=M;m++)
         if(mixes[m]>0)/*check 6.19 pp82@Daniel thesis*/
            objective += (acc1[m]*log(mixes[m]) - ( (acc2?acc2[m]:0) / C * exp(log(mixes[m]/oldMixes[m])*C)));


      /*if(fabs(objective-last_objective) < 1.0e-8*(fabs(objective)+fabs(last_objective)) || iter > 100) break; */
      if(objective < last_objective && objective < last_last_objective && !(last_last_objective==0) 
         && fabs(objective-last_objective) > fabs(objective)*0.0001 )
         HError(-1, "Objective not increasing: %f<%f,<%f", objective, last_objective, last_last_objective);
      if(iter>100) break; /*this seems to work fine.*/
    
      /*find max f_m, check 6.21 pp83@Daniel thesis*/
      fmax = 0;
      for(m=1;m<=M;m++){
         if(mixes[m]>0){
            float f = (acc2?acc2[m]:0)/(oldMixes[m]) * exp(log(mixes[m]/oldMixes[m])*(C-1));
            if(C>1)
               f *= C;
            fmax = MAX(fmax,  f);
         }
      }
      csum=0;/*check 6.20's denominator pp83@Daniel thesis*/
      for(m=1;m<=M;m++){
         if(mixes[m]>0){
            float f = (acc2?acc2[m]:0)/(oldMixes[m]) * exp(log(mixes[m]/oldMixes[m])*(C-1));
            mixes[m] = (mixes[m]*(fmax-f)) + acc1[m];
            csum+=mixes[m];
         }
      }/*check 6.20 pp83@Daniel thesis*/
      for(m=1;m<=M;m++)
         mixes[m] /= csum;
   }
   return;
}

void UpdateWeight(int s, StreamElem *ste){
   int i,n,M=0;
   WtAcc *wa1,*wa2,*wa3;
   wa1 = (WtAcc *)ste->hook;
   wa2 = (ML_MODE?NULL:wa1+1);
   wa3 = (THREEACCS?wa1+2:NULL); /*non-NULL in MPE case, where it is the ML accs. */

   if(!wa1) return; /* Already done, so set to NULL */

   switch (hsKind){
   case PLAINHS:
   case SHAREDHS:
      M=ABS(ste->nMix);
      break;
   case TIEDHS:
      M = hset.tmRecs[s].nMix;
   default: HError(1, "Unhandled hsKind.");
   }

   if( (wa3?wa3->occ:wa1->occ) > MinOccWeights ){ /* more than MinOccWeights ML stats */
      Vector NewWghts = CreateVector(&gstack, M);
      Vector OldWghts = CreateVector(&gstack, M);
    
      switch(hsKind){
      case PLAINHS: case SHAREDHS:
         for(n=1;n<=M;n++)	NewWghts[n] = OldWghts[n] = ste->spdf.cpdf[n].weight;
         break;
      case TIEDHS:
         for(n=1;n<=M;n++)	NewWghts[n] = OldWghts[n] = ste->spdf.tpdf[n];
         break;
      default: HError(1, "Unhandled hsKind.");
      }
      if(uFlagsMLE & UPMIXES && wa2) for(n=1;n<=M;n++) wa2->c[n] = 0.0;

      UpdateWeightsOrTrans(M, wa1->c, wa2?wa2->c:NULL, NewWghts, OldWghts, CWeights);
      for(i=1;i<=M;i++)
        if(NewWghts[i] == 0 && OldWghts[i] != 0)
          HError(-1, "Weights going to zero: advise setting e.g. ISMOOTHTAUW = 10 ");

      switch (hsKind){
      case PLAINHS:
      case SHAREDHS:
         for(n=1;n<=M;n++){
            ste->spdf.cpdf[n].weight=(NewWghts[n] > MINMIX ? NewWghts[n] : 0.0);
         }
         break;
      case TIEDHS:
         for(n=1;n<=M;n++){
            ste->spdf.tpdf[n]=(NewWghts[n] > MINMIX ? NewWghts[n] : 0.0);
         }
         break;
      default: HError(1, "Unhandled hsKind.");
      }
      Dispose(&gstack, NewWghts); /*disposes of both.*/
    
      if (mixWeightFloor>0.0){
         switch (hsKind){
         case PLAINHS:
         case SHAREDHS:
            FloorMixes(ste->spdf.cpdf+1,M,mixWeightFloor); 
            break;
         case TIEDHS:
            FloorTMMixes(ste->spdf.tpdf,M,mixWeightFloor);
            break;
         default: HError(1, "Unhandled hsKind.");
         }
      }
   }
   ste->hook = 0;
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
     int M=0,m=0,floored=0,equal=0;
     vsize = hset1->swidth[s];
     
     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
       M++;
     }
     EndHMMScan(&hss1); 

     varray = New(&gstack, sizeof(float*) * (vsize+1));
     for(i=1;i<=vsize;i++) varray[i] = New(&gstack, sizeof(float) * M);

     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
       int k;
       if(hss1.mp->ckind != DIAGC ) HError(1, "FloorVars expects DIAGC covariances. ");
       
       for(k=1;k<=vsize;k++){
          varray[k][m] = hss1.mp->cov.var[k];
       }
       m++;
     }
     EndHMMScan(&hss1); 
     
     for(i=1;i<=vsize;i++){
        qsort((char *) varray[i], M, sizeof(float), fltcompare);
     }
     m=0;
     if(varFloorPercent <=0 || varFloorPercent >= 100) HError(1, "varFloorPercent should be <100 and >0..");
     

     NewHMMScan(hset1,&hss1); 
     while(GoNextMix(&hss1,FALSE)){
        int k, Pos = (int)(varFloorPercent*0.01*M);
        for(k=1;k<=vsize;k++){
           if(hss1.mp->cov.var[k] < varray[k][Pos]){
              hss1.mp->cov.var[k] =  varray[k][Pos];
              floored++;
           } else  if(hss1.mp->cov.var[k] == varray[k][Pos])
              equal++;
        }
     }
     EndHMMScan(&hss1); 

     printf("Floored %d (expected to floor %d), %d were equal to floor\n", floored, 
            (int)( varFloorPercent * 0.01 * M * vsize), equal);
  }
  FixAllGConsts(hset1);
}
/*update all parameters of regression functions*/
void UpdateFproj(){
	HMMScanState hss;
	int numDev=0;

	DAcc* dAcc=(DAcc*)hset.dAcc;
	Fprojection* fproj=hset.fproj;

	int targetDim=fproj->targetDim;
	int nparm=fproj->inputDimM;
	int i,j;
	char newFn[MAXSTRLEN];
	Matrix proj;
	float stepsize,norm,devi,cij,pn,smooth,p,n;

	Vector gradient=CreateVector(&gstack,nparm);
	Vector avgdev=CreateVector(&gstack,targetDim);
	ZeroVector(avgdev);
	NewHMMScan(&hset, &hss);
	while (GoNextMix(&hss,FALSE)) {
		Vector variance=hss.mp->cov.var;
		for(i=1;i<=targetDim;i++){
			float var=hss.mp->ckind==DIAGC?variance[i]:1/variance[i];
			avgdev[i]+=sqrt(var);
		}
		numDev++;
	}
	for(i=1;i<=targetDim;i++){
		avgdev[i]/=numDev;
	}
	EndHMMScan(&hss);


	proj=fproj->xform;
	/*do the update the parameters of regression functions*/
	for(i=1;i<=targetDim;i++){/*every function has a gradient and a Hessian*/
		float chk11=fabs(dAcc->check1[i][1]);
		float chk12=fabs(dAcc->check1[i][2]);
		float chk21=fabs(dAcc->check2[i][1]);
		float chk22=fabs(dAcc->check2[i][2]);
		float chk1=dAcc->check1[i][1]+dAcc->check1[i][2];
		float chk2=dAcc->check2[i][1]+dAcc->check2[i][2];
		printf("check1(d=%e, ind=%e, sum=%e, ratio_d=%f, ratio_ind=%f\n",chk11, chk12, chk1, chk1/chk11, chk1/chk12);
		printf("check2(d=%e, ind=%e, sum=%e, ratio_d=%f, ratio_ind=%f\n",chk21, chk22, chk2, chk2/chk21, chk2/chk22);
	}
	for(i=1;i<=targetDim;i++){/*every function has a gradient and a Hessian*/
		for(j=1;j<=nparm;j++){
			gradient[j]=dAcc->pGradient[i][j]-dAcc->nGradient[i][j];
		}
		norm=P_norm(gradient,pnorm);
		/*float devi=sqrt(vFloor[1][i]);*/
		devi=avgdev[i];
		printf("Function %d, P-Norm(%d) of Gradient=%f\n",i,pnorm,norm);
		/*steepest ascent update*/
		for(j=1;j<=nparm;j++){
			if(gradient[j]==0){
				printf("no update for [%d][%d]=%e\n",i,j,proj[i][j]);
				continue;
			}
			dAcc->sumSq[i][j]*=dAcc->sumSq[i][j];
			cij=dAcc->sumSq[i][j]/dAcc->sqSum[i][j];
			pn=dAcc->pGradient[i][j]+dAcc->nGradient[i][j];
			smooth=0.5*Tau*pn/cij;
			p=dAcc->pGradient[i][j]+smooth;
			n=dAcc->nGradient[i][j]+smooth;
			stepsize=devi/(RE*(p+n));

			printf("old[%d][%d]=%e,",i,j,proj[i][j]);
			proj[i][j]=proj[i][j]+stepsize*gradient[j];
			printf("new[%d][%d]=%e,",i,j,proj[i][j]);
			printf("stepsize=%e, direction=%e, pos=%e, neg=%e, smooth=%e\n",stepsize,gradient[j],p,n,smooth);
		}
		printf("\n");
	}


	MakeFN(fproj->xformName,xfInfo.outXFormDir,xfInfo.outXFormExt,newFn);
	SaveFprojection(&hset, fproj, newFn, saveBinary);
	if (trace&T_TOP) {
		if (nFloorVar > 0)
			printf("Total %d (%.2f%%) floored variance elements in %d (%.2f%%) different mixes\n",
				nFloorVar, nFloorVar*100.0/nVar  ,nFloorVarMix, nFloorVarMix*100.0/nMix);
		if(nFloorWeight > 0) printf("Total %d (%.2f%%) floored weights\n", nFloorWeight,nFloorWeight*100.0/
								  nWeight);
		if (mmfFn == NULL)
			printf("Saving hmm's to dir %s\n",(newDir==NULL)?"Current":newDir);
		else
			printf("Saving hmm's to MMF %s\n",mmfFn);
		fflush(stdout);
	}
	ClearSeenFlags(&hset,CLR_ALL); /*?*/
	SaveHMMSet(&hset,newDir,newExt,NULL,saveBinary);
	if (trace&T_TOP) {
		printf("Saved hmm\'s with updated model parameters\n");
		PrintCriteria();
	}
}




/* UpdateModels: update all models and save them in newDir if set,
   new files have newExt if set */
void UpdateModels(void)
{
   HMMScanState hss;
   /*Measure total occupancy and update Gaussians. */
   double mlocc=0,numocc=0,denocc=0; int nMix=0,nLeft=0;
   MuAcc *ma1,*ma2, *ma3;       VaAcc *va1,*va2,*va3;

   if (parMode == -1){
      ConvDiagC(&hset,TRUE);/*invDiag -> Diag*/
      ConvExpWt(&hset);/*logwt -> wt*/
   }
   if (uFlags&UPFPROJ){
   	  UpdateFproj();
   	  return;
   }
   if (uFlags&UPSEMIT) {
	   if(xfInfo.outXForm->akind==TREE){
		   HError(1,"Tree based update is not supported");
	   }
      UpdateSemiTiedModels(&hset, &xfInfo);
      return;
   }
   if(hcrit != 1.0){
	   ScaleAccsParallel(&hset,hcrit,1);
   }

   if((uFlags&UPMEANS) || (uFlags&UPVARS)){
   FixHMMForICrit(); /*alter hset to make i-crit.*/
   RestoreAccsParallel(&hset, 0);  /*Removes mu offsets, 0 is index of numerator */
   if(!ML_MODE) RestoreAccsParallel(&hset, 1);/*Removes va offsets, 1 is index of denominator */
   if(THREEACCS) RestoreAccsParallel(&hset, 2);
   if(MMIPrior) RestoreAccsParallel(&hset, 3);

   if(hset.hsKind==PLAINHS || hset.hsKind==SHAREDHS){
      NewHMMScan(&hset, &hss);
      while(GoNextMix(&hss, FALSE)){ /* Update mixes. */
         ma1 = (MuAcc*)GetHook(hss.mp->mean);
         ma2 = (ML_MODE?NULL:ma1+1);
         ma3 = (THREEACCS?ma1+2:NULL);

         va1 = (VaAcc*)GetHook(hss.mp->cov.var);
         va2 = (ML_MODE?NULL:va1+1);
         va3 = (THREEACCS?va1+2:NULL);

         numocc+=ma1->occ;
         if(ma2) denocc+=ma2->occ;
         if(ma3) mlocc+=ma3->occ;

         if(UpdateGauss(hss.s, hss.mp)) nMix++; else nLeft++;
      }
      EndHMMScan(&hss);
   } else if(hset.hsKind==TIEDHS){
      int s;
      for(s=1;s<=S;s++){
         int m,M = hset.tmRecs[s].nMix;	          /* components */
         for(m=1;m<=M;m++){
            ma1 = (MuAcc*)GetHook(hss.mp->mean);
            ma2 = (ML_MODE?NULL:ma1+1);
            ma3 = (THREEACCS?ma1+2:NULL);

            va1 = (VaAcc*)GetHook(hss.mp->cov.var);
            va2 = (ML_MODE?NULL:va1+1);
            va3 = (THREEACCS?va1+2:NULL);

            numocc+=ma1->occ;
            if(ma2) denocc+=ma2->occ;
            if(ma3) mlocc+=ma3->occ;

            if(UpdateGauss(s, hss.mp)) nMix++; else nLeft++;
         }
      }
   } else HError(1, "Unknown hsetkind.");
   printf("Num occ=%f,Den occ=%f\n",numocc,denocc);
   }
   if(uFlags&UPTRANS){
	   UpdateTrans();
   }
   if(uFlags&UPRWGTS){
	   UpdateRWeights();
	   if(hset.tiedTVWR){
		   UpdateWeights();
	   }
   }
   if(uFlags&UPMIXES){
	   UpdateWeights();
   }

   if(varFloorPercent){
      int s;
      printf("Flooring all vars to the %f'th percentile of distribution... ", varFloorPercent);
      for(s=1;s<=hset.swidth[0];s++)
         FloorVars(&hset,s);
   }


   FixAllGConsts(&hset);

   if (trace&T_TOP) {
      if (nFloorVar > 0)
         printf("Total %d (%.2f%%) floored variance elements in %d (%.2f%%) different mixes\n",
                nFloorVar, nFloorVar*100.0/nVar  ,nFloorVarMix, nFloorVarMix*100.0/nMix);
      if(nFloorWeight > 0) printf("Total %d (%.2f%%) floored weights\n", nFloorWeight,nFloorWeight*100.0/
                                  nWeight);
      if (mmfFn == NULL)
         printf("Saving hmm's to dir %s\n",(newDir==NULL)?"Current":newDir); 
      else
         printf("Saving hmm's to MMF %s\n",mmfFn);
      fflush(stdout);
   }
   ClearSeenFlags(&hset,CLR_ALL); /*?*/
  
   SaveHMMSet(&hset,newDir,newExt,NULL,saveBinary);
   if (trace&T_TOP) {
      printf("Saved hmm\'s with updated model parameters\n");
      PrintCriteria();
   } 
}





/* ----------------------------------------------------------- */
/*                      END:  HMMIRest.c                       */
/* ----------------------------------------------------------- */
