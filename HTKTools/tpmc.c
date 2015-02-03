char *tpmc_version = "!HVER!TPMC:   3.4.1 [CUED 12/22/11]";
char *tpmc_vc_id = "$Id: TPMC.c,v 1.1.1.1 2006/12/22 09:55:01 Exp $";

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
#include "exmath.h"
/* Parameters Global storage  */
static MemHeap hmmNoiseStack;
static HMMSet hsetNoise;
static char *newDir = "";


static Boolean keepClean=FALSE;/*clean clean statistics if gmm noises model applied*/
static Boolean forceTVWR=0;
static Vector vFloor[SMAX]; /* variance floor - default is all zero */
static float minVar  = 0.0;      /* minimum variance (diagonal only) */
/* Prototype function*/
DMatrix * pmc(DMatrix *preNoiseProDMatrix,DVector cleanMeanMfcc,DVector cleanVarMfcc,double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix,DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,int numComponents,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix* postprocess(DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanMeanLog,DMatrix cleanVarLog,double alpha,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,DMatrix dctDMatrix,DMatrix cepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix * invTrajectoryHMM(DMatrix trajMeanVect,DMatrix trajVarMat,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int D,int numObsVectors,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix mfcc(DMatrix dataLog,DMatrix dctDMatrix,DMatrix cepDMatrix,int type,char* covarOpt,char* option);
DMatrix computeMeanAvg(int numObsVectors,int numFeatures,DMatrix cHat);
DMatrix computeVarAvg(int numObsVectors,int numFeatures,DMatrix PHat,DMatrix factDMatrix);
DMatrix* linear2log(DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanVarLog,double alpha,char* covarOpt,char* option);
DMatrix* preprocess(DVector noiseMeanMfcc, DVector noiseVarMfcc, int M, DMatrix Sq, int numDctDVectors, int D,DMatrix W, int windowL, DMatrix invDctDMatrix, DMatrix invCepDMatrix, char *covarOpt, char* option,DMatrix selMat, int isDebug,int t);
DMatrix invMfcc(DMatrix dataMfcc, DMatrix invDctDMatrix,DMatrix invCepDMatrix, int type, char *covarOpt, char *option,int t);
DMatrix *trajectoryHMM(DVector meanMfcc, DVector varMfcc,DMatrix Sq, DMatrix W, int windowL, int M, char* option,DMatrix selMat, int isDebug, int t);
DMatrix selectDVectors(DVector data, DMatrix selMat);
DVector moveC0Top(DVector data,int M, int D);
DVector moveC0Down(DVector data,int M, int D);
void processOptions(char **option, char **wOpt, int *numPoints, int numObsVectors,int *windowL, int *numDctDVectors);

/* The following is HTK standard function*/
void Summary(void)
{
   printf("\nTPMC Command Summary\n\n");

   Exit(0);
}

void ReportUsage(void)
{
	printf("\nUSAGE: TPMC [options] Cleanpath Noisepath alpha ... \n\n");
	printf(" Option                                       Default\n\n");
	printf(" -mmfnoise  s path of noise model.\n");
	printf(" -hmmNoiseList  s path to list phones of noise.\n");
	printf(" -alpha  i the amount of noise used to corruped speechs.\n");
	printf(" -option  s trajectory PMC - simple average.\n");
	printf(" -numPoints  i if option ends with '_r', and numPoints > 0: generate sample points so 	that later HMM retrain could be done.\n");
	printf(" -noiseAddOpt  s this option is inactive so just set to \"normal\".\n");
	printf(" -covarOpt  s set to 'diag'.\n");
	printf(" -numObsVectors  i num of observation vectors (set to 2) to be.\n");
	printf(" -wOpt  s set 'l1-wNormal' reported in the paper (here, windowL=1).\n");
	printf(" -isDebug  i 0 (no output), 1 (some summary), 2 (full debug).\n");
	PrintStdOpts("Q");
	printf("\n\n");
}
/*LoadMMF: get mean and var in hset*/
void loadMMF(HMMSet *hset, Boolean ldBinary,DVector *stateNoiseMeanData,DVector *stateNoiseVarData, DVector stateNoiseWeights, char **arNoiseInfo)
{
	HLink hmm;
	HMMScanState hss;
	int numComponents = 0; int numStates = 0; 	int numHMMs = 0;
	int numFeatures = hset->vecSize;
	MixPDF* mp;
	DVector vTemp = NULL;
	TriMat fcvar=NULL;
	vTemp = CreateDVector(&gstack,numFeatures);
	fcvar = CreateTriMat(&gstack,numFeatures);
	int ii = 1;
	NewHMMScan(hset, &hss);
	do {
		numHMMs++;
		hmm = hss.hmm;
		int i,s,m,M,N,S;
		StateElem *se;
		StreamElem *ste;
		MixtureElem *me;
		N = hmm->numStates;
		se = hmm->svec+2;
		S = hset->swidth[0];
		for (i=2; i<N; i++,se++){
			ste = se->info->pdf+1;
			numStates++;
			for (s=1;s<=S;s++,ste++){
				me = ste->spdf.cpdf + 1; M = ste->nMix;
				for (m=1;m<=M;m++,me++){
					mp=me->mpdf;
					numComponents++;
					stateNoiseWeights[numComponents]=me->weight;
					for(ii=1; ii<=numFeatures; ii++)
						vTemp[ii] = mp->mean[ii];
					stateNoiseMeanData[numComponents] = CreateDVector(&gstack,numFeatures);
					CopyDVector(vTemp,stateNoiseMeanData[numComponents]);
					if(mp->ckind==DIAGC){
						TouchV(mp->cov.var);
						for(ii=1; ii<=numFeatures; ii++){
							vTemp[ii] = mp->cov.var[ii];
						}
					}
					else if(mp->ckind==FULLC){
						TouchV(mp->cov.inv);
						CovInvert(&gstack, mp->cov.inv,fcvar);
						for(ii=1; ii<=numFeatures; ii++){
							/*vTemp[ii] = 1/hss.mp->cov.inv[ii][ii];*/
							vTemp[ii]=fcvar[ii][ii];
						}
					}
					else{
						HError(1,"Unsupported noise model");
					}
					stateNoiseVarData[numComponents] = CreateDVector(&gstack,numFeatures);
					CopyDVector(vTemp,stateNoiseVarData[numComponents]);
				}
				char* phone = malloc(15*sizeof(char));
				char *strMixture = malloc(5*sizeof(char));;
				sprintf(phone, "h-%s-s%d", hss.mac->id->name, i);
				sprintf(strMixture, "%d", M);
				arNoiseInfo[3*numStates] = phone;
				arNoiseInfo[3*numStates+1] = strMixture;
				arNoiseInfo[3*numStates+2] = strMixture;
			}
		}
	} while (GoNextHMM(&hss));
	EndHMMScan(&hss);
	/* load mean and var from tied*/

}
static HMMSet hset;                  /* Set of HMMs to be re-estimated */
static MemHeap hmmStack;           /* HMM defs and related structures */

int main(int argc, char *argv[])
{
	char *cleaninfo,*noisePrefix,*hmmNoiseList;
	int numPoints, numObsVectors, isDebug;
	double alpha;
	char *option, *noiseAddOpt, *covarOpt, *wOpt;
	DMatrix dctDMatrix,cepDMatrix,invDctDMatrix, invCepDMatrix, Sq, W, VS;
	char*s;
	setbuf(stdout, NULL); /*unbuffered output.*/
	/*InitShell(argc,argv,tpmc_version,tpmc_vc_id);*/
	if(InitShell(argc,argv,tpmc_version,tpmc_vc_id)<SUCCESS)
	      HError(2300,"HERest: InitShell failed");
	InitMem();
	InitMath();
	InitSigP();
	InitAudio();
	InitWave();
	InitVQ();
	InitModel();
	InitParm();
	InitLabel();
	InitTrain();
	InitUtil();
	InitFBLat();
	InitArc();
	InitDict();
	InitLat();
	InitNet();
	if (!InfoPrinted() && NumArgs() == 0)
	      ReportUsage();
	if (NumArgs() == 0) Exit(0);

	CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
	CreateHMMSet(&hset,&hmmStack,TRUE);
	CreateHeap(&hmmNoiseStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
	CreateHMMSet(&hsetNoise,&hmmNoiseStack,TRUE);

	while (NextArg() == SWITCHARG) {
		s = GetSwtArg();
			switch(s[0]){
				case 'a':
					if(!strcmp(s, "alpha")){
						alpha = GetChkedFlt(0,10,s);
					}else HError(1, "Unknown option %s",s);
					break;
				case 'o':
					if(!strcmp(s, "option")) {
						option=GetStrArg();
					}
					break;
				case 'c':
					if(!strcmp(s,"covarOpt")){
						covarOpt = GetStrArg();
					}
					break;
				case 'h':
					if(!strcmp(s,"hmmNoiseList")){
						hmmNoiseList = GetStrArg();
					}
					break;
				case 'H':
					if (NextArg() != STRINGARG)
						HError(2319,"tpmc: HMM macro file name expected");
					if(!strcmp(s,"H")){
						char *x = GetStrArg();
						AddMMF(&hset,x);
					}
					break;
				case 'i':
					if (!strcmp(s, "isDebug")) {
						isDebug =  GetChkedInt(0,1,s);
					}
					break;
				case 'n':
					if (!strcmp(s, "noiseAddOpt")) {
						noiseAddOpt = GetStrArg();
					}else if(!strcmp(s, "numPoints")){
						numPoints = GetChkedInt(0,100,s);
					}else if(!strcmp(s, "numObsVectors")){
						numObsVectors = GetChkedInt(0,100,s);
					} else HError(1, "Unknown option %s", s);
					break;
				case 'm':
					if(!strcmp(s,"mmfnoise")){
						char *x = GetStrArg();
						AddMMF(&hsetNoise,x);
					}
					break;
				case 's':
					if(!strcmp(s,"savePath")){
						newDir = GetStrArg();
					}
					break;
				case 'w':
					if(!strcmp(s,"wOpt")){
						wOpt = GetStrArg();
					}
					break;
			      case 'k':
			    	 keepClean=TRUE;
			    	 break;
			      case 'f':
			    	  if(!strcmp(s,"forceTVWR")){
			    		  forceTVWR=GetChkedInt(0,100,s);
			    	  }
			    	  else{
			    		  HError(1,"Unknown option %s",s);
			    	  }
			    	  break;
				default:
					HError(2319,"tpmc: Unknown switch %s",s);
					break;
			}
		} /*matches while(NextArg() == SWITCHARG)*/

	char *hmmList = GetStrArg();
	MakeHMMSet( &hset, hmmList ,NULL);
	LoadHMMSet( &hset,NULL,NULL);
	SetParmHMMSet(&hset);
	MakeHMMSet( &hsetNoise, hmmNoiseList ,NULL);
	LoadHMMSet( &hsetNoise,NULL,NULL);
	SetParmHMMSet(&hsetNoise);
	if(forceTVWR){
		hset.TVWR=TRUE;
		hset.regwSize=forceTVWR;
	}
	SetVFloor( &hset, vFloor, minVar);
	Vector minV=NULL;
	int vSize = hsetNoise.vecSize;
	int numStatesNoise = hsetNoise.numStates;
	int numMixNoise = hsetNoise.numMix;
	int M = numCepCoef + 1; /* num statics features*/
	int numPhones = hset.numPhyHMM;
	int totalNumFeatures = vSize;
	int D = totalNumFeatures / M; /* 1,...,(D-1)-order dynamic features*/
	int numDctDVectors;
	int windowL;
	int noiseMix;

	DVector *stateNoiseMeanData = malloc((numMixNoise+1)*sizeof(DVector));
	DVector *stateNoiseVarData = malloc((numMixNoise+1)*sizeof(DVector));
	DVector stateNoiseWeights = CreateDVector(&gstack,numMixNoise);
	char **arNoiseInfo = malloc((3*numStatesNoise+3)*sizeof(char*));

	/*typical setting, numObsVectors=2, windowL=1, numDctDVectors=numObsVectors+4*windowL=6*/
	processOptions(&option, &wOpt,&numPoints, numObsVectors,&windowL,&numDctDVectors);

	DMatrix *Matrices,selMat;
	/*typical setting, numChans=26, numCeps=12, cepLifter=22, M=12+1, D=3(39/13),*/
	Matrices = initializeMatrices(numObsVectors, numDctDVectors, numChans, numCepCoef, cepLifter, M, D, windowL,option, wOpt);
	dctDMatrix = Matrices[1];
	cepDMatrix = Matrices[2];
	invDctDMatrix = Matrices[3];
	invCepDMatrix = Matrices[4];
	Sq = Matrices[5];
	W = Matrices[6];
	VS = Matrices[7];
	selMat = ones(1, numObsVectors);/* % duplicate numObsVectors time for other options except tpmc. For tpmc, we will reset it later*/

	int i,j;
	DMatrix factDMatrix;/*=kron(ones(1, numObsVectors), eye(numFeatures,numFeatures));*/
	factDMatrix= CreateDMatrix(&gstack,M*D,M*D*numObsVectors);
	for(i=1;i<=M*D;i++)
		for(j=1;j<=M*D*numObsVectors;j++)
			if(i==j || M*D+i ==j )
				factDMatrix[i][j]=1;

	puts("Load noise statistics... ");
	DVector noiseMeanMfcc, noiseVarMfcc;
	loadMMF(&hsetNoise,0,stateNoiseMeanData,stateNoiseVarData,stateNoiseWeights,arNoiseInfo);
	noiseMeanMfcc = CreateDVector(&gstack,totalNumFeatures);
	noiseVarMfcc = CreateDVector(&gstack,totalNumFeatures);
	DMatrix ** preNoiseProDMatrix=New(&gstack,(numMixNoise)*sizeof(DMatrix*));
	for(noiseMix=1;noiseMix<=numMixNoise;noiseMix++){
		CopyDVector(stateNoiseMeanData[noiseMix],noiseMeanMfcc);
		CopyDVector(stateNoiseVarData[noiseMix],noiseVarMfcc);
		preNoiseProDMatrix[noiseMix] = preprocess(noiseMeanMfcc, noiseVarMfcc, M, Sq, numDctDVectors,
				D, W, windowL, invDctDMatrix, invCepDMatrix, covarOpt, option, selMat, isDebug,0);
	}

	DMatrix * noisyData;
	DVector cleanMeanMfcc,cleanVarMfcc;
	cleanMeanMfcc = CreateDVector(&gstack,totalNumFeatures);
	cleanVarMfcc = CreateDVector(&gstack,totalNumFeatures);
	puts("Start TPMC converting...");
	clock_t start = clock();
	HMMScanState hss;
	HLink hmm;
	int px;
	char type = 'q';
	char macName[255];
	char buf[255];
	DVector memMarker;
	NewHMMScan(&hset,&hss);
	float normwgt=1.0/(numMixNoise+1);
	if(keepClean){
		for (noiseMix=1;noiseMix<=numMixNoise;noiseMix++){
			stateNoiseWeights[noiseMix]*=(normwgt*numMixNoise);
		}
	}
	px=1;
	int stateId=0;
	do {
		hmm = hss.hmm;
		int i,s,m,j;
		StateElem *se;
		StreamElem *ste;
		MixtureElem *me;
		MixtureElem *newme;
		MixtureElem *noisyme;
		MixPDF* noisymp;
		se = hmm->svec+2;
		for (i=2; i<hmm->numStates; i++,se++){
			ste = se->info->pdf+1;
			if(IsSeen(se->info->nUse)){
				continue;
			}
			Touch(&se->info->nUse);
			for (s=1;s<=hset.swidth[0];s++,ste++){
				if(IsSeen(ste->nMix)){
					continue;
				}
				minV=vFloor[s];
				stateId++;
				/*printf("%d\n",stateId);*/
				if(numMixNoise>1){
					newme=CreateCME(&hset,ste->nMix*numMixNoise+(keepClean?ste->nMix:0));/*point to a new space*/
				}
				for(noiseMix=1;noiseMix<=numMixNoise;noiseMix++){
					me = ste->spdf.cpdf + 1;
					RegWgt* rwgt=NULL;
					if(forceTVWR){
						sprintf(buf,"Q_RW_%x_%x",stateId,noiseMix);
						LabId labid=GetLabId(buf,TRUE);
						rwgt=CreateRegWgt(&hset,s,-1);
						for(i=1;i<=hset.regwSize;i++){
							rwgt->regparm[i]=1.0/hset.regwSize;
						}
						IncUse(rwgt->regparm);
						MLink mac=NewMacro(&hset,0,type,labid,rwgt);
						SetMacroUse(mac,1);
					}
					for (m=1;m<=ste->nMix;m++,me++){
						for(j=1;j<=vSize;j++){
							cleanMeanMfcc[j]=me->mpdf->mean[j];
							cleanVarMfcc[j]=me->mpdf->cov.var[j];
						}
						noisymp=me->mpdf;
						noisyme=me;
						if(numMixNoise>1){
							noisyme=&newme[(noiseMix-1)*ste->nMix+m];
							noisymp=noisyme->mpdf = (MixPDF *)New(hset.hmem,sizeof(MixPDF));
							noisyme->weight=noisyme->sweight=me->weight*stateNoiseWeights[noiseMix];
							noisymp->nUse = 0; noisymp->hook = NULL; noisymp->gConst = LZERO;
							noisymp->mIdx = 0; noisymp->stream = 0; noisymp->vFloor = NULL; noisymp->info = NULL;noisymp->ckind=DIAGC;
							noisyme->rwgt=NULL;
							noisymp->mean = CreateSVector(hset.hmem,vSize);
							noisymp->cov.var=CreateSVector(hset.hmem,vSize);
							if(me->rwgt){
								noisyme->rwgt=CreateRegWgt(&hset,s,-1);
								noisyme->rwgt->nUse=me->rwgt->nUse;
								CopyVector(me->rwgt->regparm,noisyme->rwgt->regparm);
							}
							if(forceTVWR){
								noisyme->rwgt=rwgt;
							}
						}
						if(noisyme->weight<=MINMIX){/*Not only for efficiency but also for fixed mixture index*/
							noisyme->mpdf = EmptyMixPDF(&hset,hset.vecSize,s);
							continue;
						}
						memMarker=CreateDVector(&gstack,1);
						noisyData = pmc(preNoiseProDMatrix[noiseMix], cleanMeanMfcc, cleanVarMfcc,
								alpha, M, numObsVectors, numDctDVectors, D, W, VS, Sq, windowL, numChans,
								numCepCoef, dctDMatrix, cepDMatrix, invDctDMatrix, invCepDMatrix, covarOpt,
								option, 1, selMat,factDMatrix, isDebug); /*% numComponents = 1*/
						for(j=1;j<=vSize;j++){
							noisymp->mean[j]=noisyData[1][1][j];
		            		noisymp->cov.var[j]=noisyData[2][1][j]<minV[j]?minV[j]:noisyData[2][1][j];
		            	}
						FreeDVector(&gstack,memMarker);
					}
				}
				if(numMixNoise>1){
					if(keepClean){
						me=ste->spdf.cpdf + 1;
						RegWgt* rwgt=NULL;
						if(forceTVWR){
							sprintf(buf,"Q_RW_%x_%x",stateId,noiseMix);
							LabId labid=GetLabId(buf,TRUE);
							rwgt=CreateRegWgt(&hset,s,-1);
							for(i=1;i<=hset.regwSize;i++){
								rwgt->regparm[i]=1.0/hset.regwSize;
							}
							IncUse(rwgt->regparm);
							MLink mac=NewMacro(&hset,0,type,labid,rwgt);
							SetMacroUse(mac,1);
						}
						for (m=1;m<=ste->nMix;m++,me++){
							newme[ste->nMix*numMixNoise+m].weight=me->weight*normwgt;
							newme[ste->nMix*numMixNoise+m].sweight=me->weight*normwgt;
							newme[ste->nMix*numMixNoise+m].mpdf=me->mpdf;
							if(forceTVWR)
								newme[ste->nMix*numMixNoise+m].rwgt=rwgt;
						}
					}
					ste->spdf.cpdf=newme;
					ste->nMix=ste->nMix*numMixNoise+(keepClean?ste->nMix:0);
				}
				Touch(&ste->nMix);
			}
		}
		if(px%100==0){
			printf(".");fflush(stdout);
		}
		if(px%1000==0){
			printf("\nDone processing phone %d out of total %d phones...\n",px,numPhones);fflush(stdout);
		}
		px++;
	} while (GoNextHMM(&hss));
	EndHMMScan(&hss);
	FixAllGConsts(&hset);
	printf("\nDone after %f min\n", ((double)clock()-start)/(60*CLOCKS_PER_SEC));
	fflush(stdout);
	FreeDMatrix(&gstack,factDMatrix);
	printf("Save noisy MMF to %s\n",newDir);
	if(hset.xf)
		hset.xf->nUse=0;
	SaveHMMSet(&hset,newDir,NULL,NULL,0);
	return (0);          /* keep compiler happy */
}

/* return noisyMeanMfcc, noisyVarMfcc, trajMeanVect, trajVarMat*/
DMatrix * pmc(DMatrix *preNoiseProDMatrix,DVector cleanMeanMfcc,DVector cleanVarMfcc,double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix,DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,int numComponents,DMatrix selMat,DMatrix factDMatrix,int isDebug)
{
	DMatrix * noisyData;
	DMatrix  trajMeanVect, trajVarMat;
	DMatrix *preCleanProDMatrix;
	DMatrix cleanMeanLinear, cleanVarLinear, cleanMeanLog, cleanVarLog;
	DMatrix noiseMeanLinear, noiseVarLinear, noiseMeanLog, noiseVarLog;

	noiseMeanLinear = preNoiseProDMatrix[1];
	noiseVarLinear = preNoiseProDMatrix[2];
	noiseMeanLog = preNoiseProDMatrix[3];
	noiseVarLog = preNoiseProDMatrix[4];
	trajMeanVect = preNoiseProDMatrix[7];
	trajVarMat = preNoiseProDMatrix[8];

	/*obtain the linear and log domain statistics*/
	preCleanProDMatrix = preprocess(cleanMeanMfcc, cleanVarMfcc, M, Sq, numDctDVectors, D, W, windowL, invDctDMatrix, invCepDMatrix, covarOpt, option, selMat, isDebug,1);
	cleanMeanLinear = preCleanProDMatrix[1];
	cleanVarLinear = preCleanProDMatrix[2];
	cleanMeanLog = preCleanProDMatrix[3];
	cleanVarLog = preCleanProDMatrix[4];
	cleanMeanMfcc = preCleanProDMatrix[5][1];
	cleanVarMfcc = preCleanProDMatrix[6][1];
	trajMeanVect = preCleanProDMatrix[7];
	trajVarMat = preCleanProDMatrix[8];

	/*Combine model in the linear domain*/
	DMatrix noisyMeanLinear,noisyVarLinear;
	noisyMeanLinear = mulNumberwDMatrix(alpha,noiseMeanLinear);
	noisyMeanLinear = sumEDMatrix1(cleanMeanLinear,noisyMeanLinear);
	noisyVarLinear = mulNumberwDMatrix(alpha*alpha,noiseVarLinear);
	noisyVarLinear = sumEDMatrix1(cleanVarLinear,noisyVarLinear);

	/*Post process: linear -> ceps domain*/
	noisyData= postprocess(noisyMeanLinear, noisyVarLinear, noiseMeanLinear, noiseMeanLog,
			noiseVarLog, cleanMeanLinear, cleanMeanLog, cleanVarLog, alpha, cleanMeanMfcc,
			cleanVarMfcc, M, numObsVectors, numDctDVectors, D, W, VS, Sq, windowL, numChans,
			dctDMatrix, cepDMatrix, covarOpt, option, selMat,factDMatrix, isDebug);
	return noisyData;
}

/*%%% Post process: linear -> ceps domain %%%
 *  Return DMatrix with meanMfcc and varMfcc
 */

DMatrix* postprocess(DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanMeanLog,DMatrix cleanVarLog,double alpha,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,DMatrix dctDMatrix,DMatrix cepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug)
{
	 DMatrix  *noisyData = (DMatrix*)New(&gstack,3*sizeof(DMatrix));
	 DMatrix meanMfcc, varMfcc,meanLog, varLog;
	 DMatrix *logData =  (DMatrix*)New(&gstack,3*sizeof(DMatrix));;
	 DMatrix *dataMfcc = (DMatrix*)New(&gstack,3*sizeof(DMatrix));;


	 if (strcmp(option, "dpmc") == 0)
	 {
		printf("#! Option dpmc should not call postprocess method!\n");
		/*return noisyData;*/
	 }

	 /*% Convert linear -> log domain
	 */

	 logData = linear2log(meanLinear, varLinear, noiseMeanLinear, noiseVarLog, cleanMeanLinear, cleanVarLog, alpha, covarOpt, option);
	 meanLog = logData[1];
	 varLog = logData[2];

	 /*
	 % convert from log -> MFCC domain
	  */
	 meanMfcc = mfcc(meanLog, dctDMatrix, cepDMatrix,1 , covarOpt, option);

	 varMfcc = mfcc(varLog, dctDMatrix, cepDMatrix, 2, covarOpt, option);

	 /*% inverse trajectory model*/
	 if(strncmp(option, "traj", 4)==0)
	 {
		 dataMfcc = invTrajectoryHMM(meanMfcc, varMfcc, cleanMeanMfcc, cleanVarMfcc, M, D, numObsVectors, W, VS, Sq, windowL, option, selMat, factDMatrix, isDebug);
		 meanMfcc = dataMfcc[1];
		 varMfcc = dataMfcc[2];

	 }

	meanMfcc[1] = moveC0Down(meanMfcc[1], M, D);
	varMfcc[1] = moveC0Down(varMfcc[1], M, D);

	noisyData[1] = meanMfcc;
	noisyData[2] = varMfcc;

	return noisyData;
}
/*
%%%
% Input: mean row DVector, full cov
% Output: mean row DVectors, var row DVector
%%%
*/
DMatrix * invTrajectoryHMM(DMatrix trajMeanVect,DMatrix trajVarMat,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int D,int numObsVectors,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug)
{
	 DMatrix *dataMfcc = (DMatrix*)New(&gstack,3*sizeof(DMatrix));
	 DMatrix meanMfccSequence,varMfccSequence,meanMfcc,varMfcc;
	/*
	 %%% Average solution for varMfcc and meanMfcc %%%
	meanMfccSequence = W * trajMeanVect';
	varMfccSequence = W * trajVarMat * W';
	 */
	 meanMfccSequence = mulMatrix3Op6(W,trajMeanVect);


	varMfccSequence = mulMatrixOp11(W,trajVarMat);
	/* ShowDMatrix1("",varMfccSequence);*/
	 varMfccSequence = mulMatrix3Op9(varMfccSequence,W);

	 meanMfcc = computeMeanAvg(numObsVectors, M*D, meanMfccSequence);
	/* if strcmp(option, 'trajAvgNew') == 1 work on it later*/
	 varMfcc = diag(computeVarAvg(numObsVectors, M*D, varMfccSequence,factDMatrix));

	 if(strncmp(option, "tmpc", 4)!=0)
	 {
		 varMfcc = toRow(varMfcc);
		 meanMfcc = toRow(meanMfcc);

	 }
	 dataMfcc[1] = meanMfcc;
	 dataMfcc[2] = varMfcc;


	 return dataMfcc;
}

/*%%% Compute mean and var avg for trajAvg %%%*/
DMatrix computeMeanAvg(int numObsVectors,int numFeatures,DMatrix cHat)
{
	DMatrix m;
	/*int i;*/
	m= mean(reshape(cHat,numFeatures,numObsVectors),2);

	return m;
}

DMatrix computeVarAvg(int numObsVectors,int numFeatures,DMatrix PHat,DMatrix factDMatrix)
{
	DMatrix varMat;

	varMat = mulNumberwDMatrix(1/(double)numObsVectors,factDMatrix);

	varMat = mulMatrixOp6(varMat,PHat);

	varMat = mulMatrix3Op4(varMat,factDMatrix);

	return varMat;
}


/*
%%%
% mfcc
% Input: row mean DVector, var DMatrix
% Output: row mean DVector, row var (if diag) or full var
type=1:mean; type=2:var
%%%
*/
DMatrix mfcc(DMatrix dataLog,DMatrix dctDMatrix,DMatrix cepDMatrix,int type,char* covarOpt,char* option)
{
	DMatrix dataMfcc;
	if(type==1)/*mean*/
	{
		/*DCT transform dataMfcc = dctDMatrix * dataLog';*/

		dataMfcc = mulMatrix3Op7(dctDMatrix,dataLog);
		/*% cepstral liftering dataMfcc = cepDMatrix * dataMfcc;*/
		dataMfcc = mulMatrixOp5(cepDMatrix,dataMfcc);

	}else if(type==2)
	{
		/*% DCT transform dataMfcc = dctDMatrix * dataLog * dctDMatrix';
		*/
		dataMfcc = mulMatrixOp4(dctDMatrix,dataLog);

		/*dataMfcc = mulMatrix3Test(dataMfcc,dctDMatrix);*/
		dataMfcc = mulMatrix3Op3(dataMfcc,dctDMatrix);

		/*% cepstral liftering	dataMfcc = cepDMatrix * dataMfcc * cepDMatrix';
		 */
		dataMfcc = mulMatrixOp2(cepDMatrix,dataMfcc);
		dataMfcc = mulMatrix3Op(dataMfcc,cepDMatrix);

		if(strncmp(option, "traj", 4)!=0)
		 {
			dataMfcc = diag(dataMfcc);
		 }

	}

	dataMfcc = toRow(dataMfcc);


	return dataMfcc;
}
/*
%%%
% linear2log
% input: mean row DVector; var DMatrix
% output: mean row DVector; var: row DVector ('diag' case) or DMatrix ('full' case)
%%%
*/
DMatrix* linear2log(DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanVarLog,double alpha,char* covarOpt,char* option)
{
	DMatrix *logData = (DMatrix*)New(&gstack,3*sizeof(DMatrix));
	DMatrix meanLog, varLog;
	/*% convert from linear -> log domain*/
	/*  varLog = log((varLinear ./ (meanLinear'*meanLinear))+1);*/

	varLog = mulMatrix2Op1(meanLinear);

	varLog = divEDMatrix(varLinear,varLog);

	varLog = sumNumberwDMatrix(1,varLog);
	varLog = logDMatrix(varLog);

	/*%% fixing inf values for varLog (if any) %% work on it later
		[row, col, I] = find(isinf(varLog) == 1);
	*/

	/*
	%convert from linear -> log domain
	meanLog = log(meanLinear) - 0.5 * diag(varLog)';
	 */

	meanLog = mulNumberwDMatrix(0.5,toRow(diag(varLog)));
	meanLog = minusEDMatrix(logDMatrix(meanLinear),meanLog);
	logData[1] = meanLog;
	logData[2] = varLog;

	return logData;
}

/*meanLinear, varLinear, meanLog, varLog, newMeanMfcc, newVarMfcc, trajMeanVect, trajVarMat*/
DMatrix* preprocess(DVector noiseMeanMfcc, DVector noiseVarMfcc, int M, DMatrix Sq, int numDctDVectors, int D,DMatrix W, int windowL, DMatrix invDctDMatrix, DMatrix invCepDMatrix, char *covarOpt, char* option,DMatrix selMat, int isDebug, int t)
{

	DMatrix *preProDMatrix  = (DMatrix*)New(&gstack, 9*sizeof(DMatrix)) ;
	DMatrix meanLinear, varLinear, meanLog, varLog, newMeanMfcc, newVarMfcc, trajMeanVect, trajVarMat;
	DMatrix *trajHMM,temp;
	int i,n;
 /*[meanMfcc] = moveC0Top(meanMfcc, M, D);*/
	noiseMeanMfcc = moveC0Top(noiseMeanMfcc,M,D);

 /*[varMfcc] = moveC0Top(varMfcc, M, D);*/
	noiseVarMfcc = moveC0Top(noiseVarMfcc,M,D);
	n= M*D;
	newMeanMfcc = CreateDMatrix(&gstack,1,n);
	newVarMfcc = CreateDMatrix(&gstack,1,n);
	for(i=1; i<=n;i++){
		newMeanMfcc[1][i] = noiseMeanMfcc[i];
		newVarMfcc[1][i] = noiseVarMfcc[i];
	}
	preProDMatrix[5] = newMeanMfcc;
	preProDMatrix[6] = newVarMfcc;
/*
	if ~isempty(regexp(option, 'trajAvg') == 1) || ~isempty(regexp(option, 'trajNewAvg'))
		[meanMfcc, varMfcc] = trajectoryHMM(meanMfcc, varMfcc, Sq, W, windowL, M, option, selMat, isDebug);

		trajMeanVect = meanMfcc;
		trajVarMat = varMfcc;
*/

	if (strcmp(option,"trajAvg")==0){
		trajHMM = trajectoryHMM(noiseMeanMfcc, noiseVarMfcc, Sq, W, windowL, M, option, selMat, isDebug,t);
		trajMeanVect = trajHMM[1];
		trajVarMat = trajHMM[2];
 	}

	meanLog = invMfcc(trajMeanVect, invDctDMatrix, invCepDMatrix, 1, covarOpt, option,t);
	varLog = invMfcc(trajVarMat, invDctDMatrix, invCepDMatrix, 2, covarOpt, option,t);

 /*
 % cepstral -> log domain
 	[meanLog] = invMfcc(meanMfcc, invDctDMatrix, invCepDMatrix, 'mean', covarOpt, option);
  	[varLog] = invMfcc(varMfcc, invDctDMatrix, invCepDMatrix, 'var', covarOpt, option);
 */
 /*% convert from log -> linear domain*/
	temp = mulNumberwDMatrix(0.5,toRow(diag(varLog)));
	temp  = sumEDMatrix1(meanLog,temp);
	meanLinear = expDMatrix(temp);

 /*ShowDMatrix1("",meanLinear);*/

 /* % convert from log -> linear domain
 	meanLinear = exp(meanLog + diag(varLog)'/2);
 	varLinear = (meanLinear'*meanLinear) .* (exp(varLog) -1 );
*/
	varLinear = mulMatrix2Op1(meanLinear);
 /*ZeroDMatrix(temp);*/
	temp = expDMatrix(varLog);
	temp = sumNumberwDMatrix(-1,temp);
	varLinear = mulEDMatrix(varLinear,temp);

	preProDMatrix[1] = meanLinear;
	preProDMatrix[2] = varLinear;
	preProDMatrix[3] = meanLog;
	preProDMatrix[4] = varLog;
	preProDMatrix[7] = trajMeanVect;
	preProDMatrix[8] = trajVarMat;
 /*ShowDMatrix1("",varLinear);*/
 /*ShowDMatrix1("meanLinear",meanLinear);*/
	/*FreeDVector(&gstack,noiseMeanMfcc);*/

 return preProDMatrix;
}
/*
 * %%%
% invMfcc, return column DVector
% Input: row mean DVector, row var DVector / full DMatrix
% Output: row mean DVector, row var DVector / full DMatrix
type =1:mean,type =2:var
%%%
*/
DMatrix invMfcc(DMatrix dataMfcc, DMatrix invDctDMatrix,DMatrix invCepDMatrix, int type, char *covarOpt, char *option,int t)
{
	DMatrix dataLog;
	if(type==1)
	{
		dataMfcc = mulMatrix3Op8(invCepDMatrix,dataMfcc);

		if(t)
			{


				dataLog = mulMatrixOp8(invDctDMatrix, dataMfcc);
				/*ShowDMatrix1("",dataLog);*/

			}else
				dataLog = mulMatrix(invDctDMatrix, dataMfcc);

	}else if (type==2)
	{
		if(strncmp(option, "traj", 4)!=0)
		 {
			dataMfcc = diag(dataMfcc);
		 }
		dataMfcc = mulMatrixOp2(invCepDMatrix,dataMfcc);

		dataMfcc = mulMatrix3Op(dataMfcc,invCepDMatrix);
		/*dataLog = mulMatrixTest(invDctDMatrix,dataMfcc);*/
		dataLog = mulMatrixOp3(invDctDMatrix,dataMfcc);


		if(t)
		{

			/*dataLog = mulMatrix3(dataLog,invDctDMatrix);*/
			dataLog = mulMatrix3Op2(dataLog,invDctDMatrix);

		}
		else
			dataLog = mulMatrix3(dataLog,invDctDMatrix);
		/*dataMfcc = invCepDMatrix * dataMfcc * invCepDMatrix'; % inverse cepstral liftering*/
		  /*  dataLog = invDctDMatrix * dataMfcc * invDctDMatrix'; % inverse DCT*/
	}

	dataLog = toRow(dataLog);

	return dataLog;
}

/*
function [trajMeanVect, trajVarMat] = trajectoryHMM(meanMfcc, varMfcc, ...
																									Sq, W, windowL, M, option, selMat, isDebug)
*/
DMatrix *trajectoryHMM(DVector meanMfcc, DVector varMfcc,DMatrix Sq, DMatrix W, int windowL, int M, char* option,DMatrix selMat, int isDebug, int t)
{
	DMatrix *trajHMM = (DMatrix*)New(&gstack,3*sizeof(DMatrix));
	DMatrix meanVect,trajVarMat,rq,temp1;
	DMatrix trajMeanVect;
	/*DVector trajMeanVect;*/

	 /*% Transform into single-column DVectors*/

	meanVect = selectDVectors(meanMfcc, selMat);

	/*trajVarMat = inv(W' * invVarMat * W); % Pq*/

	temp1 = mulMatrix2Op(W,selectDVectors(calNumWDVector('d',1,varMfcc), selMat));

	trajVarMat = mulMatrixOp9(temp1,W);

	DMatrixInvert(trajVarMat);

	/*rq = W' * invVarMat * meanVect;*/


	rq = mulMatrixOp12(temp1,meanVect);

	/*
	trajMeanVect = trajVarMat * rq; % cq
	trajMeanVect = trajMeanVect'; % make it row DVector
	*/

	trajMeanVect = mulMatrixOp1(trajVarMat,rq);


	trajHMM[1] = trajMeanVect;
	trajHMM[2] = trajVarMat;

	/*FreeDMatrix(&gstack,trajMeanVect);*/
	return trajHMM;
}

/*
%%%
% data contain multiple row DVectors
% Select row DVectors based on selMat, then stack them vertically
%%%
*/

DMatrix selectDVectors(DVector data, DMatrix selMat)
{
	DMatrix newData;
	int i,j,col,row;
	col = NumDCols(selMat);
	row = DVectorSize(data);

	newData = CreateDMatrix(&gstack,col*row,1);

	for(i=0;i<col;i++)
	{
		for(j=1;j<=row;j++)
		newData[i*row+j][1] = data[j]*selMat[1][i+1];
	}

	return newData;
}

DVector moveC0Top(DVector data,int M, int D)
{
	DVector newData;
	newData = CreateDVector(&gstack,39);

	int i,startId,endId,k;
	for(i=1;i<=D;i++)
	{
		startId = (i-1)*M+1;
		endId	= i*M;
		/*printf("startId=%d\n endId = %d \n ans=%e",startId,endId,data[endId]);*/
		newData[startId]=data[endId];

		for(k=startId;k<endId;k++)
		{
			newData[k+1] = data[k];
		}

	}

	return newData;
}

/* newData(:, startId:endId) = [data(:, (startId+1):endId) data(:, startId)];*/
DVector moveC0Down(DVector data,int M, int D)
{
	DVector newData;
	newData = CreateDVector(&gstack,39);
	/*CopyDVector(data,newData);*/
	int i,startId,endId,k;
	for(i=1;i<=D;i++)
	{
		startId = (i-1)*M+1;
		endId	= i*M;
		/*printf("startId=%d\n endId = %d \n ans=%e",startId,endId,data[endId]);*/
		newData[endId]=data[startId];

		for(k=startId;k<endId;k++)
		{
			newData[k] = data[k+1];
		}

	}

	return newData;
}



/* --------------------- Initialisation ----------------------- */

/*
% Process 'option' and 'wOpt' and return/change configuration
% variables
*/
/*[option, wOpt, windowL, numPoints, isRetrain, noisyPrefix, numDctDVectors] = processOptions(option, wOpt, numPoints, noisyPrefix, numObsVectors)*/

void processOptions(char **option, char **wOpt, int *numPoints, int numObsVectors,int *windowL, int *numDctDVectors)
{
	char mfccDir[50];
	char *nwOpt =*wOpt;
	/*%% check for dpmc or retrain option %%*/
/*%% check to get windowL and wOpt %%	*/

		if(nwOpt[0]!='l' || nwOpt[3]!='w')
		{
			printf("#! Error:wOpt %s is not of the form l$windowL-wOpt\n",nwOpt);

		}else
		{
			*windowL=atoi(&nwOpt[1]); /* dynamic features window*/
			*wOpt = &nwOpt[3];
			/*strncpy(*wOpt,&nwOpt[3],strlen(nwOpt)-2);*/
		}
		/*%% numDctDVectors %%*/
		*numDctDVectors = getNumDctDVectors(numObsVectors, *wOpt, *windowL);

}

int getNumDctDVectors(int numObsVectors, char *wOpt, int windowL)
{
	int numDctDVectors;

	if(strcmp(wOpt,"wNormal")==0)
		numDctDVectors = numObsVectors+4*windowL;
	else
		numDctDVectors = numObsVectors; /*%% also considered numTrajDVectors*/

	return numDctDVectors;
}

