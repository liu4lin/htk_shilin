#include "HShell.h"
#include "HMem.h"
DMatrix * pmc(int phoneId,char* phone,int id,DMatrix noiseMeanLinear,DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DVector cleanMeanMfcc,DVector cleanVarMfcc,double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix,DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,int numComponents,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix* postprocess(int phoneId,char* phone,int mixtureId,DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanMeanLog,DMatrix cleanVarLog,double alpha,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,DMatrix dctDMatrix,DMatrix cepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix * invTrajectoryHMM(DMatrix trajMeanVect,DMatrix trajVarMat,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int D,int numObsVectors,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix mfcc(DMatrix dataLog,DMatrix dctDMatrix,DMatrix cepDMatrix,int type,char* covarOpt,char* option,int test);
DMatrix computeMeanAvg(int numObsVectors,int numFeatures,DMatrix cHat);
DMatrix computeVarAvg(int numObsVectors,int numFeatures,DMatrix PHat,DMatrix factDMatrix);
DMatrix* linear2log(DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanVarLog,double alpha,char* covarOpt,char* option);
DMatrix * modelPMC(int phoneId,char* phone, DMatrix noiseMeanLinear, DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog, DVector* stateCleanMeanData, DVector* stateCleanVarData, int numMixture, double alpha, int M, int numObsVectors, int numDctDVectors, int D,DMatrix W, DMatrix VS, DMatrix Sq, int windowL, int numChans,int numceps, DMatrix dctDMatrix, DMatrix cepDMatrix, DMatrix invDctDMatrix, DMatrix invCepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix, int isDebug);
DMatrix* preprocess(DVector noiseMeanMfcc, DVector noiseVarMfcc, int M, DMatrix Sq, int numDctDVectors, int D,DMatrix W, int windowL, DMatrix invDctDMatrix, DMatrix invCepDMatrix, char *covarOpt, char* option,DMatrix selMat, int isDebug,int t);
DMatrix* independentPMC(int phoneId,char* phone,DMatrix noiseMeanLinear,DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DVector* cleanMeanData, DVector* cleanVarData, int numMixture, double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS, DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix, DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug);
DMatrix invMfcc(DMatrix dataMfcc, DMatrix invDctDMatrix,DMatrix invCepDMatrix, int type, char *covarOpt, char *option,int t);
DMatrix *trajectoryHMM(DVector meanMfcc, DVector varMfcc,DMatrix Sq, DMatrix W, int windowL, int M, char* option,DMatrix selMat, int isDebug, int t);
DMatrix selectDVectors(DVector data, DMatrix selMat);
DVector moveC0Top(DVector data,int M, int D);
DVector moveC0Down(DVector data,int M, int D);
void processOptions(char **option, char **wOpt, int *numPoints, int numObsVectors,int *windowL, int *numDctDVectors);
