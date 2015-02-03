#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "exmath.h"
#include "xtpmc.h"

DMatrix * modelPMC(int phoneId,char* phone, DMatrix noiseMeanLinear, DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog, DVector* stateCleanMeanData, DVector* stateCleanVarData, int numMixture, double alpha, int M, int numObsVectors, int numDctDVectors, int D,DMatrix W, DMatrix VS, DMatrix Sq, int windowL, int numChans,int numceps, DMatrix dctDMatrix, DMatrix cepDMatrix, DMatrix invDctDMatrix, DMatrix invCepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix, int isDebug)
{
	DMatrix *stateNoisyData;

	stateNoisyData = independentPMC(phoneId, phone,noiseMeanLinear, noiseVarLinear, noiseMeanLog, noiseVarLog, stateCleanMeanData, stateCleanVarData, numMixture, alpha, M, numObsVectors, numDctDVectors, D, W, VS, Sq, windowL, numChans, numceps, dctDMatrix, cepDMatrix, invDctDMatrix, invCepDMatrix, covarOpt, option, selMat,factDMatrix, isDebug);

	/*ShowDMatrix1("",stateMeanData);*/
	return stateNoisyData;
}

DMatrix* independentPMC(int phoneId,char* phone,DMatrix noiseMeanLinear,DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DVector *cleanMeanData, DVector *cleanVarData, int numMixture, double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS, DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix, DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug)
{

	DMatrix *stateNoisyData = (DMatrix*)New(&gstack, 3*sizeof(DMatrix));
	DMatrix noisyMeanData, noisyVarData, trajMeanVectData;
	DMatrix * noisyData;
	int i,id,k,n = DVectorSize(cleanMeanData[1]);
	DVector cleanMeanMfcc,cleanVarMfcc;
	DVector	noisyMeanMfcc,noisyVarMfcc;
	noisyMeanData = CreateDMatrix(&gstack,numMixture,n);
	noisyVarData = CreateDMatrix(&gstack,numMixture,n);

	for(id=1;id<=numMixture;id++)
		{
			/*cleanMeanMfcc = CreateDVector(&gstack,n);*/
			/*cleanVarMfcc = CreateDVector(&gstack,n);*/
			cleanMeanMfcc = cleanMeanData[id];
			cleanVarMfcc = cleanVarData[id];
			/*ShowDVector("",cleanMeanMfcc,39);*/
			double cleanSum=0,cleanVar=0;
			for(k=1;k<=n;k++)
			{

				cleanSum+=fabs(cleanMeanMfcc[k]);

				cleanVar += fabs(cleanVarMfcc[k]);
			}

			if( cleanSum == 0 && cleanVar == 0)
			{

					noisyMeanMfcc = CreateDVector(&gstack,39);
				 	noisyVarMfcc = CreateDVector(&gstack,39);
			}else
			{
				/*%% Sanity check %%*/
				for(k=1;k<=n;k++)
					if(cleanVarMfcc[k] == 0)
						printf("#! Variance inputs are zeroes at:%d",k);
				/*% PMC*/
			noisyData = pmc(phoneId, phone, id,noiseMeanLinear, noiseVarLinear, noiseMeanLog, noiseVarLog, cleanMeanMfcc, cleanVarMfcc, alpha, M, numObsVectors, numDctDVectors, D, W, VS, Sq, windowL, numChans, numceps, dctDMatrix, cepDMatrix, invDctDMatrix, invCepDMatrix, covarOpt, option, 1, selMat,factDMatrix, isDebug); /*% numComponents = 1*/

		/*	% Assign results*/
			noisyMeanData[id] = noisyData[1][1];
			noisyVarData[id] = noisyData[2][1];

			}

		}


	stateNoisyData[1] = noisyMeanData;
	stateNoisyData[2] = noisyVarData;

	return stateNoisyData;
}
/* return noisyMeanMfcc, noisyVarMfcc, trajMeanVect, trajVarMat*/
DMatrix * pmc(int phoneId,char* phone,int mixtureId,DMatrix noiseMeanLinear,DMatrix noiseVarLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DVector cleanMeanMfcc,DVector cleanVarMfcc,double alpha,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,int numceps,DMatrix dctDMatrix,DMatrix cepDMatrix,DMatrix invDctDMatrix,DMatrix invCepDMatrix,char* covarOpt,char* option,int numComponents,DMatrix selMat,DMatrix factDMatrix,int isDebug)
{
 DMatrix * noisyData;
 DMatrix  trajMeanVect, trajVarMat;
 DMatrix *preProDMatrix,cleanMeanLinear, cleanVarLinear, cleanMeanLog, cleanVarLog;

 preProDMatrix = preprocess(cleanMeanMfcc, cleanVarMfcc, M, Sq, numDctDVectors, D, W, windowL, invDctDMatrix, invCepDMatrix, covarOpt, option, selMat, isDebug,1);
 cleanMeanLinear = preProDMatrix[1];
 cleanVarLinear = preProDMatrix[2];
 cleanMeanLog = preProDMatrix[3];
 cleanVarLog = preProDMatrix[4];
 cleanMeanMfcc = preProDMatrix[5][1];
 cleanVarMfcc = preProDMatrix[6][1];
 trajMeanVect = preProDMatrix[7];
 trajVarMat = preProDMatrix[8];

 /*% Combine model*/
 DMatrix noisyMeanLinear,noisyVarLinear;
 noisyMeanLinear = mulNumberwDMatrix(alpha,noiseMeanLinear);
 noisyMeanLinear = sumEDMatrix1(cleanMeanLinear,noisyMeanLinear);
 noisyVarLinear = mulNumberwDMatrix(alpha*alpha,noiseVarLinear);
 noisyVarLinear = sumEDMatrix1(cleanVarLinear,noisyVarLinear);

/*
 *
	%%% Post process: linear -> ceps domain %%%
 */

 noisyData= postprocess(phoneId, phone,mixtureId, noisyMeanLinear, noisyVarLinear, noiseMeanLinear, noiseMeanLog, noiseVarLog, cleanMeanLinear, cleanMeanLog, cleanVarLog, alpha, cleanMeanMfcc, cleanVarMfcc, M, numObsVectors, numDctDVectors, D, W, VS, Sq, windowL, numChans, dctDMatrix, cepDMatrix, covarOpt, option, selMat,factDMatrix, isDebug);

 /*noisyData[3] = trajMeanVect;*/
 /*noisyData[4] =trajVarMat;*/


return noisyData;

}

/*%%% Post process: linear -> ceps domain %%%
 *  Return DMatrix with meanMfcc and varMfcc
 */

DMatrix* postprocess(int phoneId,char* phone,int mixtureId,DMatrix meanLinear,DMatrix varLinear,DMatrix noiseMeanLinear,DMatrix noiseMeanLog,DMatrix noiseVarLog,DMatrix cleanMeanLinear,DMatrix cleanMeanLog,DMatrix cleanVarLog,double alpha,DVector cleanMeanMfcc,DVector cleanVarMfcc,int M,int numObsVectors,int numDctDVectors,int D,DMatrix W,DMatrix VS,DMatrix Sq,int windowL,int numChans,DMatrix dctDMatrix,DMatrix cepDMatrix,char* covarOpt,char* option,DMatrix selMat,DMatrix factDMatrix,int isDebug)
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
	 meanMfcc = mfcc(meanLog, dctDMatrix, cepDMatrix,1 , covarOpt, option,0);
	 int test=0;
	 if(phoneId==7 && mixtureId==5)
		 test=1;
	 varMfcc = mfcc(varLog, dctDMatrix, cepDMatrix, 2, covarOpt, option,test);

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

	/*free(logData);
	free(dataMfcc);*/

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
DMatrix mfcc(DMatrix dataLog,DMatrix dctDMatrix,DMatrix cepDMatrix,int type,char* covarOpt,char* option,int test)
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

 DMatrix *preProDMatrix  = (DMatrix*)New(&gstack,9*sizeof(DMatrix)) ;
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
 for(i=1; i<=n;i++)
 {
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

 if (strcmp(option,"trajAvg")==0)
 	{
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



