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
/*         File: HMath.h:   Math Support                       */
/* ----------------------------------------------------------- */

/* !HVER!HMath:   3.4.1 [CUED 12/03/09] */

#ifndef _HMATH_H_
#define _HMATH_H_


#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PI
#undef PI                /* PI is defined in Linux */
#endif
#define PI   3.14159265358979
#define TPI  6.28318530717959     /* PI*2 */
#define LZERO  (-1.0E10)   /* ~log(0) */
#define LSMALL (-0.5E10)   /* log values < LSMALL are set to LZERO */
#define LSMALLP (-23.0259)	   /* small log posterior value, used for sparse dnn decoding*/
#define SMALLP (1.0E-10)	/* small posterior, floor value*/
#define MINEARG (-708.3)   /* lowest exp() arg  = log(MINLARG) */
#define MINLARG 2.45E-308  /* lowest log() arg  = exp(MINEARG) */
#define BIG_FLOAT 1E20      /* Limit for float arguments */

/* NOTE: On some machines it may be necessary to reduce the
         values of MINEARG and MINLARG
*/

typedef float  LogFloat;   /* types just to signal log values */
typedef double LogDouble;

typedef enum {  /* Various forms of covariance matrix */
   DIAGC,         /* diagonal covariance */
   INVDIAGC,      /* inverse diagonal covariance */
   FULLC,         /* inverse full rank covariance */
   XFORMC,        /* arbitrary rectangular transform */
   LLTC,          /* L' part of Choleski decomposition */
   NULLC,         /* none - implies Euclidean in distance metrics */
   NUMCKIND       /* DON'T TOUCH -- always leave as final element */
} CovKind;

typedef enum {
	Gaussi,
	Euclid,
	Userdt
}DistKind;

typedef union {
   SVector var;         /* if DIAGC or INVDIAGC */
   STriMat inv;         /* if FULLC or LLTC */
   SMatrix xform;       /* if XFORMC */
} Covariance;

typedef struct _Point{
	Vector value;
	Ptr owner;
	Ptr hook;
	double dotcache;
	double occcache;
	struct _Point* next;
}*Point;
typedef struct _PntGrp{
	Point plist;/*a list of regweights belonging to this cluster*/
	Vector mean;/*cluster center*/
	Vector var;
	TriMat inv;/*Keep the inversion of full covariance matrix*/
	Vector acc1;/*1st order statistics for mean, mu=acc1/n*/
	Vector acc2;/*2nd order statistics for var, var=(acc2/n)-mu*mu*/
	TriMat acc2x;/*for full covariance*/
	float occ;
	int hardOcc;
	LogFloat gConst;/*-0.5*(vsize*log(2pi)+log(det))*/
	float weight;
	Boolean fc;
	float *floorV;
	DistKind distK;
}*Group;


/* ------------------------------------------------------------------- */

void InitMath(void);
/*
   Initialise the module
*/

/* ------------------ Vector Oriented Routines ----------------------- */

void ZeroShortVec(ShortVec v);
void ZeroIntVec(IntVec v);
void ZeroVector(Vector v);
void ZeroDVector(DVector v);
/*
   Zero the elements of v
*/

void CopyShortVec(ShortVec v1, ShortVec v2);
void CopyIntVec(IntVec v1, IntVec v2);
void CopyVector(Vector v1, Vector v2);
void CopyNVector(Vector v1, Vector v2);
void CopyiVector(Vector v1, Vector v2);
void CopyDVector(DVector v1, DVector v2);
/*
   Copy v1 into v2; sizes must be the same
*/

Boolean ReadShortVec(Source *src, ShortVec v, Boolean binary);
Boolean ReadIntVec(Source *src, IntVec v, Boolean binary);
Boolean ReadVector(Source *src, Vector v, Boolean binary);
/*
   Read vector v from source in ascii or binary
*/

void WriteShortVec(FILE *f, ShortVec v, Boolean binary);
void WriteIntVec(FILE *f, IntVec v, Boolean binary);
void WriteVector(FILE *f, Vector v, Boolean binary);
/*
   Write vector v to stream f in ascii or binary
*/

void ShowShortVec(char * title, ShortVec v,int maxTerms);
void ShowIntVec(char * title, IntVec v,int maxTerms);
void ShowVector(char * title,Vector v,int maxTerms);
void ShowDVector(char * title,DVector v,int maxTerms);
/*
   Print the title followed by upto maxTerms elements of v
*/

/* Quadratic prod of a full square matrix C and an arbitry full matrix transform A */
void LinTranQuaProd(Matrix Prod, Matrix A, Matrix C);

/* ------------------ Matrix Oriented Routines ----------------------- */

void ZeroMatrix(Matrix m);
void ZeroDMatrix(DMatrix m);
void ZeroTriMat(TriMat m);
void IdentityTriMat(TriMat m);
/*
   Zero the elements of m
*/

void CopyMatrix (Matrix m1,  Matrix m2);
void CopyNMatrix(Matrix m1, Matrix m2);
void CopyDMatrix(DMatrix m1, DMatrix m2);
void CopyTriMat (TriMat m1,  TriMat m2);
void CopyNTriMat(TriMat m1, TriMat m2);
/*
   Copy matrix m1 to m2 which must have identical dimensions
*/

void Mat2DMat(Matrix m1,  DMatrix m2);
void DMat2Mat(DMatrix m1, Matrix m2);
void Mat2Tri (Matrix m1,  TriMat m2);
void DMat2Tri (DMatrix m1,  TriMat m2);
void Tri2Mat (TriMat m1,  Matrix m2);
void Tri2DMat (TriMat m1,  DMatrix m2);
void DTri2DMat (DMatrix m1, DMatrix m2);
void Tri2Vec (TriMat m, Vector v);
void Mat2Vec (Matrix m, Vector v);
void Vec2Mat(Vector v, Matrix m);
void Vec2Tri(Vector v, TriMat m);
/*void Tri2DMat (TriMat m1, DMatrix m2);*/
/*
   Convert matrix format from m1 to m2 which must have identical 
   dimensions
*/

Boolean ReadMatrix(Source *src, Matrix m, Boolean binary);
Boolean ReadTriMat(Source *src, TriMat m, Boolean binary);
/*
   Read matrix from source into m using ascii or binary.
   TriMat version expects m to be in upper triangular form
   but converts to lower triangular form internally.
*/
   
void WriteMatrix(FILE *f, Matrix m, Boolean binary);
void WriteTriMat(FILE *f, TriMat m, Boolean binary);
/*
    Write matrix to stream in ascii or binary.  TriMat version 
    writes m in upper triangular form even though it is stored
    in lower triangular form!
*/

void ShowMatrix (char * title,Matrix m, int maxCols,int maxRows);
void ShowDMatrix(char * title,DMatrix m,int maxCols,int maxRows);
void ShowTriMat (char * title,TriMat m, int maxCols,int maxRows);
/*
   Print the title followed by upto maxCols elements of upto
   maxRows rows of m.
*/

/* ------------------- Linear Algebra Routines ----------------------- */

LogFloat CovInvert(MemHeap* x, TriMat c, Matrix invc);
/*
   Computes inverse of c in invc and returns the log of Det(c),
   c must be positive definite.
*/
LogDouble DCovInvert(TriMat c, DMatrix invc);
LogFloat CovIntp(TriMat c, Matrix invc, Boolean inv);
LogDouble DCovIntp(DMatrix c, DMatrix invc, Boolean inv);

LogFloat CovDet(MemHeap* x, TriMat c);
LogDouble DCovDet(DMatrix c);
/*
   Returns log of Det(c), c must be positive definite.
*/

/* EXPORT->MatDet: determinant of a matrix */
float MatDet(Matrix c);

/* EXPORT->DMatDet: determinant of a double matrix */
double DMatDet(DMatrix c);

/* EXPORT-> MatInvert: puts inverse of c in invc, returns Det(c) */
float MatInvert(Matrix c, Matrix invc);
double DMatInvert(MemHeap* x, DMatrix c, DMatrix invc);
 
/* DMatCofact: generates the cofactors of row r of doublematrix c */
double DMatCofact(DMatrix c, int r, DVector cofact);

/* MatCofact: generates the cofactors of row r of doublematrix c */
double MatCofact(Matrix c, int r, Vector cofact);

/* ------------- Singular Value Decomposition Routines --------------- */

void SVD(DMatrix A, DMatrix U,  DMatrix V, DVector d);
/* 
   Singular Value Decomposition (based on MESCHACH)
   A is m x n ,  U is m x n,  W is diag N x 1, V is n x n
*/

void InvSVD(DMatrix A, DMatrix U, DVector W, DMatrix V, DMatrix Result);
/* 
   Inverted Singular Value Decomposition (calls SVD)
   A is m x n ,  U is m x n,  W is diag N x 1, V is n x n, Result is m x n 
*/

/* ------------------- Log Arithmetic Routines ----------------------- */

LogDouble LAdd(LogDouble x, LogDouble y);
/*
   Return x+y where x and y are stored as logs, 
   sum < LSMALL is floored to LZERO 
*/

LogDouble LSub(LogDouble x, LogDouble y);
/*
   Return x-y where x and y are stored as logs, 
   diff < LSMALL is floored to LZERO 
*/

double L2F(LogDouble x);
/*
   Convert log(x) to real, result is floored to 0.0 if x < LSMALL 
*/

/* ------------------- Random Number Routines ------------------------ */

void RandInit(int seed);
/* 
   Initialise random number generators, if seed is -ve, then system 
   clock is used.  RandInit(-1) is called by InitMath.
*/

float RandomValue(void);
/*
   Return a random number in range 0.0->1.0 with uniform distribution
*/

float GaussDeviate(float mu, float sigma);
/*
   Return a random number with a N(mu,sigma) distribution
*/
Boolean Choleski(TriMat A, DMatrix L);
Boolean DCholeski(DMatrix A, DMatrix L);
Boolean IsPositiveDMat(DMatrix A);
Boolean IsPositiveTriM(TriMat A);
Boolean IsSingularMat(Matrix c);
void EigenDecomposition(TriMat A, Matrix U, Vector d);
LogFloat ModDiagHessian(MemHeap* x, TriMat Hessian_from, Matrix Hessian_to, Boolean inv, float *tau_out, Boolean trace);
LogDouble DModDiagHessian(DMatrix Hessian_from, DMatrix Hessian_to, Boolean inv, double *tau_out, Boolean trace);
LogFloat ModEigHessian(TriMat Hessian_from, Matrix Hessian_to, Boolean inv, int *floor_out, Boolean trace);
LogDouble DModEigHessian(DMatrix Hessian_from, DMatrix Hessian_to, Boolean inv, int *floor_out, Boolean trace);
void Identity(Matrix A);
void DIdentity(DMatrix A);
void DiagMatrix(Matrix A);
void DiagTriMat(TriMat A);
float P_norm(Vector v, int p);
float XP_norm(Matrix x, int p);
Boolean DotProduct(Vector v1,Vector v2, float *res);
Boolean MulDotProduct(Vector v1,Vector v2, short *mul, float *res);
Boolean MulHardDotProduct(Vector v1,Vector v2, short *mul, float *res);
Boolean DDotProduct(DVector v1,DVector v2, double *res);
Boolean QuaProduct(Vector v1,Vector v2, Matrix res);
Boolean Transpose(Matrix m1,Matrix m1t);
Boolean DTranspose(DMatrix m1,DMatrix m1t);
float Trace(Matrix m1);
double DTrace(DMatrix m1);
Boolean MatColVector(Matrix m, Vector v, int j);
Boolean MatMult(Matrix m1, Matrix m2, Matrix res);
Boolean SqMatMult(Matrix m1, Matrix m2, Matrix res);
Boolean DMatDiag(DMatrix A, DVector d, DMatrix Ad);
Boolean DMatDiagDMat(DMatrix A, DVector d, DMatrix B, DMatrix AdB);
Boolean DMatMult(DMatrix m1, DMatrix m2, DMatrix res);
Boolean DSqMatMult(DMatrix m1, DMatrix m2, DMatrix res);
Boolean VecAdd(Vector v1, float scale1, Vector v2, float scale2, Vector res);
Boolean DVecAdd(DVector v1, double scale1, DVector v2, double scale2, DVector res);
Boolean MatAdd(Matrix m1, float scale1, Matrix m2, float scale2, Matrix res);
Boolean DMatAdd(DMatrix m1, double scale1, DMatrix m2, double scale2, DMatrix res);
Boolean MatInterp(Matrix m,Vector v,Vector res);
Boolean rMatInterp(Vector v,Matrix m,Vector res);
Boolean VecSquare(Vector v1, Matrix m, Vector v2, float *res);
Boolean MatScale(Matrix m,float s,Matrix res);
Boolean VecScale(Vector v,float s,Vector res);
float EuclMatDist(Matrix m1, Matrix m2);
void printMatrix(Matrix mat, char* name);
void printDMatrix(DMatrix mat, char* name);
void printTriMat(TriMat tri, char* name);
void printVector(Vector vec, char* name);
void printDVector(DVector vec, char* name);

void NNMF(Matrix V, Matrix W, Matrix H, int k, int maxIter, Boolean speakout);
void NNMFpos(Matrix V, Matrix H, int k, int maxIter, Boolean speakout);

Boolean LUDecompose(Matrix a, int *perm, int *sign);
void LinSolve(Matrix a, int *perm, float *b);
/*---------------------Signal processing related--------------------*/
/*
 * linear spectral domain<-----log&exp------->log spectral domain
 * log spectral domain<--------DCT&iDCT------->cepstral domain
 * cepstral coefficients<-------cepLifter&icepLifter-------->liftered cepstral coefficients
 */
/*DMatrix DCT(int numChans, int numCeps);
DMatrix iDCT(int numChans, int numCeps);*/
DMatrix liftDCT(MemHeap* x, int numChans, int numCeps, int numLifter);
DMatrix iLiftDCT(MemHeap* x, int numChans, int numCeps, int numLifter);
DMatrix iLiftDCT2(int numChans, int numCeps, int numLifter);
DMatrix iLiftDCTx(int numChans, int numCeps, int numLifter);
/*void cepLift(int numCeps, int numLifter, DVector mu, DMatrix cov);
void icepLift(int numCeps, int numLifter, DVector mu, DMatrix cov);*/
/*DMatrix trajDCT(int trajLen, int numChans, int numCeps);
DMatrix itrajDCT(int trajLen, int numChans, int numCeps);*/
DMatrix lLiftDCT(MemHeap* x, int trajLen, int numChans, int numCeps, int numLifter);
DMatrix ilLiftDCT(MemHeap* x, int trajLen, int numChans, int numCeps, int numLifter);
void VTS(MemHeap* x, DMatrix dct, DMatrix idct, DVector ex_cep_mu, DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DVector g0);
void VtsUpVar(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DMatrix GCx, DMatrix GCxGt, DMatrix FCn, DMatrix FCnFt, DMatrix x_cep_cv, DMatrix n_cep_cv, DMatrix y_cep_cv);
void dVtsUpVar(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DVector x_cep_cv, DVector n_cep_cv, DVector y_cep_cv);
void cep2linear(DVector c_mu, DMatrix iDct, DMatrix c_cov, DVector l_mu, DMatrix l_cov);
void linear2cep(DVector l_mu, DMatrix dct, DMatrix l_cov, DVector c_mu, DMatrix c_cov);
void dvec2diag(DVector v, DMatrix m);
void DupDVector(DVector v, int numElem, DVector longv);
void DupDMatrix(DMatrix m, int numRow, int numCol, DMatrix longm);
void duplicatemdiag(DMatrix m, int numRow, int numCol, DMatrix longm);
void dvec2idiag(DVector v, DMatrix m);
void Vec2iDVec(Vector v, DVector dv);
void Vec2DVec(Vector v, DVector dv);
void DVec2Vec(DVector dv, Vector v);
void truncDVec2Vec(DVector dv, Vector v);
void truncDVec2DVec(DVector dv1, DVector dv2);
void truncDMat2DMat(DMatrix dm1, DMatrix dm2);
void truncDMat2Mat(DMatrix dm1, Matrix dm2);
void MatDiag2iVec(Matrix m, Vector v);
void MatDiag2Vec(Matrix m, Vector v);
void DMatDiag2Vec(DMatrix m, Vector v);
void kronecker(DMatrix a, DMatrix b, DMatrix res);
/*-----------------------GMM estimation related-----------------*/
Point CreatePoint(MemHeap* x, Vector value, Ptr owner, Ptr hook);
void FreePoint(MemHeap* x, Point pnt);
Group CreateGroup(MemHeap* x, int vsize, Boolean fc);
float Euclidean(Vector v1, Vector v2);
void AddIntoGroup(Group grp, Point pnt, float post, Boolean insert, Boolean callOcc);
float SearchAndInsert(MemHeap* x, Group* grps, Point pnt, int numClust, DistKind distK, Boolean insert, Boolean callOcc);
void UpdateGroup(MemHeap* x, Group grp, int numPnts, Boolean upVar);
void UpdateGroups(MemHeap *x, Group* grps, int numClust, int numPnts, Boolean upVar);
void ResetGroups(Group* grps, int numClust);
Group* InitKmeansMU(Point* pnts, int numClust, int vSize, int numPnts, Boolean fc, Boolean tiedCov);

void Interior_point_method(Vector c, Matrix A, Vector b, Vector x, char *ident);
#ifdef __cplusplus
}
#endif

#endif  /* _HMATH_H_ */

/* ------------------------- End of HMath.h -------------------------- */
