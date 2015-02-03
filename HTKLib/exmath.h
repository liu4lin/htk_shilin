DMatrix reshape(DMatrix A, int m, int n);
char* stringcat(char *des, char * from );
char** TextRead(char *filepath,DVector info);
DMatrix toRow(DMatrix m);
void initDMatrix(int start,int step, int end,DMatrix m);
DMatrix cosDMatrix(DMatrix m);
DMatrix sinDMatrix(DMatrix m);
DMatrix expDMatrix(DMatrix m);
DMatrix logDMatrix(DMatrix m);
DMatrix powDMatrix(double num,DMatrix m);
DMatrix mulDMatrixRC(DMatrix m1,DMatrix m2);
DMatrix mulDMatrix(DMatrix m1,DMatrix m2,int type);
DMatrix mulDMatrix1(DMatrix m1,DMatrix m2,int type);
DMatrix calDMatrix(DMatrix m1,DMatrix m2,int type);
DMatrix sumEDMatrix1(DMatrix m1,DMatrix m2);
DMatrix minusEDMatrix(DMatrix m1,DMatrix m2);
DMatrix mulEDMatrix(DMatrix m1,DMatrix m2);
DMatrix divEDMatrix(DMatrix m1,DMatrix m2);

DMatrix mulNumberwDMatrix(double num, DMatrix m);
DMatrix sumNumberwDMatrix(double num, DMatrix m);
double sumEDMatrix(DMatrix m);
DMatrix divNumberwDMatrix(double num, DMatrix m);
DMatrix calNumWDMatrix(char type,double num,DMatrix m);
DVector calNumWDVector(char type,double num,DVector v);
DMatrix eye(int nr, int nc);
DMatrix ones(int nr, int nc);
void ShowDMatrix1(char * title,DMatrix m);
DMatrix kron(DMatrix a, DMatrix b);
DMatrix cutColDMatrix(DMatrix m, int col);
DMatrix diag(DMatrix m);
DMatrix*  initializeDctMatrices(int numDctDVectors, int numChans, int numceps, int cepLifter);
DMatrix* initializeMatrices(int numObsDVectors, int numDctDVectors,int numChans,int numceps,int cepLifter,int M,int D,int windowL,char *option, char* wOpt);
DMatrix* getDynamicFeatureDMatrix(int M, int D, int windowL, int T, char* wOpt);
DMatrix getVSMat(int M,int D,int numDctDVectors,DMatrix W,DMatrix W0);
double strtodouble(char *s);
DMatrix mulMatrix(DMatrix m1,DMatrix m2);
DMatrix mulMatrixTest(DMatrix m1,DMatrix m2);
/* Mul Matrix m1 with m2, type=1: m1*m2 with nc1 = nc2*/
DMatrix mulMatrixOp(DMatrix m1,DMatrix m2);
/* Mul Matrix m1 with m2, type=1: m1*m2 with nc1 =nr1*/
DMatrix mulMatrixOp1(DMatrix m1,DMatrix m2);
/* Mul Matrix m1 with m2, type=1: m1*m2
DMatrix mulMatrixOp2(DMatrix m1,DMatrix m2);
/*  this function only for a tri-diagonal matrix m1 mul with m2 and the same size*/
DMatrix mulMatrix2Op_t(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp3(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp4(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp5(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp6(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp7(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp8(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp9(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp10(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp11(DMatrix m1,DMatrix m2);
DMatrix mulMatrixOp12(DMatrix m1,DMatrix m2);
DMatrix mulMatrix2(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3(DMatrix m1,DMatrix m2);
/* m1 = m2*/
DMatrix mulMatrix2Op1(DMatrix m1);
/* Mul Matrix m1 with m2:m1*m2' nc1 = 2*nr2*/
DMatrix mulMatrix3Op1(DMatrix m1,DMatrix m2);
/* Mul Matrix m1 with m2:m1*m2' nc1 = nr2*/
DMatrix mulMatrix3Op(DMatrix m1,DMatrix m2);
/* Mul DMatrix m1 with m2:m1'*m2 with nr1=nc2*/
DMatrix mulMatrix2Op(DMatrix m1,DMatrix m2);
/* Mul DMatrix m1 with m2:m1*m2' j=1*/
DMatrix mulMatrix3Op2(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3Op3(DMatrix m1,DMatrix m2);
/* only for m2 is factmatrix*/
DMatrix mulMatrix3Op4(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3Op5(DMatrix m1,DMatrix m2);
/* Mul Matrix test*/
DMatrix mulMatrix3Op6(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3Op7(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3Op8(DMatrix m1,DMatrix m2);
DMatrix mulMatrix3Op9(DMatrix m1,DMatrix m2);
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Test(DMatrix m1,DMatrix m2);
DMatrix mulMatrix2Test(DMatrix m1,DMatrix m2);


/* Inverting a double matrix */
void DMatrixInvert(DMatrix c);
void DLinSolve1(DMatrix a, int *perm, double *b);
Boolean DLUDecompose1(DMatrix a, int *perm, int *sign);
