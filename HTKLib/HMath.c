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
/*         File: HMath.c   Math Support Module                 */
/* ----------------------------------------------------------- */

char *hmath_version = "!HVER!HMath:   3.4.1 [CUED 12/03/09]";
char *hmath_vc_id = "$Id: HMath.c,v 1.1.1.1 2006/10/11 09:54:58 jal58 Exp $";

/*
   This library provides math support in the following three areas
   
      a) Vector/Matrix Operators
      b) Log Arithmetic
      c) Random Deviates
      
   It relies on the basic vector and matrix types defined by HMem.
   The separation of functionality betwee HMem and HMath is such
   that no routine in this module requires any knowledge of the
   hidden fields in these types.  Thus, if a change of representation
   is required, it should only be necessary to change HMem.
*/

#include "HShell.h"        /* HTK Libraries */
#include "HMem.h"
#include "HMath.h"

/* ----------------------------- Trace Flags ------------------------- */

static int trace = 0;

/* -------------------- Configuration Parameters --------------------- */

static ConfParam *cParm[MAXGLOBS];       /* config parameters */
static int numParm = 0;
static float min_eigen=1E-3;

/* ------------------ Vector Oriented Routines ----------------------- */

/*
   ShortVecs and IntVecs are pointers to arrays of short/int.  Vectors
   are pointers to arrays of float (ie float*).  All indexing is
   v[1..n].  The size is stored in v[0].
*/

/* EXPORT->ZeroShortVec: Zero the elements of v */
void ZeroShortVec(ShortVec v)
{
   int i,n;
   
   n=ShortVecSize(v);
   for (i=1;i<=n;i++) v[i]=0;
}

/* EXPORT->ZeroIntVec: Zero the elements of v */
void ZeroIntVec(IntVec v)
{
   int i,n;
   
   n=IntVecSize(v);
   for (i=1;i<=n;i++) v[i]=0;
}

/* EXPORT->ZeroVector: Zero the elements of v */
void ZeroVector(Vector v)
{
   int i,n;
   
   n=VectorSize(v);
   for (i=1;i<=n;i++) v[i]=0.0;
}

/* EXPORT->ZeroDVector: Zero the elements of v */
void ZeroDVector(DVector v)
{
   int i,n;
   
   n=DVectorSize(v);
   for (i=1;i<=n;i++) v[i]=0.0;
}

/* EXPORT->CopyShortVec: copy v1 into v2 */
void CopyShortVec(ShortVec v1, ShortVec v2)
{
   int i,size; 
   
   size = ShortVecSize(v1);
   if (size != ShortVecSize(v2))
      HError(5270,"CopyShortVec: sizes differ %d vs %d",
             size,ShortVecSize(v2));
   for (i=1; i<=size; i++) 
      v2[i] = v1[i];
}

/* EXPORT->CopyIntVec: copy v1 into v2 */
void CopyIntVec(IntVec v1, IntVec v2)
{
   int i,size; 
   
   size = IntVecSize(v1);
   if (size != IntVecSize(v2))
      HError(5270,"CopyIntVec: sizes differ %d vs %d",
             size,IntVecSize(v2));
   for (i=1; i<=size; i++) 
      v2[i] = v1[i];
}

/* EXPORT->CopyVector: copy v1 into v2 */
void CopyVector(Vector v1, Vector v2)
{
   int i,size; 
   
   size = VectorSize(v1);
   if (size != VectorSize(v2))
      HError(5270,"CopyVector: sizes differ %d vs %d",
             size,VectorSize(v2));
   for (i=1; i<=size; i++) 
      v2[i] = v1[i];
}

/* EXPORT->CopyNVector: copy negative v1 into v2 */
void CopyNVector(Vector v1, Vector v2)
{
   int i,size;

   size = VectorSize(v1);
   if (size != VectorSize(v2))
      HError(5270,"CopyVector: sizes differ %d vs %d",
             size,VectorSize(v2));
   for (i=1; i<=size; i++)
      v2[i] = -v1[i];
}

/* EXPORT->CopyiVector: copy the inverse of v1 into v2 */
void CopyiVector(Vector v1, Vector v2)
{
   int i,size;

   size = VectorSize(v1);
   if (size != VectorSize(v2))
      HError(5270,"CopyVector: sizes differ %d vs %d",
             size,VectorSize(v2));
   for (i=1; i<=size; i++)
      v2[i] = 1/v1[i];
}

/* EXPORT->CopyDVector: copy v1 into v2 */
void CopyDVector(DVector v1, DVector v2)
{
   int i,size; 
   
   size = DVectorSize(v1);
   if (size != DVectorSize(v2))
      HError(5270,"CopyDVector: sizes differ %d vs %d",
             size,DVectorSize(v2));
   for (i=1; i<=size; i++) 
      v2[i] = v1[i];
}

/* EXPORT->ReadShortVec: read vector from src in ascii or binary */ 
Boolean ReadShortVec(Source *src, ShortVec v, Boolean binary)
{
   return ReadShort(src,v+1,ShortVecSize(v),binary);
}

/* EXPORT->ReadIntVec: read vector from src in ascii or binary */ 
Boolean ReadIntVec(Source *src, IntVec v, Boolean binary)
{
   return ReadInt(src,v+1,IntVecSize(v),binary);
}

/* EXPORT->ReadVector: read vector from src in ascii or binary */ 
Boolean ReadVector(Source *src, Vector v, Boolean binary)
{
   return ReadFloat(src,v+1,VectorSize(v),binary);
}

/* EXPORT->WriteShortVec: write vector v to stream f */
void WriteShortVec(FILE *f, ShortVec v, Boolean binary)
{
   WriteShort(f,v+1,ShortVecSize(v),binary);
   if (!binary) fputc('\n',f);
}

/* EXPORT->WriteIntVec: write vector v to stream f */
void WriteIntVec(FILE *f, IntVec v, Boolean binary)
{
   WriteInt(f,v+1,IntVecSize(v),binary);
   if (!binary) fputc('\n',f);
}

/* EXPORT->WriteVector: write vector v to stream f */
void WriteVector(FILE *f, Vector v, Boolean binary)
{
   WriteFloat(f,v+1,VectorSize(v),binary);
   if (!binary) fputc('\n',f);
}

/* Export->ShowShortVec: show the short vector v preceded by title */
void ShowShortVec(char * title, ShortVec v,int maxTerms)
{
   int i, size, maxi;
   
   size = maxi = ShortVecSize(v);
   if (maxi>maxTerms) maxi=maxTerms;
   printf("%s\n   ",title);   
   for (i=1;i<=maxi;i++)  printf("%3d ",v[i]);
   if (maxi<size) printf("...");
   printf("\n");
}

/* Export->ShowIntVec: show the int vector v preceded by title */
void ShowIntVec(char * title, IntVec v,int maxTerms)
{
   int i, size, maxi;
   
   size = maxi = IntVecSize(v);
   if (maxi>maxTerms) maxi=maxTerms;
   printf("%s\n   ",title);   
   for (i=1;i<=maxi;i++)  printf("%5d ",v[i]);
   if (maxi<size) printf("...");
   printf("\n");
}

/* Export->ShowVector: show the vector v preceded by title */
void ShowVector(char * title,Vector v,int maxTerms)
{
   int i, size, maxi;
   
   size = maxi = VectorSize(v);
   if (maxi>maxTerms) maxi=maxTerms;
   printf("%s\n   ",title);   
   for (i=1;i<=maxi;i++)  printf("%8.2f ",v[i]);
   if (maxi<size) printf("...");
   printf("\n");
}

/* Export->ShowDVector: show the vector v preceded by title */
void ShowDVector(char * title, DVector v,int maxTerms)
{
   int i, size, maxi;
   
   size = maxi = DVectorSize(v);
   if (maxi>maxTerms) maxi=maxTerms;
   printf("%s\n   ",title);   
   for (i=1;i<=maxi;i++)  printf("%10.4f ",v[i]);
   if (maxi<size) printf("...");
   printf("\n");
}

/* ------------------ Matrix Oriented Routines ----------------------- */

/*
  Matrices are pointers to an array of vectors (ie float**).  Matrices
  are indexed m[1..r][1..c].  The rows of matrix are genuine vectors
  and can be treated as such except that matrices must be disposed of
  in their entirety.  If the rows of a matrix are substituted by
  other vectors then it should not be disposed.  TriMats are indexed
  m[1..r][1..i] where i is the row number ie only the lower triangle
  is stored.
*/

/* EXPORT->ZeroMatrix: Zero the elements of m */
void ZeroMatrix(Matrix m)
{
   int i,j,nr,nc;
   
   nr=NumRows(m); nc=NumCols(m);
   for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++) m[i][j]=0.0;
}

/* EXPORT->ZeroDMatrix: Zero the elements of m */
void ZeroDMatrix(DMatrix m)
{
   int i,j,nr,nc;
   
   nr=NumDRows(m); nc=DVectorSize(m[1]);
   for (i=1;i<=nr;i++)
      for (j=1;j<=nc;j++) m[i][j]=0.0;
}

/* EXPORT->ZeroTriMat: Zero the elements of m */
void ZeroTriMat(TriMat m)
{
   int i,j,size;
   
   size = TriMatSize(m);
   for (i=1;i<=size;i++)
      for (j=1;j<=i;j++) m[i][j]=0.0;
}

void IdentityTriMat(TriMat m)
{
   int i,j,size;

   size = TriMatSize(m);
   for (i=1;i<=size;i++)
      for (j=1;j<=i;j++) m[i][j]=i==j?1:0;
}

/* EXPORT->CopyMatrix: copy matrix m1 to m2 */
void CopyMatrix(Matrix m1, Matrix m2)
{
   int i,nrows;
   
   nrows = NumRows(m1);
   if (nrows != NumRows(m2))
      HError(5270,"CopyMatrix: row sizes differ %d vs %d",
             nrows,NumRows(m2));
   for (i=1; i<=nrows; i++)
      CopyVector(m1[i],m2[i]);
}

/* EXPORT->CopyNMatrix: copy negative matrix m1 to m2 */
void CopyNMatrix(Matrix m1, Matrix m2)
{
   int i,nrows;

   nrows = NumRows(m1);
   if (nrows != NumRows(m2))
      HError(5270,"CopyMatrix: row sizes differ %d vs %d",
             nrows,NumRows(m2));
   for (i=1; i<=nrows; i++)
      CopyNVector(m1[i],m2[i]);
}

/* EXPORT->CopyDMatrix: copy matrix m1 to m2 */
void CopyDMatrix(DMatrix m1, DMatrix m2)
{
   int i,nrows;
   
   nrows = NumDRows(m1);
   if (nrows != NumDRows(m2))
      HError(5270,"CopyDMatrix: row sizes differ %d vs %d",
             nrows,NumDRows(m2));
   for (i=1; i<=nrows; i++)
      CopyDVector(m1[i],m2[i]);
}

/* EXPORT->CopyTriMat: copy triangular matrix m1 to m2 */
void CopyTriMat(TriMat m1, TriMat m2)
{
   int i,size;
   
   size = TriMatSize(m1);
   if (size != TriMatSize(m2))
      HError(5270,"CopyTriMat: sizes differ %d vs %d",
             size,TriMatSize(m2));
   for (i=1; i<=size; i++)
      CopyVector(m1[i],m2[i]);
}

/* EXPORT->CopyNTriMat: copy negative triangular matrix m1 to m2 */
void CopyNTriMat(TriMat m1, TriMat m2)
{
   int i,size;

   size = TriMatSize(m1);
   if (size != TriMatSize(m2))
      HError(5270,"CopyTriMat: sizes differ %d vs %d",
             size,TriMatSize(m2));
   for (i=1; i<=size; i++)
      CopyNVector(m1[i],m2[i]);
}

/* EXPORT->Mat2DMat: convert matrix m1 to double matrix m2 */
void Mat2DMat(Matrix m1,  DMatrix m2)
{
   int i,j,nrows,ncols;

   nrows = NumRows(m1);
   if (nrows != NumDRows(m2))
      HError(5270,"Mat2DMat: row sizes differ %d vs %d",
             nrows,NumDRows(m2));
   ncols = NumCols(m1);
   if (ncols != NumDCols(m2))
      HError(5270,"Mat2DMat: col sizes differ %d vs %d",
             ncols,NumDCols(m2));
   for (i=1; i<=nrows; i++)
      for (j=1; j<=ncols; j++) 
         m2[i][j] = m1[i][j];
}

/* EXPORT->DMat2Mat: convert double matrix m1 to matrix m2 */
void DMat2Mat(DMatrix m1, Matrix m2)
{
   int i,j,nrows,ncols;

   nrows = NumDRows(m1);
   if (nrows != NumRows(m2))
      HError(5270,"DMat2Mat: row sizes differ %d vs %d",
             nrows,NumRows(m2));
   ncols = NumDCols(m1);
   if (ncols != NumCols(m2))
      HError(5270,"DMat2Mat: col sizes differ %d vs %d",
             ncols,NumCols(m2));
   for (i=1; i<=nrows; i++)
      for (j=1; j<=ncols; j++) 
         m2[i][j] = m1[i][j];
}

/* EXPORT->Mat2Tri: convert matrix m1 to tri matrix m2 */
void Mat2Tri (Matrix m1,  TriMat m2)
{
   int i,j,nrows,ncols;

   nrows = NumRows(m1); ncols = NumCols(m1);
   if (nrows != ncols)
      HError(5270,"Mat2Tri: source matrix not square %d vs %d",
             nrows,ncols);   
   if (ncols != TriMatSize(m2))
      HError(5270,"Mat2Tri: sizes differ %d vs %d",
             ncols,TriMatSize(m2));
   for (i=1; i<=nrows; i++)
      for (j=1; j<=i; j++) 
         m2[i][j] = m1[i][j];
}

/* EXPORT->DMat2Tri: convert matrix m1 to tri matrix m2 */
void DMat2Tri (DMatrix m1,  TriMat m2)
{
   int i,j,nrows,ncols;

   nrows = NumDRows(m1); ncols = NumDCols(m1);
   if (nrows != ncols)
      HError(5270,"DMat2Tri: source matrix not square %d vs %d",
             nrows,ncols);
   if (ncols != TriMatSize(m2))
      HError(5270,"DMat2Tri: sizes differ %d vs %d",
             ncols,TriMatSize(m2));
   for (i=1; i<=nrows; i++)
      for (j=1; j<=i; j++)
         m2[i][j] = m1[i][j];
}

void Tri2Vec (TriMat m, Vector v){
	int i,j,nrows,neles,k;
	nrows = TriMatSize(m);
	neles=VectorSize(v);
	if (nrows*(1+nrows)/2 != neles)
		HError(5270,"Tri2Vec: target vector size %d vs matrix size %d",neles,nrows*(1+nrows)/2);
	k=1;
	for (i=1; i<=nrows; i++){
		for (j=1; j<=i; j++) {
			v[k++]=m[i][j];
		}
	}
}

void Vec2Tri(Vector v, TriMat m){
	int i,j,nrows,neles,k;
	nrows = TriMatSize(m);
	neles=VectorSize(v);
	if (nrows*(1+nrows)/2 != neles)
		HError(5270,"Vec2Tri: target matrix size %d vs vector size %d",nrows*(1+nrows)/2, neles);
	k=1;
	for (i=1; i<=nrows; i++){
		for (j=1; j<=i; j++) {
			m[i][j]=v[k++];
		}
	}
}

void Mat2Vec (Matrix m, Vector v){
	int i,j,nrows,ncols,neles;
	nrows =NumRows(m);
	ncols= NumCols(m);
	neles=VectorSize(v);
	if (nrows*ncols != neles)
		HError(5270,"Tri2Vec: target vector size %d vs matrix size %d",neles,nrows*ncols);
	for (i=1; i<=nrows; i++){
		for (j=1; j<=ncols; j++) {
			v[(i-1)*ncols+j]=m[i][j];
		}
	}
}

void Vec2Mat(Vector v, Matrix m){
	int i,j,nrows,ncols,neles;
	nrows =NumRows(m);
	ncols= NumCols(m);
	neles=VectorSize(v);
	if (nrows*ncols != neles)
		HError(5270,"Vec2Mat: target matrix size %d vs vector size %d",nrows*ncols, neles);
	for (i=1; i<=nrows; i++){
		for (j=1; j<=ncols; j++) {
			m[i][j]=v[(i-1)*ncols+j];
		}
	}
}

/* EXPORT->Tri2Mat: convert tri matrix m1 to matrix m2 */
void Tri2Mat (TriMat m1, Matrix m2)
{
   int i,j,nrows,ncols;

   nrows = NumRows(m2); ncols = NumCols(m2);
   if (nrows != ncols)
      HError(5270,"Tri2Mat: target matrix not square %d vs %d",
             nrows,ncols);   
   if (ncols != TriMatSize(m1))
      HError(5270,"Tri2Mat: sizes differ %d vs %d",
             TriMatSize(m1),ncols);
   for (i=1; i<=nrows; i++)
      for (j=1; j<=i; j++) {
         m2[i][j] = m1[i][j];
         if (i!=j) m2[j][i] = m1[i][j];
      }
}

/* EXPORT->Tri2DMat: convert tri matrix m1 to dmatrix m2 */
void Tri2DMat (TriMat m1, DMatrix m2)
{
   int i,j,nrows,ncols;

   nrows = NumDRows(m2); ncols = NumDCols(m2);
   if (nrows != ncols)
      HError(5270,"Tri2DMat: target matrix not square %d vs %d",
             nrows,ncols);
   if (ncols != TriMatSize(m1))
      HError(5270,"Tri2DMat: sizes differ %d vs %d",
             TriMatSize(m1),ncols);
   for (i=1; i<=nrows; i++)
      for (j=1; j<=i; j++) {
         m2[i][j] = m1[i][j];
         if (i!=j) m2[j][i] = m1[i][j];
      }
}

/* EXPORT->DTri2DMat: convert tri matrix m1 to dmatrix m2 */
void DTri2DMat (DMatrix m1, DMatrix m2)
{
   int i,j,nrows,ncols;

   nrows = NumDRows(m2); ncols = NumDCols(m2);
   if (nrows != ncols)
      HError(5270,"DTri2DMat: target matrix not square %d vs %d",
             nrows,ncols);
   if (ncols != NumDRows(m1))
      HError(5270,"DTri2DMat: sizes differ %d vs %d",
             NumDRows(m1),ncols);
   if (ncols != NumDCols(m1))
      HError(5270,"DTri2DMat: sizes differ %d vs %d",
             NumDCols(m1),ncols);
   for (i=1; i<=nrows; i++)
      for (j=1; j<=i; j++) {
         m2[i][j] = m1[i][j];
         if (i!=j) m2[j][i] = m1[i][j];
      }
}

/* EXPORT->ReadMatrix: read matrix from source into m */
Boolean ReadMatrix(Source *src, Matrix m, Boolean binary)
{
   int i,nrows;
   
   nrows = NumRows(m);
   for (i=1; i<=nrows; i++)
      if (!ReadVector(src,m[i],binary)) 
         return FALSE;
   return TRUE;
}

/* EXPORT->ReadTriMat: read symmetric matrix in lower triangular
                       form from source into m */
Boolean ReadTriMat(Source *src, TriMat m, Boolean binary)
{
   int i,j,size;
   
   size = TriMatSize(m);
   for (j=1; j<=size; j++) {
      for (i=j; i<=size; i++)
         if (!ReadFloat(src,&(m[i][j]),1,binary))
            return FALSE;
   }
   return TRUE;
}

/* EXPORT->WriteMatrix: write matrix to f */
void WriteMatrix(FILE *f, Matrix m, Boolean binary)
{
   int i,nrows;
   
   nrows = NumRows(m);
   for (i=1; i<=nrows; i++)
      WriteVector(f,m[i],binary);
}

/* EXPORT->WriteTriMat: write symmetric matrix to stream f in
                        upper triangular form */
void WriteTriMat(FILE *f, TriMat m, Boolean binary)
{
   int i,j,size;
   
   size = TriMatSize(m);
   for (j=1; j<=size; j++) {
      for (i=j; i<=size; i++)
         WriteFloat(f,&(m[i][j]),1,binary);
      if (!binary) fputc('\n',f);
   }
}

/* Export->ShowMatrix: show the matrix m preceded by title */
void ShowMatrix(char * title,Matrix m,int maxCols,int maxRows)
{
   int i,j;
   int maxi,maxj,nrows,ncols;
   
   maxi = nrows = NumRows(m);
   if (maxi>maxRows) maxi = maxRows;
   maxj = ncols = NumCols(m);
   if (maxj>maxCols) maxj = maxCols;
   printf("%s\n",title);
   for (i=1;i<=maxi;i++) {
      printf("   ");
      for (j=1;j<=maxj;j++)
         printf("%8.2f ",m[i][j]);
      if (maxj<ncols) printf("...");
      printf("\n");
   }
   if (maxi<nrows)
      printf("   ...\n");
}

/* Export->ShowDMatrix: show the matrix m preceded by title */
void ShowDMatrix(char * title,DMatrix m,int maxCols,int maxRows)
{
   int i,j;
   int maxi,maxj,nrows,ncols;
   
   maxi = nrows = NumDRows(m);
   if (maxi>maxRows) maxi = maxRows;
   maxj = ncols = DVectorSize(m[1]);
   if (maxj>maxCols) maxj = maxCols;
   printf("%s\n",title);
   for (i=1;i<=maxi;i++) {
      printf("   ");
      for (j=1;j<=maxj;j++)
         printf("%10.4f ",m[i][j]);
      if (maxj<ncols) printf("...");
      printf("\n");
   }
   if (maxi<nrows)
      printf("   ...\n");
}

/* Export->ShowTriMat: show the matrix m preceded by title */
void ShowTriMat(char * title,TriMat m,int maxCols,int maxRows)
{
   int i,j;
   int maxi,maxj,size;
   
   size = TriMatSize(m);
   maxi = size;
   if (maxi>maxRows) maxi = maxRows;
   printf("%s\n",title);
   for (i=1;i<=maxi;i++) {
      printf("   ");
      maxj = i;
      if (maxj>maxCols) maxj = maxCols;
      for (j=1;j<=maxj;j++)
         printf("%8.2f ",m[i][j]);
      if (maxj<i) printf("...");
      printf("\n");
   }
   if (maxi<size)
      printf("   ...\n");
}

/* -------------------- Matrix Operations ---------------------- */

/* Choleski: Place lower triangular choleski factor of A in L.*/
/*           Return FALSE if matrix singular or not +definite */
Boolean Choleski(TriMat A, DMatrix L)
{
   int size,i,j,k;
   double sum;

   size = TriMatSize(A);
   for (i=1; i<=size; i++)
      for (j=1; j<=i; j++) {
         sum=A[i][j];
         for (k=1; k<j; k++)
            sum -= (L[i][k]*L[j][k]);
         if ((i==j)&&(sum<=0.0)) 
            return FALSE;
         else if (i==j)
            sum = sqrt(sum);
         else if (L[j][j]==0.0)
            return FALSE;
         else
            sum /= L[j][j];
         L[i][j] = sum;
      }
   for (i=1; i<=size; i++) 
      for (j=i+1; j<=size; j++) 
         L[i][j] = 0.0;
   return TRUE;
}

Boolean DCholeski(DMatrix A, DMatrix L)
{
   int size,i,j,k;
   double sum;

   size = NumDRows(A);
   for (i=1; i<=size; i++)
      for (j=1; j<=i; j++) {
         sum=A[i][j];
         for (k=1; k<j; k++)
            sum -= (L[i][k]*L[j][k]);
         if ((i==j)&&(sum<=0.0))
            return FALSE;
         else if (i==j)
            sum = sqrt(sum);
         else if (L[j][j]==0.0)
            return FALSE;
         else
            sum /= L[j][j];
         L[i][j] = sum;
      }
   for (i=1; i<=size; i++)
      for (j=i+1; j<=size; j++)
         L[i][j] = 0.0;
   return TRUE;
}
Boolean IsPositiveDMat(DMatrix A){
	int size=NumDRows(A);
	DMatrix L=CreateDMatrix(&gstack,size,size);
	Boolean inv=DCholeski(A, L);
	FreeDMatrix(&gstack,L);
	return inv;
}
Boolean IsPositiveTriM(TriMat A){
	int size=TriMatSize(A);
	DMatrix L=CreateDMatrix(&gstack,size,size);
	Boolean inv=Choleski(A, L);
	FreeDMatrix(&gstack,L);
	return inv;
}
static DMatrix V;
static DVector d;
static DVector e;
static int n;
	static void tred2 () {
		int i,j,k;

   /*  This is derived from the Algol procedures tred2 by
     Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
     Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
     Fortran subroutine in EISPACK.*/

      for (j = 1; j <= n; j++) {
         d[j] = V[n][j];
      }

      /* Householder reduction to tridiagonal form.*/

      for (i = n; i > 1; i--) {

         /* Scale to avoid under/overflow.*/

         double scale = 0.0;
         double h = 0.0;
         for (k = 1; k < i; k++) {
            scale = scale + fabs(d[k]);
         }
         if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 1; j < i; j++) {
               d[j] = V[i-1][j];
               V[i][j] = 0.0;
               V[j][i] = 0.0;
            }
         } else {

            /* Generate Householder vector.*/

            for (k = 1; k < i; k++) {
               d[k] /= scale;
               h += d[k] * d[k];
            }
            double f = d[i-1];
            double g = sqrt(h);
            if (f > 0) {
               g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 1; j < i; j++) {
               e[j] = 0.0;
            }

            /*Apply similarity transformation to remaining columns.*/

            for (j = 1; j < i; j++) {
               f = d[j];
               V[j][i] = f;
               g = e[j] + V[j][j] * f;
               for (k = j+1; k <= i-1; k++) {
                  g += V[k][j] * d[k];
                  e[k] += V[k][j] * f;
               }
               e[j] = g;
            }
            f = 0.0;
            for (j = 1; j < i; j++) {
               e[j] /= h;
               f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (j = 1; j < i; j++) {
               e[j] -= hh * d[j];
            }
            for (j = 1; j < i; j++) {
               f = d[j];
               g = e[j];
               for (k = j; k <= i-1; k++) {
                  V[k][j] -= (f * e[k] + g * d[k]);
               }
               d[j] = V[i-1][j];
               V[i][j] = 0.0;
            }
         }
         d[i] = h;
      }

      /* Accumulate transformations.*/

      for (i = 1; i < n; i++) {
         V[n][i] = V[i][i];
         V[i][i] = 1.0;
         double h = d[i+1];
         if (h != 0.0) {
            for (k = 1; k <= i; k++) {
               d[k] = V[k][i+1] / h;
            }
            for (j = 1; j <= i; j++) {
               double g = 0.0;
               for (k = 1; k <= i; k++) {
                  g += V[k][i+1] * V[k][j];
               }
               for (k = 1; k <= i; k++) {
                  V[k][j] -= g * d[k];
               }
            }
         }
         for (k = 1; k <= i; k++) {
            V[k][i+1] = 0.0;
         }
      }
      for (j = 1; j <= n; j++) {
         d[j] = V[n][j];
         V[n][j] = 0.0;
      }
      V[n][n] = 1.0;
      e[0] = 0.0;
   }
   static void tql2 () {
	   int i,j,k, l;
   /* This is derived from the Algol procedures tql2, by
     Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
     Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
     Fortran subroutine in EISPACK.*/

      for (i = 2; i <= n; i++) {
         e[i-1] = e[i];
      }
      e[n] = 0.0;

      double f = 0.0;
      double tst1 = 0.0;
      double eps = pow(2.0,-52.0);
      for (l = 1; l <= n; l++) {

         /* Find small subdiagonal element*/
    	 double temp=fabs(d[l]) + fabs(e[l]);
         tst1 = tst1>temp?tst1:temp;
         int m = l;
         while (m < n) {
            if (fabs(e[m]) <= eps*tst1) {
               break;
            }
            m++;
         }

         /* If m == l, d[l] is an eigenvalue,
          otherwise, iterate.*/

         if (m > l) {
            int iter = 0;
            do {
               iter = iter + 1;  /* (Could check iteration count here.)*/

               /* Compute implicit shift*/

               double g = d[l];
               double p = (d[l+1] - g) / (2.0 * e[l]);
               double r = sqrt(p*p+1.0);
               if (p < 0) {
                  r = -r;
               }
               d[l] = e[l] / (p + r);
               d[l+1] = e[l] * (p + r);
               double dl1 = d[l+1];
               double h = g - d[l];
               for (i = l+2; i <= n; i++) {
                  d[i] -= h;
               }
               f = f + h;

               /* Implicit QL transformation.*/

               p = d[m];
               double c = 1.0;
               double c2 = c;
               double c3 = c;
               double el1 = e[l+1];
               double s = 0.0;
               double s2 = 0.0;
               for (i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c * e[i];
                  h = c * p;
                  r = sqrt(p*p+e[i]*e[i]);
                  e[i+1] = s * r;
                  s = e[i] / r;
                  c = p / r;
                  p = c * d[i] - s * g;
                  d[i+1] = h + s * (c * g + s * d[i]);

                  /* Accumulate transformation.*/

                  for (k = 1; k <= n; k++) {
                     h = V[k][i+1];
                     V[k][i+1] = s * V[k][i] + c * h;
                     V[k][i] = c * V[k][i] - s * h;
                  }
               }
               p = -s * s2 * c3 * el1 * e[l] / dl1;
               e[l] = s * p;
               d[l] = c * p;

               /* Check for convergence.*/

            } while (fabs(e[l]) > eps*tst1);
         }
         d[l] = d[l] + f;
         e[l] = 0.0;
      }

      /* Sort eigenvalues and corresponding vectors.*/

      for (i = 1; i < n; i++) {
         int k = i;
         double p = d[i];
         for (j = i+1; j <= n; j++) {
            if (d[j] < p) {
               k = j;
               p = d[j];
            }
         }
         if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; j++) {
               p = V[j][i];
               V[j][i] = V[j][k];
               V[j][k] = p;
            }
         }
      }
   }

/*the columns of fV as eigenvectors, eigenvalues are already sorted from the smallest. Note: A must be +definite*/
void EigenDecomposition(TriMat A, Matrix fV, Vector fd){
	int i,j;
	n=TriMatSize(A);
	V=CreateDMatrix(&gstack,n,n);
	d=CreateDVector(&gstack,n);
	e=CreateDVector(&gstack,n);
	for (i = 1; i <= n; i++) {
	   for (j = 1; j <=i; j++) {
		  V[i][j] = A[i][j];
		  V[j][i] = V[i][j];
	   }
	}

	/* Tridiagonalize.*/
	tred2();

	/* Diagonalize.*/
	tql2();

	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++){
			fV[i][j]=V[i][j];
		}
		fd[i]=d[i];
	}
	FreeDMatrix(&gstack,V);
}

void DEigenDecomposition(DMatrix A, DMatrix fV, DVector fd){
	int i,j;
	n=NumDRows(A);
	V=CreateDMatrix(&gstack,n,n);
	d=CreateDVector(&gstack,n);
	e=CreateDVector(&gstack,n);
	for (i = 1; i <= n; i++) {
	   for (j = 1; j <=i; j++) {
		  V[i][j] = A[i][j];
		  V[j][i] = V[i][j];
	   }
	}

	/* Tridiagonalize.*/
	tred2();

	/* Diagonalize.*/
	tql2();

	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++){
			fV[i][j]=V[i][j];
		}
		fd[i]=d[i];
	}
	FreeDMatrix(&gstack,V);
}

/*Modified the diagonal elements of Hessian so that it is positive definite, return log(Det(ModH))*/
LogFloat ModDiagHessian(MemHeap* x, TriMat Hessian_from, Matrix Hessian_to, Boolean inv, float* tau_out, Boolean trace){
	float beta=0.001;
	float tau=0;
	float minaii=Hessian_from[1][1];
	int i,k;
	LogFloat lgdet;
	int size=TriMatSize(Hessian_from);
	DMatrix L=CreateDMatrix(x,size,size);
	TriMat ModH=CreateTriMat(x,size);
	CopyTriMat(Hessian_from,ModH);
	for(i=2;i<=size;i++){
		if(Hessian_from[i][i]<minaii){
			minaii=Hessian_from[i][i];
		}
	}
	if(minaii<=0)
		tau=-minaii+beta;
	for(k=1;k;k++){
		for(i=1;i<=size;i++)
			ModH[i][i]=Hessian_from[i][i]+tau;
		if(Choleski(ModH, L))
			break;
		else
			tau=2*tau>beta?2*tau:beta;
	}
	if(tau>0&&trace)
		printf("Hessian was modified by adding Identiy*tau=%f\n",tau);
	if(inv){
		lgdet=CovInvert(x, ModH, Hessian_to);
	}
	else{
		lgdet=CovDet(x, ModH);
		if(IsTriMat(Hessian_to)){
			CopyTriMat(ModH,Hessian_to);
		}
		else{
			Tri2Mat(ModH,Hessian_to);
		}
	}
	FreeDMatrix(x,L);    /* cut back stack to entry state */
	if(tau_out)
		*tau_out=tau;
	return lgdet;
}

LogDouble DModDiagHessian(DMatrix Hessian_from, DMatrix Hessian_to, Boolean inv, double *tau_out, Boolean trace){
	double beta=0.001;
	double tau=0;
	double minaii=Hessian_from[1][1];
	int i,k;
	LogDouble lgdet;
	int size=TriMatSize(Hessian_from);
	DMatrix L=CreateDMatrix(&gstack,size,size);
	DMatrix ModH=CreateDMatrix(&gstack,size,size);
	CopyDMatrix(Hessian_from,ModH);
	for(i=2;i<=size;i++){
		if(Hessian_from[i][i]<minaii){
			minaii=Hessian_from[i][i];
		}
	}
	if(minaii<=0)
		tau=-minaii+beta;
	for(k=1;k;k++){
		for(i=1;i<=size;i++)
			ModH[i][i]=Hessian_from[i][i]+tau;
		if(DCholeski(ModH, L))
			break;
		else
			tau=2*tau>beta?2*tau:beta;
	}
	if(tau>0&&trace)
		printf("Hessian was modified by adding Identiy*tau=%f\n",tau);
	if(inv){
		lgdet=DCovInvert(ModH, Hessian_to);
	}
	else{
		lgdet=DCovDet(ModH);
		CopyDMatrix(ModH,Hessian_to);
	}
	FreeDMatrix(&gstack,L);    /* cut back stack to entry state */
	if(tau_out)
		*tau_out=tau;
	return lgdet;
}

/*Modified the eigen-values of Hessian so that it is positive definite, return log(Det(ModH))*/
LogFloat ModEigHessian(TriMat Hessian_from, Matrix Hessian_to, Boolean inv, int *floor_out, Boolean trace){
	float tau=0;
	int i;
	int size=TriMatSize(Hessian_from);
	int floored=0;
	LogFloat lgdet=0;
	Boolean isTri = IsTriMat(Hessian_to);
	Matrix V=CreateMatrix(&gstack,size,size);
	Matrix Vt=CreateMatrix(&gstack,size,size);
	Matrix D=CreateMatrix(&gstack,size,size);
	Matrix temp=CreateMatrix(&gstack,size,size);
	Vector d=CreateVector(&gstack,size);
	Matrix tempX=CreateMatrix(&gstack,size,size);
	ZeroMatrix(D);
	EigenDecomposition(Hessian_from, V, d);
	for(i=1;i<=size;i++){
		if(d[i]<min_eigen){
			tau=min_eigen-d[i];
			floored++;
		}
		else{
			tau=0;
		}
		d[i]=d[i]+tau;
		D[i][i]=inv?1/d[i]:d[i];/*invert*/
		lgdet+=log(d[i]);
	}
	Transpose(V,Vt);
	MatMult(V,D,temp);
	if(isTri){
		MatMult(temp,Vt,tempX);
		Mat2Tri(tempX,Hessian_to);
	}
	else{
		MatMult(temp,Vt,Hessian_to);
	}
	FreeMatrix(&gstack,V);
	if(floored>0&&trace){
		printf("%d eigen-values floored out of %d \n",floored,size);fflush(stdout);
	}
	if(floor_out){
		*floor_out=floored;
	}
	return lgdet;
}

LogDouble DModEigHessian(DMatrix Hessian_from, DMatrix Hessian_to, Boolean inv, int *floor_out, Boolean trace){
	double tau=0;
	int i;
	int size=NumDRows(Hessian_from);
	int floored=0;
	LogDouble lgdet=0;
	DMatrix V=CreateDMatrix(&gstack,size,size);
	DMatrix Vt=CreateDMatrix(&gstack,size,size);
	DMatrix D=CreateDMatrix(&gstack,size,size);
	DMatrix temp=CreateDMatrix(&gstack,size,size);
	DVector d=CreateDVector(&gstack,size);
	ZeroDMatrix(D);
	DEigenDecomposition(Hessian_from, V, d);
	for(i=1;i<=size;i++){
		if(d[i]<min_eigen){
			tau=min_eigen-d[i];
			floored++;
		}
		else{
			tau=0;
		}
		d[i]=d[i]+tau;
		D[i][i]=inv?1/d[i]:d[i];/*invert*/
		lgdet+=log(d[i]);
	}
	DTranspose(V,Vt);
	DMatMult(V,D,temp);
	DMatMult(temp,Vt,Hessian_to);
	FreeDMatrix(&gstack,V);
	if(floored>0&&trace){
		printf("%d eigen-values floored out of %d \n",floored,size);fflush(stdout);
	}
	if(floor_out){
		*floor_out=floored;
	}
	return lgdet;
}


/*LogFloat ModifEigHessian(TriMat Hessian, Matrix invHessian){
	float tau=0;
	int i;
	int size=TriMatSize(Hessian);
	LogFloat lgdet=0;
	DMatrix L=CreateDMatrix(&gstack,size,size);
	Matrix V=CreateMatrix(&gstack,size,size);
	Matrix Vt=CreateMatrix(&gstack,size,size);
	Matrix D=CreateMatrix(&gstack,size,size);
	Matrix temp=CreateMatrix(&gstack,size,size);
	Vector d=CreateVector(&gstack,size);
	ZeroMatrix(D);
	if(!Choleski(Hessian,L)){
		EigenDecomposition(Hessian, V, d);
		for(i=1;i<=size;i++){
			if(d[i]<delta){
				tau=delta-d[i];
			}
			else{
				tau=0;
			}
			d[i]=d[i]+tau;
			D[i][i]=1/d[i];
			lgdet+=log(d[i]);
		}
		Transpose(V,Vt);
		MultiplyMM(V,D,temp);
		MultiplyMM(temp,Vt,invHessian);
	}
	else{
		lgdet=CovInvert(Hessian, invHessian);
	}
	FreeDMatrix(&gstack,L);
	return lgdet;
}*/


/* MSolve: solve Ly=e^i and L^t x = y, where e^i is a unit vector */
static void MSolve(DMatrix L, int i, DVector x, DVector y)
{
   int nr,j,k;
   double sum;
   
   nr=NumDRows(L);
   for (j=1; j<i; j++) y[j] = 0.0;  /* forward sub */
   y[i] = 1.0/L[i][i];
   for (j=i+1; j<=nr; j++){
      sum = 0.0;
      for (k=i; k<j; k++)
         sum -= L[j][k]*y[k];
      y[j] = sum /L[j][j];
   }
   x[nr] = y[nr]/L[nr][nr];         /* backward sub */
   for (j=nr-1; j>=1; j--){
      sum = y[j];
      for (k=j+1; k<=nr; k++)
         sum -= L[k][j]*x[k];
      x[j] = sum / L[j][j];
   }
}

/* EXPORT->CovInvert: puts inverse of c in invc, returns log(Det(c)) */
/*          Note that c must be positive definite */
LogFloat CovInvert(MemHeap* xx, TriMat c, Matrix invc)
{
   DMatrix l;     /* Lower Tri Choleski Matrix */
   DVector x,y;   /* for f/b substitution */
   LogFloat ldet = 0.0;
   int i,j,n;
   Boolean isTri;
   
   n = TriMatSize(c); isTri = IsTriMat(invc);
   l = CreateDMatrix(xx,n,n);
   x = CreateDVector(xx,n);
   y = CreateDVector(xx,n);
   if (Choleski(c,l)){
      for (j=1; j<=n; j++){
         MSolve(l,j,x,y);
         for (i=isTri?j:1; i<=n; i++)
            invc[i][j] = x[i];
         ldet += log(l[j][j]);
      }
   } else
      HError(5220,"CovInvert: [%f ...] not invertible",c[1][1]);
   FreeDMatrix(xx,l);    /* cut back stack to entry state */
   if(isnan(ldet))
		HError(1, "nan lgdet...");
   return 2.0*ldet;
}

/* EXPORT->CovInvert: puts inverse of c in invc, returns log(Det(c)) */
/*          Note that c must be positive definite */
LogDouble DCovInvert(TriMat c, DMatrix invc)
{
   DMatrix l;     /* Lower Tri Choleski Matrix */
   DVector x,y;   /* for f/b substitution */
   LogDouble ldet = 0.0;
   int i,j,n;
   Boolean isTri;

   n = TriMatSize(c); isTri = IsTriMat(invc);
   l = CreateDMatrix(&gstack,n,n);
   x = CreateDVector(&gstack,n);
   y = CreateDVector(&gstack,n);
   if (Choleski(c,l)){
      for (j=1; j<=n; j++){
         MSolve(l,j,x,y);
         for (i=isTri?j:1; i<=n; i++)
            invc[i][j] = x[i];
         ldet += log(l[j][j]);
      }
   } else
      HError(5220,"DCovInvert: [%f ...] not invertible",c[1][1]);
   FreeDMatrix(&gstack,l);    /* cut back stack to entry state */
   return 2.0*ldet;
}

/* EXPORT->CovInvert: puts inverse of c^=alpha*c+(1-alpha)*diag in invc^, returns log(Det(c^)) */
LogFloat CovIntp(TriMat c, Matrix invc, Boolean inv)
{
   DMatrix l;     /* Lower Tri Choleski Matrix */
   DVector x,y;   /* for f/b substitution */
   LogFloat ldet = 0.0;
   int i,j,n;
   Boolean isTri;
   float alpha;
   TriMat newc;
   n = TriMatSize(c); isTri = IsTriMat(invc);
   l = CreateDMatrix(&gstack,n,n);
   x = CreateDVector(&gstack,n);
   y = CreateDVector(&gstack,n);
   newc=CreateTriMat(&gstack,n);
   for(alpha=1;alpha>=0;alpha/=2){
	   for(i=1;i<=n;i++){
		   for(j=1;j<=i;j++){/*Convex combination of full matrix and diagonal matrix*/
			   if(i!=j){
				   newc[i][j]=alpha*c[i][j];
			   }
			   else{
				   newc[i][j]=alpha*c[i][j]+(1-alpha)*c[i][j];
			   }
		   }
	   }
	   if (Choleski(newc,l)){
		  if(inv){
			  for (j=1; j<=n; j++){
				  MSolve(l,j,x,y);
				  for (i=isTri?j:1; i<=n; i++)
					  invc[i][j] = x[i];
				  ldet += log(l[j][j]);
			  }
		  }
		  else{
			  if(isTri)
				  CopyTriMat(newc, invc);
			  else
				  CopyMatrix(newc, invc);
			  ldet=CovDet(&gstack, newc);
		  }
	      break;
	   }
   }
   FreeDMatrix(&gstack,l);    /* cut back stack to entry state */
   if(alpha<1){
 	  printf("Invertible Full Covariance matrix is made with alpha=%f\n",alpha);
 	  fflush(stdout);
   }
   return 2.0*ldet;
}
/* EXPORT->DCovIntp: puts inverse of c^=alpha*c+(1-alpha)*diag in invc^, returns log(Det(c^)) */
LogDouble DCovIntp(DMatrix c, DMatrix invc, Boolean inv)
{
   DMatrix l;     /* Lower Tri Choleski Matrix */
   DVector x,y;   /* for f/b substitution */
   LogDouble ldet = 0.0;
   int i,j,n;
   float alpha;
   DMatrix newc;
   n = NumDRows(c);
   l = CreateDMatrix(&gstack,n,n);
   x = CreateDVector(&gstack,n);
   y = CreateDVector(&gstack,n);
   newc=CreateDMatrix(&gstack,n,n);
   for(alpha=1;alpha>=0;alpha-=0.1){
	   for(i=1;i<=n;i++){
		   for(j=1;j<=n;j++){/*Convex combination of full matrix and diagonal matrix*/
			   if(i!=j){
				   newc[i][j]=alpha*c[i][j];
			   }
			   else{
				   newc[i][j]=alpha*c[i][j]+(1-alpha)*c[i][j];
			   }
		   }
	   }
	   if (DCholeski(newc,l)){
		   if(inv){
			   for (j=1; j<=n; j++){
				   MSolve(l,j,x,y);
				   for (i=1; i<=n; i++)
					   invc[i][j] = x[i];
				   ldet += log(l[j][j]);
			   }
		   }
		   else{
			   CopyDMatrix(newc, invc);
			   ldet=DCovDet(newc);
		   }
		   break;
	   }
   }
   FreeDMatrix(&gstack,l);    /* cut back stack to entry state */
   if(alpha<1){
 	  printf("Invertible Full Covariance matrix is made with alpha=%f\n",alpha);
 	  fflush(stdout);
   }
   return 2.0*ldet;
}

/* EXPORT->CovDet: Returns log(Det(c)), c must be positive definite */
LogFloat CovDet(MemHeap* x, TriMat c)
{
   DMatrix l;  /* Lower Tri Choleski Matrix */
   LogFloat ldet = 0.0;
   int j,n;
   
   n = TriMatSize(c);
   l = CreateDMatrix(x,n,n);
   if (Choleski(c,l)){
      for (j=1; j<=n; j++)
         ldet += log(l[j][j]);
   } else
      HError(5220,"CovDet: [%f ...] not invertible",c[1][1]);
   FreeDMatrix(x,l);
   return 2.0*ldet;
}

LogDouble DCovDet(DMatrix c)
{
   DMatrix l;  /* Lower Tri Choleski Matrix */
   LogDouble ldet = 0.0;
   int j,n;

   n = NumDRows(c);
   l = CreateDMatrix(&gstack,n,n);
   if (DCholeski(c,l)){
      for (j=1; j<=n; j++)
         ldet += log(l[j][j]);
   } else
      HError(5220,"CovDet: [%f ...] not invertible",c[1][1]);
   FreeDMatrix(&gstack,l);
   return 2.0*ldet;
}

/* Quadratic prod of a full square matrix C and an arbitry full matrix transform A */
void LinTranQuaProd(Matrix Prod, Matrix A, Matrix C)
{
   int i,j,k;
   float tempElem;
   Matrix tempMatrix_A_mult_C;
   
   if (NumRows(C) != NumCols(C)){
      HError(999,"HMath: LinTranQuaProd: Matrix C is not square!\n");
   }
   else {
      tempMatrix_A_mult_C = CreateMatrix(&gstack,NumRows(A),NumCols(C));
      ZeroMatrix(tempMatrix_A_mult_C);
      
      for (i=1;i<=NumRows(tempMatrix_A_mult_C);i++){
         for (j=1;j<=NumCols(tempMatrix_A_mult_C);j++){
            tempElem = 0.0;
            for (k=1;k<=NumCols(A);k++){
               tempElem += A[i][k]*C[j][k];
            }
            tempMatrix_A_mult_C[i][j] = tempElem;
         }
      }
      
      for (i=1;i<=NumRows(Prod);i++){
         for (j=1;j<=i;j++){
            tempElem = 0.0;
            for (k=1;k<=NumCols(tempMatrix_A_mult_C);k++){
               tempElem += tempMatrix_A_mult_C[i][k]*A[j][k];
            }
            Prod[i][j] = tempElem;
         }
      }
      
      for (i=1;i<=NumRows(Prod);i++){
         for (j=1;j<i;j++){
            Prod[j][i] = Prod[i][j];
         }
      }    
      
      FreeMatrix(&gstack,tempMatrix_A_mult_C);
   }
}


/* ------------- Singular Value Decomposition --------------- */
/**************************************************************************
 **
 ** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
 **
 **                          Meschach Library
 ** 
 ** This Meschach Library is provided "as is" without any express 
 ** or implied warranty of any kind with respect to this software. 
 ** In particular the authors shall not be liable for any direct, 
 ** indirect, special, incidental or consequential damages arising 
 ** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
**  Modifications made to conform with HTK formats by 
**  Daniel Kershaw, Entropic Ltd, Cambridge, England.
**
***************************************************************************/
#define MACHEPS 2.22045e-16
#define FZERO 1.0e-6
#define sgn(x)  ((x) >= 0 ? 1 : -1)
#define minab(a,b) ((a) > (b) ? (b) : (a))
#define MAX_STACK       100

/* Givens -- returns c,s parameters for Givens rotation to
   eliminate y in the vector [ x y ]' */
static void Givens(double x, double y, double *c, double *s)
{
   double norm;
  
   norm = sqrt(x*x+y*y);
   if ( norm == 0.0 ) {
      *c = 1.0;
      *s = 0.0;       
   }       /* identity */
   else {
      *c = x/norm;
      *s = y/norm;
   }
}


/* RotRows -- premultiply mat by givens rotation described by c,s */
static void RotRows(DMatrix M, int i, int k, 
                    double c, double s)
{
   int   j, n;
   double temp;
  
   n = NumDRows(M);
  
   if (i > n || k > n)
      HError(1, "RotRows: Index tooo big i=%d k=%d\n", i, k);
  
  
   for ( j=1; j<=n; j++ ) {
      temp = c*M[i][j] + s*M[k][j];
      M[k][j] = -s*M[i][j] + c*M[k][j];
      M[i][j] = temp;
   }
  
}

/* FixSVD -- fix minor details about SVD make singular values non-negative
   -- sort singular values in decreasing order */
static void FixSVD(DVector d, DMatrix U, DMatrix V)
{
  
   int  i, j, n;
  
   n = DVectorSize(d);


   /* make singular values non-negative */
   for (i = 1; i <= n; i++) {
      if ( d[i] < 0.0 ) {
         d[i] = - d[i];
         for ( j = 1; j <= NumDRows(U); j++ )
            U[i][j] = - U[i][j];
      }
   }

   return;

#if 0 /* #### ge: what is this code after return supposed to do here? */
   {
      int  k, l, r, stack[MAX_STACK], sp;
      double tmp, v;

   /* sort singular values */
   sp = -1;
   l = 1;       
   r = n;
   for ( ; ; ) {
      while ( r >= l ) {
         /* i = partition(d->ve,l,r) */
         v = d[r];
         i = l-1;           
         j = r;
         for ( ; ; ) {
            /* inequalities are "backwards" for **decreasing** order */
            while ( d[++i] > v );
            while ( d[--j] < v );
            if ( i >= j )
               break;
            /* swap entries in d->ve */
            tmp = d[i];   
            d[i] = d[j];
            d[j] = tmp;
            /* swap rows of U & V as well */
            for ( k = 1; k <= DVectorSize(U[1]); k++ ) {
               tmp = U[i][k];
               U[i][k] = U[j][k];
               U[j][k] = tmp;
            }
            for ( k = 1; k <= DVectorSize(V[1]); k++ ) {
               tmp = V[i][k];
               V[i][k] = V[j][k];
               V[j][k] = tmp;
            }
         }
         tmp = d[i];
         d[i] = d[r];
         d[r] = tmp;
         for ( k = 1; k <= DVectorSize(U[1]); k++ ) {
            tmp = U[i][k];
            U[i][k] = U[r][k];
            U[r][k] = tmp;
         }
         for ( k = 1; k <= DVectorSize(V[1]); k++ ) {
            tmp = V[i][k];
            V[i][k] = V[r][k];
            V[r][k] = tmp;
         }
         /* end i = partition(...) */
         if ( i - l > r - i ) {
            stack[++sp] = l;    
            stack[++sp] = i-1;
            l = i+1;    
         }
         else {
            stack[++sp] = i+1;
            stack[++sp] = r;
            r = i-1;    
         }
      }
      if ( sp < 0 )
         break;
      r = stack[sp--];
      l = stack[sp--];
   }
   }
#endif
}

/* BiSvd -- svd of a bidiagonal m x n matrix represented by d (diagonal) and
   f (super-diagonals) */
static void BiSVD(DVector d, DVector f, DMatrix U, DMatrix V)
{
   int i, j, n;
   int i_min, i_max, split;
   double c, s, shift, size, z;
   double d_tmp, diff, t11, t12, t22;

   if ( ! d || ! f )
      HError(1,"BiSVD: Vectors are null!");
   if ( DVectorSize(d) != DVectorSize(f) + 1 )
      HError(1, "BiSVD: Error with the vector sizes!");

   n = DVectorSize(d);

   if ( ( U && DVectorSize(U[1]) < n ) || ( V && NumDRows(V) < n ) )
      HError(1, "BiSVD: Error Matrix sizes!");
   if ( ( U && NumDRows(U) != DVectorSize(U[1])) ||
        ( V && NumDRows(V) != DVectorSize(V[1])) )
      HError(1, "BiSVD: One of the matrices must be square");


   if ( n == 1 )
      return;

   s = 0.0;
   for ( i = 1; i <= n; i++)
      s += d[i]*d[i];
   size = sqrt(s);
   s = 0.0;
   for ( i = 1; i < n; i++)
      s += f[i]*f[i];
   size += sqrt(s);
   s = 0.0;

   i_min = 1;
   while ( i_min <= n ) {   /* outer while loop */
      /* find i_max to suit;
         submatrix i_min..i_max should be irreducible */
      i_max = n;
      for ( i = i_min; i < n; i++ )
         if ( d[i] == 0.0 || f[i] == 0.0 ) {
            i_max = i;
            if ( f[i] != 0.0 ) {
               /* have to ``chase'' f[i] element out of matrix */
               z = f[i];
               f[i] = 0.0;
               for ( j = i; j < n && z != 0.0; j++ ) {
                  Givens(d[j+1],z, &c, &s);
                  s = -s;
                  d[j+1] =  c*d[j+1] - s*z;
                  if ( j+1 < n ) {
                     z      = s*f[j+1];
                     f[j+1] = c*f[j+1];
                  }
                  RotRows(U,i,j+1,c,s);
               }
            }
            break;
         }


      if ( i_max <= i_min ) {
         i_min = i_max + 1;
         continue;
      }

      split = FALSE;
      while ( ! split ) {
         /* compute shift */
         t11 = d[i_max-1]*d[i_max-1] +
            (i_max > i_min+1 ? f[i_max-2]*f[i_max-2] : 0.0);
         t12 = d[i_max-1]*f[i_max-1];
         t22 = d[i_max]*d[i_max] + f[i_max-1]*f[i_max-1];
         /* use e-val of [[t11,t12],[t12,t22]] matrix
            closest to t22 */
         diff = (t11-t22)/2;
         shift = t22 - t12*t12/(diff +
                                sgn(diff)*sqrt(diff*diff+t12*t12));

         /* initial Givens' rotation */
         Givens(d[i_min]*d[i_min]-shift,
                d[i_min]*f[i_min], &c, &s);

         /* do initial Givens' rotations */
         d_tmp      = c*d[i_min] + s*f[i_min];
         f[i_min]   = c*f[i_min] - s*d[i_min];
         d[i_min]   = d_tmp;
         z          = s*d[i_min+1];
         d[i_min+1] = c*d[i_min+1];
         RotRows(V,i_min,i_min+1,c,s);

         /* 2nd Givens' rotation */
         Givens(d[i_min],z, &c, &s);
         d[i_min]   = c*d[i_min] + s*z;
         d_tmp      = c*d[i_min+1] - s*f[i_min];
         f[i_min]   = s*d[i_min+1] + c*f[i_min];
         d[i_min+1] = d_tmp;
         if ( i_min+1 < i_max ) {
            z          = s*f[i_min+1];
            f[i_min+1] = c*f[i_min+1];
         }
         RotRows(U,i_min,i_min+1,c,s);

         for ( i = i_min+1; i < i_max; i++ ) {
            /* get Givens' rotation for zeroing z */
            Givens(f[i-1],z, &c, &s);
            f[i-1] = c*f[i-1] + s*z;
            d_tmp  = c*d[i] + s*f[i];
            f[i]   = c*f[i] - s*d[i];
            d[i]   = d_tmp;
            z      = s*d[i+1];
            d[i+1] = c*d[i+1];
            RotRows(V,i,i+1,c,s);

            /* get 2nd Givens' rotation */
            Givens(d[i],z, &c, &s);
            d[i]   = c*d[i] + s*z;
            d_tmp  = c*d[i+1] - s*f[i];
            f[i]   = c*f[i] + s*d[i+1];
            d[i+1] = d_tmp;
            if ( i+1 < i_max ) {
               z      = s*f[i+1];
               f[i+1] = c*f[i+1];
            }
            RotRows(U,i,i+1,c,s);
         }
         /* should matrix be split? */
         for ( i = i_min; i < i_max; i++ )
            if ( fabs(f[i]) <
                 MACHEPS*(fabs(d[i])+fabs(d[i+1])) )
               {
                  split = TRUE;
                  f[i] = 0.0;
               }
            else if ( fabs(d[i]) < MACHEPS*size )
               {
                  split = TRUE;
                  d[i] = 0.0;
               }
      }
   }
}

/* HholdVec -- calulates Householder vector to eliminate all entries after the
   i0 entry of the vector vec. It is returned as out. May be in-situ */
static void HholdVec(DVector tmp, int i0, int size,
                     double *beta, double *newval)
{
   int i;
   double norm = 0.0;

   for (i = i0; i <= size; i++) {
      norm += tmp[i]*tmp[i];
   }
   norm = sqrt(norm);

   if ( norm <= 0.0 ) {
      *beta = 0.0;
   }
   else {
      *beta = 1.0/(norm * (norm+fabs(tmp[i0])));
      if ( tmp[i0] > 0.0 )
         *newval = -norm;
      else
         *newval = norm;
      tmp[i0] -= *newval;
   }

}


/* HholdTrRows -- transform a matrix by a Householder vector by rows
   starting at row i0 from column j0 -- in-situ */
static void HholdTrRows(DMatrix M, int i0, int j0, DVector hh, double beta)
{
   double ip, scale;
   int i, j;
   int m,n;

   m = NumDRows(M);
   n = DVectorSize(M[1]);

   if ( M==NULL || hh==NULL )
      HError(1,"HholdTrRows: matrix or vector is NULL!");
   if ( DVectorSize(hh) != n )
      HError(1,"HholdTrRows: hh vector size must = number of M columns");
   if ( i0 > m+1 || j0 > n+1 )
      HError(1,"HholdTrRows: Bounds matrix/vec size error i=%d j=%d m=%d n=%d",
             i0, j0, m, n);
  
   if ( beta != 0.0 ) {
      /* for each row ... */
      for ( i = i0; i <= m; i++ )
         {  
            /* compute inner product */
            /* ip = __ip__(&(M->me[i][j0]),&(hh->ve[j0]),(int)(M->n-j0));*/
            ip = 0.0;
            for ( j = j0; j <= n; j++ )
               ip += M[i][j]*hh[j];
            scale = beta*ip;
            if ( scale == 0.0 )
               continue;
            /* __mltadd__(&(M->me[i][j0]),&(hh->ve[j0]),-scale,
               (int)(M->n-j0)); */
            for ( j = j0; j <= n; j++ )
               M[i][j] -= scale*hh[j];
         }
   }
}

/* HholdTrCols -- transform a matrix by a Householder vector by columns
   starting at row i0 from column j0 -- in-situ */
static void HholdTrCols(DMatrix M, int i0, int j0, 
                        DVector hh, double beta, DVector w)
{
   int i, j;
   int n;

   n = NumDRows(M);

   if ( M==NULL || hh==NULL )
      HError(1,"HholdTrCols: matrix or vector is NULL!");
   if ( DVectorSize(hh) != n )
      HError(1,"HholdTrCols: hh vector size must = number of M columns");
   if ( i0 > n+1 || j0 > n+1 )
      HError(1,"HholdTrCols: Bounds matrix/vec size error i=%d j=%d n=%d",
             i0, j0, n);

   ZeroDVector(w);

   if ( beta != 0.0 ) {

      for ( i = i0; i <= n; i++ )
         if ( hh[i] != 0.0 )
            for ( j = j0; j <= n; j++ )
               w[j] += M[i][j]*hh[i];

      for ( i = i0; i <= n; i++ )
         if ( hh[i] != 0.0 )
            for ( j = j0; j <= n; j++ )
               M[i][j] -= w[j]*beta*hh[i];

   }
}


/* copy a row from a matrix into  a vector */
static void CopyDRow(DMatrix M, int k, DVector v) 
{
   int i, size;
   DVector w;

   if (v == NULL)
      HError(1, "CopyDRow: Vector is NULL");

   size = DVectorSize(v);
   w = M[k];

   for (i = 1; i <= size; i++)
      v[i] = w[i];
}

/* copy a column from a matrix into  a vector */
static void CopyDColumn(DMatrix M, int k, DVector v) 
{
   int i, size;

   if (v == NULL)
      HError(1, "CopyDColumn: Vector is NULL");

   size = DVectorSize(v);
   for (i = 1; i <= size; i++)
      v[i] = M[i][k];
}

/* BiFactor -- perform preliminary factorisation for bisvd
   -- updates U and/or V, which ever is not NULL */
static void BiFactor(DMatrix A, DMatrix U, DMatrix V)
{
   int n, k;
   DVector tmp1, tmp2, tmp3;
   double beta;

   n = NumDRows(A);

   tmp1 = CreateDVector(&gstack, n);
   tmp2 = CreateDVector(&gstack, n);
   tmp3 = CreateDVector(&gstack, n);

   for ( k = 1; k <= n; k++ ) {
      CopyDColumn(A,k,tmp1);
      HholdVec(tmp1,k,n,&beta,&(A[k][k]));
      HholdTrCols(A,k,k+1,tmp1,beta,tmp3);
      if ( U )
         HholdTrCols(U,k,1,tmp1,beta,tmp3);
      if ( k+1 > n )
         continue;
      CopyDRow(A,k,tmp2);
      HholdVec(tmp2,k+1,n,&beta,&(A[k][k+1]));
      HholdTrRows(A,k+1,k+1,tmp2,beta);
      if ( V )
         HholdTrCols(V,k+1,1,tmp2,beta,tmp3);
   }

   FreeDVector(&gstack, tmp1);
}

/* mat_id -- set A to being closest to identity matrix as possible
        -- i.e. A[i][j] == 1 if i == j and 0 otherwise */
static void InitIdentity(DMatrix A) 
{
   int     i, size;
  
   ZeroDMatrix(A);
   size = minab(NumDRows(A), DVectorSize(A[1]));
   for ( i = 1; i <= size; i++ )
      A[i][i] = 1.0;
}

/* EXPORT->SVD: Calculate the decompostion of matrix A.
   NOTE: on return that U and V hold U' and V' respectively! */
void SVD(DMatrix A, DMatrix U, DMatrix V, DVector d)
{
   DVector f=NULL;
   int i, n;
   DMatrix A_tmp;

   /* do initial size checks! */
   if ( A == NULL )
      HError(1,"svd: Matrix A is null");

   n = NumDRows(A);

   if (U == NULL || V == NULL || d == NULL)
      HError(1, "SVD: The svd matrices and vector must be initialised b4 call");
 
   A_tmp = CreateDMatrix(&gstack, n, n);
   CopyDMatrix(A, A_tmp);
   InitIdentity(U);
   InitIdentity(V);
   f = CreateDVector(&gstack,n-1);

   BiFactor(A_tmp,U,V);
   for ( i = 1; i <= n; i++ ) {
      d[i] = A_tmp[i][i];
      if ( i+1 <= n )
         f[i] = A_tmp[i][i+1];
   }

   BiSVD(d,f,U,V);
   FixSVD(d,U,V);
   FreeDMatrix(&gstack, A_tmp);
}

/* EXPORT->InvSVD: Inverted Singular Value Decomposition (calls SVD)
   and inverse of A is returned in Result */
void InvSVD(DMatrix A, DMatrix U, DVector W, DMatrix V, DMatrix Result)
{
   int m, n, i, j, k;
   double wmax, wmin;
   Boolean isSmall = FALSE;
   DMatrix tmp1;

   m = NumDRows(U);
   n = DVectorSize(U[1]);

   if (m != n)
      HError(1, "InvSVD: Matrix inversion only for symmetric matrices!\n");

   SVD(A, U, V, W);
   /* NOTE U and V actually now hold U' and V' ! */

   tmp1 = CreateDMatrix(&gstack,m, n);

   wmax = 0.0;
   for (k = 1; k <= n; k ++)
      if (W[k] > wmax)
         wmax = W[k];
   wmin = wmax * 1.0e-8;
   for (k = 1; k <= n; k ++)
      if (W[k] < wmin) {
         /* A component of the diag matrix 'w' of the SVD of 'a'
            was smaller than 1.0e-6 and consequently set to zero. */
         if (trace>0) {
            printf("%d (%e) ", k, W[k]); 
            fflush(stdout);
         }
         W[k] = 0.0;
         isSmall = TRUE;
      }
   if (trace>0 && isSmall) {
      printf("\n"); 
      fflush(stdout);
   }
   /* tmp1 will be the product of matrix v and the diagonal 
      matrix of singular values stored in vector w. tmp1 is then
      multiplied by the transpose of matrix u to produce the 
      inverse which is returned */
   for (j = 1; j <= m; j++)
      for (k = 1; k <= n; k ++)
         if (W[k] > 0.0)
            /* Only non-zero elements of diag matrix w are 
               used to compute the inverse. */
            tmp1[j][k] = V[k][j] / W[k];
         else
            tmp1[j][k] = 0.0;

   ZeroDMatrix(Result);
   for (i=1;i<=m;i++)
      for (j=1;j<=m;j++)
         for (k=1;k<=n;k++)
            Result[i][j] += tmp1[i][k] * U[k][j];
   FreeDMatrix(&gstack,tmp1);
}

/* LUDecompose: perform LU decomposition on Matrix a, the permutation
       of the rows is returned in perm and sign is returned as +/-1
       depending on whether there was an even/odd number of row 
       interchanges */
Boolean LUDecompose(Matrix a, int *perm, int *sign)
{
   int i,imax,j,k,n;
   double scale,sum,xx,yy;
   Vector vv,tmp;
   
   n = NumRows(a); imax = 0;
   vv = CreateVector(&gstack,n);
   *sign = 1;
   for (i=1; i<=n; i++) {
      scale = 0.0;
      for (j=1; j<=n; j++)
         if ((xx = fabs(a[i][j])) > scale )
            scale = xx;
      if (scale == 0.0) {
         HError(-1,"LUecompose: Matrix is Singular");
	 return(FALSE);
      }
      vv[i] = 1.0/scale;
   }
   for (j=1; j<=n; j++) {
      for (i=1; i<j; i++) {
         sum = a[i][j];
         for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      scale=0.0;
      for (i=j; i<=n; i++) {
         sum = a[i][j];
         for (k=1; k<j; k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (yy=vv[i]*fabs(sum)) >=scale) {
            scale = yy; imax = i;
         }
      }
      if (j != imax) {
         tmp = a[imax]; a[imax] = a[j]; a[j] = tmp;
         *sign = -(*sign);
         vv[imax]=vv[j];
      }
      perm[j]=imax;
      if (a[j][j] == 0.0) {
         HError(-1,"LUDecompose: Matrix is Singular");
	 return(FALSE);
      }
      if (j != n) {
         yy = 1.0/a[j][j];
         for (i=j+1; i<=n;i++) a[i][j] *= yy;
      }
   }
   FreeVector(&gstack,vv);
   return(TRUE);
}


/* EXPORT->MatDet: determinant of a matrix */
float MatDet(Matrix c)
{
   Matrix a;
   float det;
   int n,perm[1600],i,sign;
   
   n=NumRows(c);
   a=CreateMatrix(&gstack,n,n);
   CopyMatrix(c,a);                /* Make a copy of c */
   LUDecompose(a,perm,&sign);      /* Do LU Decomposition */
   det = sign;                     /* Calc Det(c) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   FreeMatrix(&gstack,a);
   return det;
}

/* mat_id -- set A to being closest to identity matrix as possible
        -- i.e. A[i][j] == 1 if i == j and 0 otherwise */
void Identity(Matrix A)
{
   int     i, j;

   ZeroMatrix(A);
   int row=NumRows(A);
   int col=NumCols(A);
   if(row!=col){
	   HError(-1,"Identity: A must be square");
   }
   for ( i = 1; i <= row; i++ ){
	   for( j=1; j<=col;j++){
		   if(i!=j){
			   A[i][j] = 0;
		   }
		   else{
			   A[i][j] = 1.0;
		   }
	   }
   }
}
void DIdentity(DMatrix A)
{
   int     i, j;

   ZeroDMatrix(A);
   int row=NumDRows(A);
   int col=NumDCols(A);
   if(row!=col){
	   HError(-1,"DIdentity: A must be square");
   }
   for ( i = 1; i <= row; i++ ){
	   for( j=1; j<=col;j++){
		   if(i!=j){
			   A[i][j] = 0;
		   }
		   else{
			   A[i][j] = 1.0;
		   }
	   }
   }
}
void DiagMatrix(Matrix A)
{
   int     i, j;
   int row=NumRows(A);
   int col=NumCols(A);
   if(row!=col){
	   HError(-1,"Identity: A must be square");
   }
   for ( i = 1; i <= row; i++ ){
	   for( j=1; j<=col;j++){
		   if(i!=j){
			   A[i][j] = 0;
		   }
	   }
   }
}

void DiagTriMat(TriMat A)
{
   int     i, j;
   int size=TriMatSize(A);
   for ( i = 1; i <= size; i++ ){
	   for( j=1; j<=i;j++){
		   if(i!=j){
			   A[i][j] = 0;
		   }
	   }
   }
}

float P_norm(Vector v, int p){
	int vs=VectorSize(v);
	int i;
	float norm=p==0?-1E10:0;
	for(i=1;i<=vs;i++){
		float temp=fabs(v[i]);
		if(p==0){
			norm=(norm>=temp)?norm:temp;
		}
		else{
			norm+=pow(temp,p);
		}
	}
	if(p>=1)
		norm=pow(norm,1.0/p);
	return norm;
}

float XP_norm(Matrix x, int p){
	int i,j;
	Boolean isTri=IsTriMat(x);
	int row,col;
	if(isTri){
		row=col=TriMatSize(x);
	}
	else{
		row=NumRows(x);
		col=NumCols(x);
	}
	float norm=p==0?-1E10:0;
	for(i=1;i<=row;i++){
		for(j=1;j<=(isTri?i:col);j++){
			float temp=fabs(x[i][j]);
			if(p==0){
				norm=(norm>=temp)?norm:temp;
			}
			else{
				norm+=pow(temp,p);
			}
		}
	}
	if(p>=1)
		norm=pow(norm,1.0/p);
	return norm;
}

Boolean DotProduct(Vector v1,Vector v2, float *res){
	int s1=VectorSize(v1);
	int s2=VectorSize(v2);
	if(s1!=s2){
		HError(-1,"DotProduct: v1 and v2 doesn't match");
		return FALSE;
	}
	int i;
	*res=0;
	for(i=1;i<=s1;i++){
		*res+=v1[i]*v2[i];
	}
	return TRUE;
}

Boolean MulDotProduct(Vector v1,Vector v2, short *mul, float *res){
	int s1=VectorSize(v1);
	int s2=VectorSize(v2);
	if(s1!=s2){
		HError(-1,"MulDotProduct: v1 and v2 doesn't match");
		return FALSE;
	}
	int s,i,S=mul[0],offset=0;
	res[0]=1;
	for(s=1;s<=S;s++){
		res[s]=0;
		for(i=1;i<=mul[s];i++){
			res[s]+=v1[offset+i]*v2[offset+i];
		}
		res[0]*=res[s];/*plus for mixture model (init to 0), multiply for product model (init to 1)*/
		offset+=mul[s];
	}
	return TRUE;
}

Boolean MulHardDotProduct(Vector v1,Vector v2, short *mul, float *res){
	int s1=VectorSize(v1);
	int s2=VectorSize(v2);
	if(s1!=s2){
		HError(-1,"MulDotProduct: v1 and v2 doesn't match");
		return FALSE;
	}
	int s,i,S=mul[0],offset=0,maxIdx2;
	float maxV2;
	res[0]=1;
	for(s=1;s<=S;s++){
		maxV2=0;
		maxIdx2=1;
		for(i=1;i<=mul[s];i++){
			if(v2[offset+i]>maxV2){
				maxV2=v2[offset+i];
				maxIdx2=offset+i;
			}
		}
		res[s]=v1[maxIdx2];
		res[0]*=res[s];/*plus for mixture model (init to 0), multiply for product model (init to 1)*/
		offset+=mul[s];
	}
	return TRUE;
}


Boolean VecAdd(Vector v1, float scale1, Vector v2, float scale2, Vector res){
	int s1=VectorSize(v1);
	int s2=VectorSize(v2);
	int s=VectorSize(res);
	if(s1!=s2){
		HError(-1,"VecAdd: v1 and v2 doesn't match");
		return FALSE;
	}
	int i;
	for(i=1;i<=s1;i++){
		res[i]=v1[i]*scale1+v2[i]*scale2;
	}
	return TRUE;
}
Boolean DVecAdd(DVector v1, double scale1, DVector v2, double scale2, DVector res){
	int s1=DVectorSize(v1);
	int s2=DVectorSize(v2);
	int s=DVectorSize(res);
	if(s1!=s2){
		HError(-1,"DVecAdd: v1 and v2 doesn't match");
		return FALSE;
	}
	int i;
	for(i=1;i<=s1;i++){
		res[i]=v1[i]*scale1+v2[i]*scale2;
	}
	return TRUE;
}


Boolean DDotProduct(DVector v1,DVector v2, double *res){
	int s1=DVectorSize(v1);
	int s2=DVectorSize(v2);
	if(s1!=s2){
		HError(-1,"DotProduct: v1 and v2 doesn't match");
		return FALSE;
	}
	int i;
	*res=0;
	for(i=1;i<=s1;i++){
		*res+=v1[i]*v2[i];
	}
	return TRUE;
}

Boolean QuaProduct(Vector v1,Vector v2, Matrix res){
	int s1=VectorSize(v1);
	int s2=VectorSize(v2);
	int row=NumRows(res);
	int col=NumCols(res);
	if(row!=s1||col!=s2){
		HError(-1,"QuaProduct: res doesn't match with v1*v2'");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			res[i][j]=v1[i]*v2[j];
		}
	}
	return TRUE;
}
float Trace(Matrix m1){
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int i,j;
	float trace=0;
	for(i=1;i<=row1&&i<=col1;i++){
		trace+=m1[i][i];
	}
	return trace;
}
double DTrace(DMatrix m1){
	int row1=NumDRows(m1);
	int col1=NumDCols(m1);
	int i,j;
	double trace=0;
	for(i=1;i<=row1&&i<=col1;i++){
		trace+=m1[i][i];
	}
	return trace;
}
Boolean Transpose(Matrix m1,Matrix m1t){
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int row2=NumRows(m1t);
	int col2=NumCols(m1t);
	if(row1!=col2||col1!=row2){
		HError(-1,"Transpose: m1t doesn't match with m1");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row2;i++){
		for(j=1;j<=col2;j++){
			m1t[i][j]=m1[j][i];
		}
	}
	return TRUE;
}
Boolean DTranspose(DMatrix m1,DMatrix m1t){
	int row1=NumDRows(m1);
	int col1=NumDCols(m1);
	int row2=NumDRows(m1t);
	int col2=NumDCols(m1t);
	if(row1!=col2||col1!=row2){
		HError(-1,"Transpose: m1t doesn't match with m1");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row2;i++){
		for(j=1;j<=col2;j++){
			m1t[i][j]=m1[j][i];
		}
	}
	return TRUE;
}
Boolean MatColVector(Matrix m, Vector v, int j){
	int row=NumRows(m);
	int vs=VectorSize(v);
	if(vs!=row){
		HError(-1,"MatColVector: v doesn't match with m[*][j]");
		return FALSE;
	}
	int i;
	for(i=1;i<=row;i++){
		v[i]=m[i][j];
	}
	return TRUE;
}
Boolean MatMult(Matrix m1, Matrix m2, Matrix res){
	int row=NumRows(res);
	int col=NumCols(res);
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int row2=NumRows(m2);
	int col2=NumCols(m2);
	if(col1!=row2||row1!=row||col2!=col){
		HError(-1,"MatMult: res doesn't match with m1*m2");
		return FALSE;
	}
	int i,j;
	Matrix m2t=CreateMatrix(&gstack,col2,row2);
	Transpose(m2,m2t);
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			Vector v1=m1[i];
			Vector v2=m2t[j];
			DotProduct(v1,v2,&res[i][j]);
		}
	}
	FreeMatrix(&gstack,m2t);
	return TRUE;
}
Boolean SqMatMult(Matrix m1, Matrix m2, Matrix res){
	int row=NumRows(res);
	int col=NumCols(res);
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int row2=NumRows(m2);
	int col2=NumCols(m2);
	if(col1!=row2||row1!=row||col2!=col){
		HError(-1,"MatMult: res doesn't match with m1*m2");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			Vector v1=m1[i];
			Vector v2=m2[j];
			DotProduct(v1,v2,&res[i][j]);
		}
	}
	return TRUE;
}
Boolean MatAdd(Matrix m1, float scale1, Matrix m2, float scale2, Matrix res){
	int row=NumRows(res);
	int col=NumCols(res);
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int row2=NumRows(m2);
	int col2=NumCols(m2);
	if(col1!=col2||row1!=row2||col1!=col||row1!=row){
		HError(-1,"MatAdd: res doesn't match with m1*m2");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			res[i][j]=m1[i][j]*scale1+m2[i][j]*scale2;
		}
	}
	return TRUE;
}
Boolean DMatAdd(DMatrix m1, double scale1, DMatrix m2, double scale2, DMatrix res){
	int row=NumDRows(res);
	int col=NumDCols(res);
	int row1=NumDRows(m1);
	int col1=NumDCols(m1);
	int row2=NumDRows(m2);
	int col2=NumDCols(m2);
	if(col1!=col2||row1!=row2||col1!=col||row1!=row){
		HError(-1,"DMatAdd: res doesn't match with m1*m2");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			res[i][j]=m1[i][j]*scale1+m2[i][j]*scale2;
		}
	}
	return TRUE;
}

Boolean DMatMult(DMatrix m1, DMatrix m2, DMatrix res){
	int row=NumDRows(res);
	int col=NumDCols(res);
	int row1=NumDRows(m1);
	int col1=NumDCols(m1);
	int row2=NumDRows(m2);
	int col2=NumDCols(m2);
	if(col1!=row2||row1!=row||col2!=col){
		HError(-1,"MatMult: res doesn't match with m1*m2");
		return FALSE;
	}
	int i,j,k;
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			res[i][j]=0;
			for(k=1;k<=col1;k++){
				res[i][j]+=m1[i][k]*m2[k][j];
			}
		}
	}
	return TRUE;
}



float EuclMatDist(Matrix m1, Matrix m2){
	int row1=NumRows(m1);
	int col1=NumCols(m1);
	int row2=NumRows(m2);
	int col2=NumCols(m2);
	if(row1!=row2||col1!=col2){
		HError(-1,"EuclMatDist: m1 doesn't match with m2");
		return FALSE;
	}
	int i,j;
	float res=0;
	for(i=1;i<=row1;i++){
		for(j=1;j<=col1;j++){
			float dif=m1[i][j]-m2[i][j];
			res+=dif*dif;
		}
	}
	return sqrt(res);
}
Boolean MatScale(Matrix m,float s,Matrix res){
	int row=NumRows(m);
	int col=NumCols(m);
	int row2=NumRows(m);
	int col2=NumCols(m);
	if(col!=col2||row!=row2){
		HError(-1,"MatScale: res doesn't match with m*s");
		return FALSE;
	}
	int i;
	for(i=1;i<=row;i++){
		VecScale(m[i],s,res[i]);
	}
	return TRUE;
}
Boolean VecScale(Vector v,float s,Vector res){
	int vs=VectorSize(v);
	int vs2=VectorSize(res);
	if(vs!=vs2){
		HError(-1,"VecScale: res doesn't match with v*s");
		return FALSE;
	}
	int i;
	for(i=1;i<=vs;i++){
		res[i]=s*v[i];
	}
	return TRUE;
}
Boolean MatInterp(Matrix m,Vector v,Vector res){
	int row=NumRows(m);
	int col=NumCols(m);
	int vs=VectorSize(v);
	int rs=VectorSize(res);
	if(col!=vs||row!=rs){
		HError(-1,"MatInterp: res doesn't match with m*v");
		return FALSE;
	}
	int i;
	for(i=1;i<=row;i++){
		Vector v1=m[i];
		DotProduct(v1,v,&res[i]);
	}
	return TRUE;
}
Boolean rMatInterp(Vector v,Matrix m,Vector res){
	int row=NumRows(m);
	int col=NumCols(m);
	int vs=VectorSize(v);
	int rs=VectorSize(res);
	if(row!=vs||col!=rs){
		HError(-1,"MatInterp: res doesn't match with v*m");
		return FALSE;
	}
	int i,j;
	for(i=1;i<=col;i++){
		res[i]=0;
		for(j=1;j<=row;j++){
			res[i]+=v[j]*m[j][i];
		}
	}
	return TRUE;
}
Boolean VecSquare(Vector v1, Matrix m, Vector v2, float *res){
	int v1s=VectorSize(v1);
	int v2s=VectorSize(v2);
	int row=NumRows(m);
	int col=NumCols(m);
	if(v1s!=row||v2s!=col){
		HError(-1,"VecSquare: res doesn't match with v1'*m*v2");
		return FALSE;
	}
	int i,j;
	Vector v3=CreateVector(&gstack,row);
	MatInterp(m,v2,v3);
	DotProduct(v1,v3,res);
	FreeVector(&gstack,v3);
	return TRUE;
}
void printDMatrix(DMatrix mat, char* name){
	int row=NumDRows(mat);
	int col=NumDCols(mat);
	int i,j;
	printf("%s=[\n",name);
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			printf("%f\t",mat[i][j]);
		}
		if(i<row)
			printf(";\n");
		else
			printf("];");
	}
	printf("\n");
}
void printMatrix(Matrix mat, char* name){
	int row=NumRows(mat);
	int col=NumCols(mat);
	int i,j;
	printf("%s=[\n",name);
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			printf("%f\t",mat[i][j]);
		}
		if(i<row)
			printf(";\n");
		else
			printf("];");
	}
	printf("\n");
}

void printTriMat(TriMat tri, char* name){
	int row=TriMatSize(tri);
	int col=row;
	int i,j;
	printf("%s=[",name);
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			if(j<=i)
				printf("%f ",tri[i][j]);
			else
				printf("%f ",tri[j][i]);
		}
		if(i<row)
			printf(";\n");
		else
			printf("];");
	}
	printf("\n");
}
void printVector(Vector vec, char* name){
	int row=VectorSize(vec);
	int i;
	printf("%s=[",name);
	for(i=1;i<=row;i++){
		printf("%+6.3f ",vec[i]);
	}
	printf("];\n");
}
void printDVector(DVector vec, char* name){
	int row=DVectorSize(vec);
	int i;
	printf("%s=[",name);
	for(i=1;i<=row;i++){
		printf("%+6.3f ",vec[i]);
	}
	printf("];\n");
}
/* DLUDecompose: perform LU decomposition on Matrix a, the permutation
       of the rows is returned in perm and sign is returned as +/-1
       depending on whether there was an even/odd number of row 
       interchanges */
static Boolean DLUDecompose(MemHeap* x, DMatrix a, int *perm, int *sign)
{
   int i,imax,j,k,n;
   double scale,sum,xx,yy;
   DVector vv,tmp;
   
   n = NumDRows(a); imax = 0;
   vv = CreateDVector(x,n);
   *sign = 1;
   for (i=1; i<=n; i++) {
      scale = 0.0;
      for (j=1; j<=n; j++)
         if ((xx = fabs(a[i][j])) > scale )
            scale = xx;
      if (scale == 0.0) {
         /*HError(-1,"LUDecompose: Matrix is Singular");*/
         return(FALSE);
      }
      vv[i] = 1.0/scale;
   }
   for (j=1; j<=n; j++) {
      for (i=1; i<j; i++) {
         sum = a[i][j];
         for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
      }
      scale=0.0;
      for (i=j; i<=n; i++) {
         sum = a[i][j];
         for (k=1; k<j; k++) sum -= a[i][k]*a[k][j];
         a[i][j]=sum;
         if ( (yy=vv[i]*fabs(sum)) >=scale) {
            scale = yy; imax = i;
         }
      }
      if (j != imax) {
         tmp = a[imax]; a[imax] = a[j]; a[j] = tmp;
         *sign = -(*sign);
         vv[imax]=vv[j];
      }
      perm[j]=imax;
      if (a[j][j] == 0.0) {
         HError(-1,"LUDecompose: Matrix is Singular");
         return(FALSE);
      }
      if (j != n) {
         yy = 1.0/a[j][j];
         for (i=j+1; i<=n;i++) a[i][j] *= yy;
      }
   }
   FreeDVector(x,vv);
   return(TRUE);
}


/* EXPORT->DMatDet: determinant of a double matrix */
double DMatDet(DMatrix c)
{
   DMatrix a;
   double det;
   int n,perm[1600],i,sign;
   
   n=NumDRows(c);
   a=CreateDMatrix(&gstack,n,n);
   CopyDMatrix(c,a);                /* Make a copy of c */
   DLUDecompose(&gstack, a,perm,&sign);      /* Do LU Decomposition */
   det = sign;                     /* Calc Det(c) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   FreeDMatrix(&gstack,a);
   return det;
}



/* LinSolve: solve the set of linear equations Ax = b, returning
        the result x in  b */
void LinSolve(Matrix a, int *perm, float *b)
{
   int i,ii=0,ip,j,n;
   double sum;
   
   n=NumRows(a);
   for (i=1;i<=n;i++) {
      ip=perm[i]; sum=b[ip]; b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum -=a[i][j]*b[j];
      else
         if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n; i>=1; i--) {
      sum=b[i];
      for (j=i+1; j<=n; j++)
         sum -=a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}        

Boolean IsSingularMat(Matrix c){
	Matrix a;
	Boolean isSing;
	int sign;
	IntVec perm;
	Boolean tri=IsTriMat(c);
	n=tri?TriMatSize(c):NumRows(c);
	a=CreateMatrix(&gstack,n,n);
	perm=CreateIntVec(&gstack,n);
	if(tri){
		Tri2Mat(c,a);
	}
	else{
		CopyMatrix(c,a);           /* Make a copy of c */
	}
	if(!LUDecompose(a,perm,&sign)){      /* Do LU Decomposition */
		isSing=TRUE;
	}
	else{
		isSing=FALSE;
	}
	FreeMatrix(&gstack,a);
	return isSing;
}


/* EXPORT-> MatInvert: puts inverse of c in invc, returns Det(c) */
float MatInvert(Matrix c, Matrix invc)
{
   Matrix a;
   Vector col;
   float det;
   int sign;
   int n,i,j;
   IntVec perm;
   Boolean tri=IsTriMat(c);
   n=tri?TriMatSize(c):NumRows(c);
   a=CreateMatrix(&gstack,n,n);
   col=CreateVector(&gstack,n);
   perm=CreateIntVec(&gstack,n);
   if(tri){
	   Tri2Mat(c,a);
   }
   else{
	   CopyMatrix(c,a);           /* Make a copy of c */
   }
   if(!LUDecompose(a,perm,&sign)){      /* Do LU Decomposition */
	   if(tri)printTriMat(c,"tri");
	   else printMatrix(c,"mat");
	   HError(+1,"MatInvert: Matrix is Singular");
	   return 0;
   }
   for (j=1; j<=n; j++) {     /* Invert matrix */
      for (i=1; i<=n; i++)
         col[i]=0.0;
      col[j]=1.0;
      LinSolve(a,perm,col);
      for (i=1; i<=n; i++)
         invc[i][j] = col[i];
   }
   det = sign;                /* Calc log(det(c)) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   FreeMatrix(&gstack,a);
   return det;
}
 
/* DLinSolve: solve the set of linear equations Ax = b, returning
        the result x in  b */
static void DLinSolve(DMatrix a, int *perm, double *b)
{
   int i,ii=0,ip,j,n;
   double sum;
   
   n=NumDRows(a);
   for (i=1;i<=n;i++) {
      ip=perm[i]; sum=b[ip]; b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum -=a[i][j]*b[j];
      else
         if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n; i>=1; i--) {
      sum=b[i];
      for (j=i+1; j<=n; j++)
         sum -=a[i][j]*b[j];
      b[i]=sum/a[i][i];
   }
}       

/* Inverting a double matrix */
double DMatInvert(MemHeap *x, DMatrix c, DMatrix invc)
{
   DMatrix a;
   DVector col;
   double det;
   int sign;
   int n,i,j;
   IntVec perm;
   
   n=NumDRows(c);
   a=CreateDMatrix(x,n,n);
   col=CreateDVector(x,n);
   perm=CreateIntVec(x,n);
   CopyDMatrix(c,a);           /* Make a copy of c */
   if(!DLUDecompose(x, a,perm,&sign)){      /* Do LU Decomposition */
	   HError(1,"DMatInvert: Matrix is Singular");
	   /*printDMatrix(c,"X");*/
	   return 0;
   }
   for (j=1; j<=n; j++) {     /* Invert matrix */
      for (i=1; i<=n; i++)
         col[i]=0.0;
      col[j]=1.0;
      DLinSolve(a,perm,col);
      for (i=1; i<=n; i++)
         invc[i][j] = col[i];
   }  
   det = sign;                /* Calc log(det(c)) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   FreeDMatrix(x,a);
   return det;
}

/* EXPORT-> DMatCofact: generates the cofactors of row r of matrix c */
double DMatCofact(DMatrix c, int r, DVector cofact)
{
   DMatrix a;
   /*old code*/
   /*double col[100];*/
   double col[400];
   double det;
   int sign;
   int n,i;
   /*old code*/
   /*int perm[100];*/
   int perm[400];
   
   n=NumDRows(c);
   a=CreateDMatrix(&gstack,n,n);
   CopyDMatrix(c,a);                      /* Make a copy of c */
   if (! DLUDecompose(&gstack, a,perm,&sign))      /* Do LU Decomposition */
     return 0;
   det = sign;                         /* Calc det(c) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   for (i=1; i<=n; i++)
     col[i]=0.0;
   col[r]=1.0;
   DLinSolve(a,perm,col);
   for (i=1; i<=n; i++)
     cofact[i] = col[i]*det;
   FreeDMatrix(&gstack,a);
   return det;
}

/* EXPORT-> MatCofact: generates the cofactors of row r of matrix c */
double MatCofact(Matrix c, int r, Vector cofact)
{
   DMatrix a;
   DMatrix b;
   double col[100];
   float det;
   int sign;
   int n,i,perm[100];
 
   n=NumRows(c);
   a=CreateDMatrix(&gstack,n,n);
   b=CreateDMatrix(&gstack,n,n);
   Mat2DMat(c,b);
   CopyDMatrix(b,a);                      /* Make a copy of c */
   if (! DLUDecompose(&gstack, a,perm,&sign))      /* Do LU Decomposition */
     return 0;
   det = sign;                         /* Calc det(c) */
   for (i=1; i<=n; i++) {
      det *= a[i][i];
   }
   for (i=1; i<=n; i++)
     col[i]=0.0;
   col[r]=1.0;
   DLinSolve(a,perm,col);
   for (i=1; i<=n; i++)
     cofact[i] = col[i]*det;
   
   FreeDMatrix(&gstack,b);
   FreeDMatrix(&gstack,a);
   return det;
}

/* -------------------- Log Arithmetic ---------------------- */

/*
  The types LogFloat and LogDouble are used for representing
  real numbers on a log scale.  LZERO is used for log(0) 
  in log arithmetic, any log real value <= LSMALL is 
  considered to be zero.
*/

static LogDouble minLogExp;

/* EXPORT->LAdd: Return sum x + y on log scale, 
                sum < LSMALL is floored to LZERO */
/* z=LAdd(x,y)<-->z=log(e^x+e^y) */
LogDouble LAdd(LogDouble x, LogDouble y)
{
   LogDouble temp,diff,z;
   
   if (x<y) {
      temp = x; x = y; y = temp;
   }
   diff = y-x;
   if (diff<minLogExp) 
      return  (x<LSMALL)?LZERO:x;
   else {
      z = exp(diff);
      return x+log(1.0+z);
   }
}

/* EXPORT->LSub: Return diff x - y on log scale, 
                 diff < LSMALL is floored to LZERO */
LogDouble LSub(LogDouble x, LogDouble y)
{
   LogDouble diff,z;
   
   if (x<y)    
      HError(5271,"LSub: result -ve");
   diff = y-x;
   if (diff<minLogExp) 
      return  (x<LSMALL)?LZERO:x;
   else {
      z = 1.0 - exp(diff);
      return (z<MINLARG) ? LZERO : x+log(z);
   }
}

/* EXPORT->L2F: Convert log(x) to double, result is
                floored to 0.0 if x < LSMALL */
double   L2F(LogDouble x)
{
   return (x<LSMALL) ? 0.0 : exp(x);
}

/* -------------------- Random Numbers ---------------------- */


#ifdef UNIX
/* Prototype for C Library functions drand48 and srand48 */
double drand48(void);
void srand48(long);
#define RANDF() drand48()
#define SRAND(x) srand48(x)
#else
/* if not unix use ANSI C defaults */
#define RANDF() ((float)rand()/RAND_MAX)
#define SRAND(x) srand(x)
#endif

/* EXPORT->RandInit: Initialise random number generators 
           if seed is -ve, then system clock is used */
void RandInit(int seed)
{
   if (seed<0) seed = (int)time(NULL)%257;
   SRAND(seed);
}

/* EXPORT->RandomValue:  */
float RandomValue(void)
{
   return RANDF();
}

/* EXPORT->GaussDeviate: random number with a N(mu,sigma) distribution */
float GaussDeviate(float mu, float sigma)
{
   double fac,r,v1,v2,x;
   static int gaussSaved = 0; /* GaussDeviate generates numbers in pairs */
   static float gaussSave;    /* 2nd of pair is remembered here */


   if (gaussSaved) {
      x = gaussSave; gaussSaved = 0;
   }
   else {
      do {
         v1 = 2.0*(float)RANDF() - 1.0;
         v2 = 2.0*(float)RANDF() - 1.0;
         r = v1*v1 + v2*v2;
      }
      while (r>=1.0);
      fac = sqrt(-2.0*log(r)/r);
      gaussSaved = 1;
      gaussSave = v1*fac;
      x = v2*fac;
   }
   return x*sigma+mu;
}
static float eps=2.2204e-16;
static float threshold=0.01;
/*Randomise a non-negative matrix*/
static void randNNMatrix(Matrix X){
	int n=NumRows(X);
	int m=NumCols(X);
	int i,j;
	RandInit(-1);
	for(i=1;i<=n;i++){
		for(j=1;j<=m;j++){
			X[i][j]=RandomValue();
		}
	}
}
/*V:n*m, W=n*k, H=k*m*/
void NNMF(Matrix V, Matrix W, Matrix H, int k, int maxIter, Boolean speakout){
	int n=NumRows(V);
	int m=NumCols(V);
	randNNMatrix(W);
	randNNMatrix(H);
	float old_dist,new_dist;
	int i,j,iter;
	Matrix WH=NULL;
	Matrix Wt=CreateMatrix(&gstack,k,n);
	Matrix Ht=CreateMatrix(&gstack,m,k);
	Matrix WtV=CreateMatrix(&gstack,k,m);
	Matrix WtW=CreateMatrix(&gstack,k,k);
	Matrix WtWH=CreateMatrix(&gstack,k,m);
	Matrix VHt=CreateMatrix(&gstack,n,k);
	Matrix HHt=CreateMatrix(&gstack,k,k);
	Matrix WHHt=CreateMatrix(&gstack,n,k);
	WH=CreateMatrix(&gstack,n,m);
	MatMult(W,H,WH);
	old_dist=EuclMatDist(V,WH);
	Transpose(H,Ht);
	for(iter=1;iter<=maxIter;iter++){
		Transpose(W,Wt);
		MatMult(Wt,V,WtV);
		MatMult(Wt,W,WtW);
		MatMult(WtW,H,WtWH);
		for(i=1;i<=k;i++){
			for(j=1;j<=m;j++){
				H[i][j]*=WtV[i][j]/(WtWH[i][j]+eps);
			}
		}
		Transpose(H,Ht);
		MatMult(V,Ht,VHt);
		MatMult(H,Ht,HHt);
		MatMult(W,HHt,WHHt);
		for(i=1;i<=n;i++){
			for(j=1;j<=k;j++){
				W[i][j]*=VHt[i][j]/(WHHt[i][j]+eps);
			}
		}
		MatMult(W,H,WH);
		new_dist=EuclMatDist(V,WH);
		/*do something*/
		if(old_dist-new_dist<threshold){
			break;
		}
		if(speakout)
			printf("%f to %f with dif=%f @iter=%d\n",old_dist, new_dist, old_dist-new_dist,iter);
		old_dist=new_dist;
	}
	FreeMatrix(&gstack,Wt);
}
/*V:n*m, W=n*k, H=k*m*/
void NNMFpos(Matrix V, Matrix H, int k, int maxIter, Boolean speakout){
	int n=NumRows(V);
	int m=NumCols(V);
	int i,j;
	float colAvrSum=0,colSum;
	Matrix W=CreateMatrix(&gstack,n,k);
	Vector D=CreateVector(&gstack,k);
	NNMF(V,W,H,k,maxIter,speakout);
	ZeroVector(D);
	for(i=1;i<=n;i++){
		for(j=1;j<=k;j++){
			D[j]+=W[i][j];
		}
	}
	for(j=1;j<=m;j++){
		colSum=0;
		for(i=1;i<=k;i++){
			H[i][j]*=D[i];
			colSum+=H[i][j];
		}
		colAvrSum+=colSum;
		if(speakout){
			printf("chksum=%f @column=%d\n",colSum,j);
		}
	}
	if(speakout){
		printf("chksum average=%f \n",colAvrSum/m);
	}
}
/*-----------------------Signal processing related----------------- */
/*% DCT DMatrix*/
/*dctDMatrix = sqrt(2.0/numChans) * cos ( (PI/(2*numChans)) * (0:1:numceps)' * (1:2:(2*numChans-1)) );*/
/*numChans: num of filter bank, numCeps: num of MFCC coefficients, always including 0th coefficient*/
DMatrix DCT(int numChans, int numCeps){
	int i,j;
	int shifted_i;
	double mfnorm=sqrt(2.0/numChans);
	DMatrix dctMatrix=CreateDMatrix(&gstack,numCeps,numChans);
	for(i=0;i<numCeps;i++){
		shifted_i=i+(i==0?numCeps:0);/*since typical MFCC feature placed C0 at location 13 instead of 1, [13, 1..12]*/
		for(j=1;j<=numChans;j++){
			dctMatrix[shifted_i][j]=mfnorm*cos(PI*i*(j-0.5)/numChans);/*Starting from 0th coefficient*/
		}
	}
	return dctMatrix;
}
DMatrix liftDCT(MemHeap* x, int numChans, int numCeps, int numLifter){
	int i,j;
	int shifted_i;
	double mfnorm=sqrt(2.0/numChans);
	double lifter_i;
	DMatrix dctMatrix=CreateDMatrix(x,numCeps,numChans);
	for(i=0;i<numCeps;i++){/*[0..12]*/
		lifter_i=1+numLifter*0.5*sin(PI*i/numLifter);
		shifted_i=i+(i==0?numCeps:0);/*since typical MFCC feature placed C0 at location 13 instead of 1, [13, 1..12]*/
		for(j=1;j<=numChans;j++){
			dctMatrix[shifted_i][j]=mfnorm*cos(PI*i*(j-0.5)/numChans)*lifter_i;/*Starting from 0th coefficient*/
		}
	}
	return dctMatrix;
}
/*% inv DCT DMatrix
  invDctDMatrix = inv(sqrt(2.0/numChans) * cos ( (pi/(2*numChans)) * (0:1:(numChans-1))' * (1:2:(2*numChans-1)) ));
  invDctDMatrix = invDctDMatrix(:, 1:(numceps+1)); % cut off unnecessary cols
 */
/*DMatrix iDCT(int numChans, int numCeps){
	int i,j;
	DMatrix dct=DCT(numChans, numChans);/*Get a square DCT first*
	DMatrix iDct_square=CreateDMatrix(&gstack, numChans, numChans);
	DMatrix iDct_trunct=CreateDMatrix(&gstack, numChans, numCeps);
	DMatInvert(dct, iDct_square);
	for(i=1;i<=numChans;i++){
		for(j=1;j<=numCeps;j++){
			iDct_trunct[i][j]=iDct_square[i][j];
		}
	}
	return iDct_trunct;
}*/
/*DMatrix iLiftDCTx(int numChans, int numCeps, int numLifter){
	int i,j;
	DMatrix dct=liftDCT(numChans, numChans, numLifter);/*Get a square DCT first*
	DMatrix iDct_square=CreateDMatrix(&gstack, numChans, numChans);
	DMatrix iDct_trunct=CreateDMatrix(&gstack, numChans, numCeps);
	DMatInvert(dct, iDct_square);
	for(i=1;i<=numChans;i++){
		for(j=1;j<=numCeps;j++){
			iDct_trunct[i][j]=iDct_square[i][j];
		}
	}
	return iDct_trunct;
}*/
DMatrix iLiftDCT(MemHeap* x, int numChans, int numCeps, int numLifter){/*right pseudo-inverse matrix, Note that left pseudo-inverse does not exist since C:m by n with m<n*/
	int i,j;
	DMatrix Ct=CreateDMatrix(x, numChans, numCeps);
	DMatrix CCt=CreateDMatrix(x, numCeps, numCeps);
	DMatrix iCCt=CreateDMatrix(x, numCeps, numCeps);
	DMatrix C=liftDCT(x, numChans, numCeps, numLifter);
	DMatrix iC=CreateDMatrix(x, numChans, numCeps);
	DTranspose(C,Ct);
	DMatMult(C,Ct,CCt);
	DMatInvert(x, CCt, iCCt);
	DMatMult(Ct,iCCt,iC);
	return iC;
}


/*% cepstral liftering coefficient
	cepDMatrix = diag(1+cepLifter*0.5*sin((0:numceps)'*pi/cepLifter));
*/
/*void cepLift(int numCeps, int numLifter, DVector mu, DMatrix cov){
	int i,j;
	int row=NumDRows(cov);
	int col=NumDCols(cov);
	int cepIdx_i,cepIdx_j;
	double lifter_i,lifter_j;
	for(i=1;i<=row;i++){
		cepIdx_i=(i-1)%numCeps;/*[0..12]*
		lifter_i=1+numLifter*0.5*sin(PI*cepIdx_i/numLifter);
		mu[i]=mu[i]*lifter_i;
		for(j=1;j<=col;j++){
			cepIdx_j=(j-1)%numCeps;/*[0..12]*
			lifter_j=1+numLifter*0.5*sin(PI*cepIdx_j/numLifter);
			cov[i][j]=lifter_i*cov[i][j]*lifter_j;
		}
	}
}*/
/*% inv cepstral liftering coefficient
  invCepDMatrix = diag(1./(1+cepLifter*0.5*sin((0:numceps)*pi/cepLifter)));
*/
/*void icepLift(int numCeps, int numLifter, DVector mu, DMatrix cov){
	int i,j;
	int row=NumDRows(cov);
	int col=NumDCols(cov);
	int cepIdx_i,cepIdx_j;
	double ilifter_i,ilifter_j;
	for(i=1;i<=row;i++){
		cepIdx_i=(i-1)%numCeps;/*[0..12]*
		ilifter_i=1/(1+numLifter*0.5*sin(PI*cepIdx_i/numLifter));
		mu[i]=mu[i]*ilifter_i;
		for(j=1;j<=col;j++){
			cepIdx_j=(j-1)%numCeps;/*[0..12]*
			ilifter_j=1/(1+numLifter*0.5*sin(PI*cepIdx_j/numLifter));
			cov[i][j]=ilifter_i*cov[i][j]*ilifter_j;
		}
	}
}*/

/*numChans: num of filter bank, numCeps: num of MFCC coefficients, always including 0th coefficient*/
/*dctDMatrix = kron(eye(numDctDVectors), dctDMatrix);/* % for trajectory model if numDctDVectors > 1*/
/*DMatrix trajDCT(int trajLen, int numChans, int numCeps){
	DMatrix dct=DCT(numChans, numCeps);
	DMatrix trajDct=CreateDMatrix(&gstack,trajLen*numCeps,trajLen*numChans);
	duplicatem(dct,numCeps,numChans,trajDct);
	return trajDct;
}
/*numChans: num of filter bank, numCeps: num of MFCC coefficients, always including 0th coefficient*/
/*invDctDMatrix = kron(eye(numDctDVectors), invDctDMatrix); % for trajectory model if numDctDVectors > 1*/
/*DMatrix itrajDCT(int trajLen, int numChans, int numCeps){
	DMatrix iDct=iDCT(numChans,numCeps);
	DMatrix trajiDct=CreateDMatrix(&gstack,trajLen*numChans,trajLen*numCeps);
	duplicatem(iDct,numChans, numCeps,trajiDct);
	return trajiDct;
}*/


DMatrix lLiftDCT(MemHeap* x, int numFrms, int numChans, int numCeps, int numLifter){
	DMatrix dct=liftDCT(x, numChans, numCeps, numLifter);
	DMatrix lDct=CreateDMatrix(x,numFrms*numCeps,numFrms*numChans);
	DupDMatrix(dct,numCeps,numChans,lDct);
	return lDct;
}

DMatrix ilLiftDCT(MemHeap* x, int numFrms, int numChans, int numCeps, int numLifter){
	DMatrix iDct=iLiftDCT(x, numChans,numCeps,numLifter);
	DMatrix trajiDct=CreateDMatrix(x,numFrms*numChans,numFrms*numCeps);
	DupDMatrix(iDct,numChans, numCeps,trajiDct);
	return trajiDct;
}
/**
 * Covert the cepstral domain statistics to linear spectral domain
 * input: c_mu, cepstral mean; c_cov, cepstral covariance matrix; iDct, inverse DCT transform
 * output: l_mu, linear mean; l_cov, linear covariance matrix
 */
void cep2linear(DVector c_mu, DMatrix iDct, DMatrix c_cov, DVector l_mu, DMatrix l_cov){
	int i,j;
	int numChans=NumDRows(iDct);
	int numCeps=NumDCols(iDct);
	double temp;
	DMatrix iCSig=CreateDMatrix(&gstack,numChans,numCeps);
	DMatrix iCSigCt=CreateDMatrix(&gstack,numChans,numChans);
	DMatrix iDctT=CreateDMatrix(&gstack,numCeps,numChans);
	DMatMult(iDct,c_cov,iCSig);
	DTranspose(iDct,iDctT);
	DMatMult(iCSig,iDctT,iCSigCt);
	for(i=1;i<=numChans;i++){
		DDotProduct(iDct[i],c_mu,&temp);
		l_mu[i]=exp(temp+0.5*iCSigCt[i][i]);
	}
	for(i=1;i<=numChans;i++){
		for(j=1;j<=numChans;j++){
			l_cov[i][j]=l_mu[i]*(exp(iCSigCt[i][j])-1)*l_mu[j];
		}
	}
	FreeDMatrix(&gstack,iCSig);
}
void VTS(MemHeap* x, DMatrix dct, DMatrix idct, DVector exPoint, DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DVector g0){/*right-pseudo inverse dct*/
	double dot,pexpdot;
	DVector log1pexp, f;
	int i,j,k;
	int numCeps=NumDRows(dct);
	int numChans=NumDCols(dct);
	f=CreateDVector(x,numChans);
	log1pexp=CreateDVector(x,numChans);
	for(i=1;i<=numChans;i++){
		DDotProduct(idct[i],exPoint,&dot);
		pexpdot=1+exp(dot>=50?50:dot);/*just in case of large value*/
		f[i]=1.0/pexpdot;/*F(ex)=1/(1+exp(C^(-1)*ex))*/
		log1pexp[i]=log(pexpdot);/*log(1+exp(C^(-1)*ex))*/
	}
	for(i=1;i<=numCeps;i++){
		for(j=1;j<=numCeps;j++){
			G[i][j]=0;
			for(k=1;k<=numChans;k++){
				G[i][j]+=f[k]*dct[i][k]*idct[k][j];/*G=C*diag(.)*C^(-1)*/
			}
			Gt[j][i]=G[i][j];/*Note that G is not a symmetric matrix, this is important for consistent results*/
			F[i][j]=(i==j)?(1-G[i][j]):-G[i][j];/*F=I-G*/
			Ft[j][i]=F[i][j];
		}
	}
	if(g0){
		for(i=1;i<=numCeps;i++){
			DDotProduct(dct[i],log1pexp,&g0[i]);/*g(z)=C*log(1+exp(C^(-1)*mu))*/
		}
	}
	FreeDVector(x,f);
}



void VtsUpVar(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DMatrix GCx, DMatrix GCxGt, DMatrix FCn, DMatrix FCnFt, DMatrix x_cep_cv, DMatrix n_cep_cv, DMatrix y_cep_cv){
	int i,j;
	int vSize=NumDRows(G);
	DMatMult(G,x_cep_cv,GCx);
	DMatMult(GCx,Gt,GCxGt);/*A*Cov_x*A^T*/
	DMatMult(F,n_cep_cv,FCn);
	DMatMult(FCn,Ft,FCnFt);/*(I-A)*Cov_n*(I-A)^T*/
	for(i=1;i<=vSize;i++){
		for(j=1;j<=vSize;j++){
			y_cep_cv[i][j]=GCxGt[i][j]+FCnFt[i][j];/*Cov_y=A*Cov_x*A^T+A*Cov_h*A^T+(I-A)*Cov_n*(I-A)^T*/
		}
	}
}

/**
 * A multiply diagonal matrix
 */
Boolean DMatDiag(DMatrix A, DVector d, DMatrix Ad){
	int i,j;
	int row=NumDRows(A);
	int col=NumDCols(A);
	int vSize=DVectorSize(d);
	if(vSize!=col){
		HError(1,"Matrix size A|B doesn't match diagonal matrix size, A(%d %d) d(%d)",row,col,vSize);
		return FALSE;
	}
	for(i=1;i<=row;i++){
		for(j=1;j<=col;j++){
			Ad[i][j]=A[i][j]*d[j];
		}
	}
	return TRUE;
}
/**
 * A multiply diagonal matrix and another full matrix
 */
Boolean DMatDiagDMat(DMatrix A, DVector d, DMatrix B, DMatrix AdB){
	int i,j,k;
	int row=NumDRows(A);
	int col=NumDCols(A);
	int vSize=VectorSize(d);
	int row2=NumDRows(B);
	int col2=NumDCols(B);
	if(vSize!=col||vSize!=row2){
		HError(1,"Matrix size A|B doesn't match diagonal matrix size, A(%d %d) d(%d) B(%d %d)",row,col,vSize,row2,col2);
		return FALSE;
	}
	for(i=1;i<=row;i++){
		for(j=1;j<=col2;j++){
			AdB[i][j]=0;
			for(k=1;k<=vSize;k++){
				AdB[i][j]+=A[i][k]*d[k]*B[k][j];
			}
		}
	}
	return TRUE;
}
/*diagonal covariance case*/
void dVtsUpVar(DMatrix G, DMatrix Gt, DMatrix F, DMatrix Ft, DVector x_cep_cv, DVector n_cep_cv, DVector y_cep_cv){
	int i,k;
	int vSize=NumDRows(G);
	double GCxGt_ii,FCnFt_ii;
	for(i=1;i<=vSize;i++){
		GCxGt_ii=0;
		FCnFt_ii=0;
		for(k=1;k<=vSize;k++){
			GCxGt_ii+=G[i][k]*x_cep_cv[k]*Gt[k][i];/*A*Cov_x*A^T*/
			FCnFt_ii+=F[i][k]*n_cep_cv[k]*Ft[k][i];/*(I-A)*Cov_n*(I-A)^T*/
		}
		y_cep_cv[i]=GCxGt_ii+FCnFt_ii;
	}
}
/**
 * Covert the linear spectral domain statistics to cepstral domain
 * input: l_mu, linear mean; l_cov, linear covariance matrix; Dct, DCT transform
 * output: c_mu, cepstral mean; c_cov, cepstral covariance matrix
 */
void linear2cep(DVector l_mu, DMatrix dct, DMatrix l_cov, DVector c_mu, DMatrix c_cov){
	int i,j;
	int numChans=NumDCols(dct);
	int numCeps=NumDRows(dct);
	DVector tempV=CreateDVector(&gstack,numChans);
	DMatrix V=CreateDMatrix(&gstack,numChans,numChans);
	DMatrix ClogV=CreateDMatrix(&gstack,numCeps,numChans);
	DMatrix dctT=CreateDMatrix(&gstack,numChans,numCeps);
	for(i=1;i<=numChans;i++){
		for(j=1;j<=numChans;j++){
			V[i][j]=(1/l_mu[i])*l_cov[i][j]*(1/l_mu[j]);
		}
	}
	for(i=1;i<=numCeps;i++){
		for(j=1;j<=numChans;j++){
			tempV[j]=log(l_mu[j])-0.5*log(V[j][j]+1);
		}
		DDotProduct(dct[i],tempV,&c_mu[i]);
	}
	for(i=1;i<=numChans;i++){
		for(j=1;j<=numChans;j++){
			V[i][j]=log(V[i][j]+1);
		}
	}
	/*printf("%d\n",IsPositiveDef(V));*/
	DMatMult(dct,V,ClogV);
	DTranspose(dct,dctT);
	DMatMult(ClogV,dctT,c_cov);
	FreeDVector(&gstack,tempV);
}

/**
 * Move C0 to location 1
 */
void moveC0topx(int numCeps, DVector mu, DMatrix cov){return;
	int i,j,newi,newj;
	int row, col;
	row=col=DVectorSize(mu);
	if(cov&&NumDRows(cov)!=row){
		HError(1,"Mean and Variance size inconsistent");
	}
	DVector tempV=CreateDVector(&gstack,row);
	DMatrix tempX=CreateDMatrix(&gstack,row,col);
	for(i=1;i<=row;i++){
		newi=(i%numCeps+1)+((i-1)/numCeps)*numCeps;/*[1:12, 0]+1=>[2:13, 1]+offset*/
		tempV[newi]=mu[i];
		if(cov){
			for(j=1;j<=col;j++){
				newj=(j%numCeps+1)+((j-1)/numCeps)*numCeps;
				tempX[newi][newj]=cov[i][j];
			}
		}
	}
	CopyDVector(tempV,mu);
	if(cov){
		CopyDMatrix(tempX,cov);
	}
	FreeDVector(&gstack,tempV);
}
/**
 * Move C0 to location 13
 */
void moveC0btmx(int numCeps, DVector mu, DMatrix cov){return;
	int i,j,newi,newj;
	int row, col;
	row=col=DVectorSize(mu);
	if(cov&&NumDRows(cov)!=row){
		HError(1,"Mean and Variance size inconsistent");
	}
	DVector tempV=CreateDVector(&gstack,row);
	DMatrix tempX=CreateDMatrix(&gstack,row,col);
	for(i=1;i<=row;i++){
		newi=(i-1)%numCeps+((i-1)/numCeps)*numCeps+((i%numCeps==1)?numCeps:0);/*[0:12]+offset+first?13:0=>[13, 1:12]+offset*/
		tempV[newi]=mu[i];
		if(cov){
			for(j=1;j<=col;j++){
				newj=(j-1)%numCeps+((j-1)/numCeps)*numCeps+((j%numCeps==1)?numCeps:0);
				tempX[newi][newj]=cov[i][j];
			}
		}
	}
	CopyDVector(tempV,mu);
	if(cov){
		CopyDMatrix(tempX,cov);
	}
	FreeDVector(&gstack,tempV);
}
void Vec2iDVec(Vector v, DVector dv){
	int i;
	int vsize=VectorSize(v);
	for(i=1;i<=vsize;i++){
		dv[i]=1/v[i];
	}
}
void Vec2DVec(Vector v, DVector dv){
	int i;
	int vsize=VectorSize(v);
	for(i=1;i<=vsize;i++){
		dv[i]=v[i];
	}
}
void DVec2Vec(DVector dv, Vector v){
	int i;
	int vsize=DVectorSize(dv);
	for(i=1;i<=vsize;i++){
		v[i]=dv[i];
	}
}
void truncDVec2Vec(DVector dv, Vector v){
	int i;
	int v1size=DVectorSize(dv);
	int v2size=VectorSize(v);
	int trim=(v1size-v2size)/2;
	for(i=trim+1;i<=v1size-trim;i++){
		v[i-trim]=dv[i];
	}
}
void truncDVec2DVec(DVector dv1, DVector dv2){
	int i;
	int v1size=DVectorSize(dv1);
	int v2size=DVectorSize(dv2);
	int trim=(v1size-v2size)/2;
	for(i=trim+1;i<=v1size-trim;i++){
		dv2[i-trim]=dv1[i];
	}
}
void truncDMat2DMat(DMatrix dm1, DMatrix dm2){
	int i,j;
	int m1r=NumDRows(dm1);
	int m1c=NumDCols(dm1);
	int m2r=NumDRows(dm2);
	int m2c=NumDCols(dm2);
	int trimr=(m1r-m2r)/2;
	int trimc=(m1c-m2c)/2;
	for(i=trimr+1;i<=m1r-trimr;i++){
		for(j=trimc+1;j<=m1c-trimc;j++){
			dm2[i-trimr][j-trimc]=dm1[i][j];
		}
	}
}
void truncDMat2Mat(DMatrix dm1, Matrix dm2){
	int i,j;
	int m1r=NumDRows(dm1);
	int m1c=NumDCols(dm1);
	int m2r=NumRows(dm2);
	int m2c=NumCols(dm2);
	int trimr=(m1r-m2r)/2;
	int trimc=(m1c-m2c)/2;
	for(i=trimr+1;i<=m1r-trimr;i++){
		for(j=trimc+1;j<=m1c-trimc;j++){
			dm2[i-trimr][j-trimc]=dm1[i][j];
		}
	}
}
void dvec2diag(DVector v, DMatrix m){
	int i,j;
	int vsize=DVectorSize(v);
	for(i=1;i<=vsize;i++){
		for(j=1;j<=vsize;j++){
			m[i][j]=(i==j)?v[i]:0;
		}
	}
}
void dvec2idiag(DVector v, DMatrix m){
	int i,j;
	int vsize=DVectorSize(v);
	for(i=1;i<=vsize;i++){
		for(j=1;j<=vsize;j++){
			if(v[i]==0)
				HError(1, "dvec2idiag: Vector element is not invertible");
			m[i][j]=(i==j)?1.0/v[i]:0;
		}
	}
}
void MatDiag2iVec(Matrix m, Vector v){
	int i;
	int vsize=VectorSize(v);
	for(i=1;i<=vsize;i++){
		if(m[i][i]==0)
			HError(1, "diag2ivec:Matrix diagonal element is not invertible");
		v[i]=1/m[i][i];
	}
}
void MatDiag2Vec(Matrix m, Vector v){
	int i;
	int vsize=VectorSize(v);
	for(i=1;i<=vsize;i++){
		v[i]=m[i][i];
	}
}
void DMatDiag2Vec(DMatrix m, Vector v){
	int i;
	int vsize=VectorSize(v);
	for(i=1;i<=vsize;i++){
		v[i]=m[i][i];
	}
}
void DupDVector(DVector v, int numElem, DVector longv){
	int i,j;
	int vsize0=DVectorSize(v);
	int vsizel=DVectorSize(longv);
	if(numElem>0)
		vsize0=numElem;
	int dupsize=vsizel/vsize0;
	for(i=1;i<=dupsize;i++){
		for(j=1;j<=vsize0;j++){
			longv[(i-1)*vsize0+j]=v[j];
		}
	}
}
void DupDMatrix(DMatrix m, int numRow, int numCol, DMatrix longm){
	int i,j,k;
	int offsetr,offsetc;
	int row0=NumDRows(m);
	int col0=NumDCols(m);
	int rowl=NumDRows(longm);
	int coll=NumDCols(longm);
	if(numRow>0)
		row0=numRow;
	if(numCol>0)
		col0=numCol;
	int dupsize=rowl/row0;
	ZeroDMatrix(longm);
	for(i=1;i<=dupsize;i++){
		offsetr=(i-1)*row0;
		offsetc=(i-1)*col0;
		for(j=1;j<=row0;j++){
			for(k=1;k<=col0;k++){
				longm[offsetr+j][offsetc+k]=m[j][k];
			}
		}
	}
}
/* Kronecker tensor product*/
void kronecker(DMatrix a, DMatrix b, DMatrix res)
{
	int i,j,ii,jj,nc,nr,nra,nca,nrb,ncb;
	int k,q;
	nra=NumDRows(a); nrb=NumDRows(b);nca= NumDCols(a); ncb =NumDCols(b);
	nr=nra*nrb; nc=nca*ncb;
	for (i=1;i<=nra;i++){
		for (j=1;j<=nca;j++){
			for (ii=1;ii<=nrb;ii++){
				for (jj=1;jj<=ncb;jj++){
					k = nrb*(i-1)+ii;
					q = ncb*(j-1)+jj;
					res[k][q]=a[i][j]*b[ii][jj];
				}
			}
		}
	}
}
/*------------------GMM estimation related----------------------*/
static float floorV[500];
Point CreatePoint(MemHeap* x, Vector value, Ptr owner, Ptr hook){
	Point pnt=(Point)New(x,sizeof(struct _Point));
	pnt->owner=owner;
	pnt->value=value;
	pnt->hook=hook;
	pnt->next=NULL;
	return pnt;
}
void FreePoint(MemHeap* x, Point pnt){
	Dispose(x,pnt);
}
Group CreateGroup(MemHeap* x, int vsize, Boolean fc){
	int i;
	Group grp=(Group)New(x,sizeof(struct _PntGrp));
	Point plist=CreatePoint(x,NULL,NULL,NULL);
	grp->plist=plist;
	grp->mean=CreateVector(x,vsize);
	grp->acc1=CreateVector(x,vsize);
	ZeroVector(grp->acc1);
	grp->fc=fc;
	if(fc){
		grp->inv=CreateTriMat(x,vsize);
		grp->acc2x=CreateTriMat(x,vsize);
		ZeroMatrix(grp->acc2x);
	}
	else{
		grp->var=CreateVector(x,vsize);
		grp->acc2=CreateVector(x,vsize);
		ZeroVector(grp->acc2);
	}
	for(i=1;i<=vsize;i++){
		grp->mean[i]=0;
		if(fc){
			grp->inv[i][i]=1;
		}
		else{
			grp->var[i]=1;
		}
	}
	grp->gConst=-0.5*(vsize*log(2*PI)+log(1.0));
	grp->occ=0;
	grp->floorV=floorV;
	grp->hardOcc=0;
	return grp;
}
void AddIntoGroup(Group grp, Point pnt, float post, Boolean insert, Boolean callOcc){
	int i,j;
	float x;
	int vsize=VectorSize(pnt->value);
	if(insert){/*In the case of GMM, only closest cluster is chosen to be inserted*/
		Point temp=grp->plist->next;
		grp->plist->next=pnt;
		pnt->next=temp;
		grp->hardOcc++;
	}
	if(callOcc){
		Vector hjm=(Vector) (pnt->hook);
		for(i=1;i<=vsize;i++){
			grp->acc1[i]+=hjm[i];
			grp->occ+=hjm[i];
		}
		return;
	}
	for(i=1;i<=vsize;i++){
		x=pnt->value[i];
		grp->acc1[i]+=post*x;
		if(grp->fc){
			for(j=1;j<=i;j++){
				grp->acc2x[i][j]+=post*pnt->value[i]*pnt->value[j];
			}
		}
		else{
			grp->acc2[i]+=post*x*x;
		}
	}
	grp->occ+=post;
}
float Euclidean(Vector v1, Vector v2){
	int i;
	int vsize=VectorSize(v1);
	float sum=0;
	for(i=1;i<=vsize;i++){
		float x=v1[i]-v2[i];
		sum+=x*x;
	}
	return sqrt(sum);
}

LogFloat pdf(Group grp, Vector v1){
	int i;
	int vsize=VectorSize(v1);
	float sum=0;
	float res;
	for(i=1;i<=vsize;i++){
		float x=v1[i]-grp->mean[i];
		sum+=x*x/grp->var[i];
	}
	res=grp->gConst-0.5*sum+log(grp->weight);
	return res;
}
LogFloat fpdf(MemHeap* x, Group grp, Vector vec){
	float sum;
	int i,j;
	int vecSize=VectorSize(vec);
	Vector xmm;
	TriMat m = grp->inv;

	xmm = CreateVector(x,vecSize);
	for (i=1;i<=vecSize;i++)
		xmm[i] = vec[i] - grp->mean[i];
	sum = 0.0;
	for (j=1;j<vecSize;j++)
		for (i=j+1;i<=vecSize;i++)
			sum += xmm[i]*xmm[j]*m[i][j];
	sum *= 2;
	for (i=1;i<=vecSize;i++)
		sum += xmm[i] * xmm[i] * m[i][i];
	FreeVector(x,xmm);
	return grp->gConst-0.5*sum+log(grp->weight);
}
float SearchAndInsert(MemHeap* x, Group* grps, Point pnt, int numClust, DistKind distK, Boolean insert, Boolean callOcc){
	int i,j,vsize;
	double curMinDist=BIG_FLOAT, dist, dist2;
	Group closestGrp=NULL;
	if(distK==Gaussi){
		Vector poses=CreateVector(x,numClust);
		float sum=LZERO,like,bestlike=LZERO;
		int besti=0;
		for(i=0;i<numClust;i++){
			poses[i+1]=like=grps[i]->fc?fpdf(x,grps[i],pnt->value):pdf(grps[i],pnt->value);
			sum=LAdd(poses[i+1],sum);
			if(like>bestlike){
				bestlike=like;
				besti=i;
			}
		}
		for(i=0;i<numClust;i++){
			dist=-poses[i+1];/*distance is the inverse of likelihood*/
			poses[i+1]=exp(poses[i+1]-sum);
			AddIntoGroup(grps[i],pnt,poses[i+1],besti==i&&insert,FALSE);
		}
		FreeVector(x,poses);
		return sum;
	}
	else if(distK==Euclid){
		for(i=0;i<numClust;i++){
			dist=Euclidean(grps[i]->mean,pnt->value);
			if(dist<curMinDist){
				closestGrp=grps[i];
				curMinDist=dist;
			}
		}
		dist2=curMinDist*curMinDist;
		AddIntoGroup(closestGrp,pnt,1.0,insert,FALSE);
		return dist2;
	}
	else if(distK==Userdt){/*criteria based distance*/
		vsize=VectorSize(pnt->value);
		Vector hjm=(Vector) (pnt->hook);
		for(i=0;i<numClust;i++){/*search a cluster, offering the maximum objective value, given the incoming gamma_jmi*/
			dist=0;
			for(j=1;j<=vsize;j++){/*maximise sum_i gamma_jmi*log(wjm_i), or minimise negative summation*/
				dist-=hjm[j]*(grps[i]->mean[j]<=0?LZERO:log(grps[i]->mean[j]));
			}
			if(dist<curMinDist){
				closestGrp=grps[i];
				curMinDist=dist;
			}
		}
		dist2=curMinDist;
		AddIntoGroup(closestGrp,pnt,1.0,insert,callOcc);
		return dist2;
	}
	else{
		return 0;
	}
}

void UpdateGroup(MemHeap* x, Group grp, int numPnts, Boolean upVar){
	int i,j;
	int vsize=VectorSize(grp->mean);
	float ldet=0;
	TriMat var=CreateTriMat(x,vsize);
	float* floorV=grp->floorV;
	float floorv=1E-4;
	for(i=1;i<=vsize;i++){/*update to regular mean*/
		grp->mean[i]=grp->occ>=1?(grp->acc1[i]/grp->occ):0;/*If no sample, then reset to 0 mean*/
		if(upVar){
			if(grp->fc){
				for(j=1;j<=i;j++){
					var[i][j]=grp->occ>=1?((1/grp->occ)*grp->acc2x[i][j]-grp->mean[i]*grp->mean[j]):(i==j?1:0);/*If no sample, then reset to unit variance*/
				}
				if(var[i][i]<floorV[i])
					var[i][i]=floorV[i];/*In case of zero variance or too small variance*/
			}
			else{
				grp->var[i]=grp->occ>=1?(grp->acc2[i]/grp->occ-grp->mean[i]*grp->mean[i]):1;/*If no sample, then reset to unit variance*/
				if(grp->var[i]<(floorV[i]>0?floorV[i]:floorv))
					grp->var[i]=(floorV[i]>0?floorV[i]:floorv);/*In case of zero variance or too small variance*/
				ldet+=log(grp->var[i]);
			}
		}
	}
	if(upVar){
		if(grp->fc){
			ldet=CovInvert(&gstack, var,grp->inv);
		}
		grp->gConst=-0.5*(vsize*log(2*PI)+ldet);/*Fix gConst*/
	}
	grp->weight=grp->occ/numPnts;
	FreeTriMat(x,var);
}

void UpdateGroups(MemHeap *x, Group* grps, int numClust, int numPnts, Boolean upVar){
	int i;
	for(i=0;i<numClust;i++){
		UpdateGroup(x, grps[i], numPnts, upVar);
	}
}
void ResetGroups(Group* grps, int numClust){
	int i;
	for(i=0;i<numClust;i++){
		grps[i]->occ=0;
		grps[i]->hardOcc=0;
		grps[i]->plist->next=NULL;
		ZeroVector(grps[i]->acc1);
		if(grps[i]->fc)
			ZeroTriMat(grps[i]->acc2x);
		else
			ZeroVector(grps[i]->acc2);
	}
}

Group* InitKmeansMU(Point* pnts, int numClust, int vSize, int numPnts, Boolean fc, Boolean tiedCov){
	int i,j,grpCount=1,tempCount=1,k;
	double x;
	double measure;
	const float pertDepth = 0.2;
	float normWgt;
	Group *grps,*pGrps;
	TriMat fcvar;
	pGrps=(Group*)New(&gstack,sizeof(Group)*numClust);
	grps=(Group*)New(&gstack,sizeof(Group)*numClust);
	fcvar=CreateTriMat(&gstack,vSize);
	float *floorV;
	int toSplitCount;
	/*for(i=0;i<numPnts;i++){
		for(j=1;j<=vSize;j++){
			pnts[i]->value[j]=1000*pnts[i]->value[j];
		}
	}*/
	for(i=0;i<numClust;i++){
		grps[i]=CreateGroup(&gstack,vSize,fc);
		pGrps[i]=CreateGroup(&gstack,vSize,fc);
	}
	floorV=pGrps[0]->floorV;
	pGrps[0]->weight=1;
	for(i=0;i<numPnts;i++){
		SearchAndInsert(&gstack,pGrps, pnts[i], grpCount, Gaussi, FALSE, FALSE);
	}
	UpdateGroups(&gstack,pGrps, grpCount, numPnts, TRUE);
	if(fc){
		CovInvert(&gstack, pGrps[0]->inv,fcvar);/*Note that you can't just invert the diagonal element to get the variance*/
	}
	for(j=1;j<=vSize;j++){
		floorV[j]=fc?(0.01*fcvar[j][j]):(0.01*pGrps[0]->var[j]);
	}
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
		for(j=0;j<4;j++){/*do 4 iteration EM estimation*/
			measure=0;
			ResetGroups(grps,grpCount);/*reset the accumulators*/
			for(i=0;i<numPnts;i++){
				measure+=SearchAndInsert(&gstack,grps, pnts[i], grpCount, Gaussi,FALSE,FALSE);
			}
			UpdateGroups(&gstack,grps, grpCount, numPnts, tiedCov?FALSE:TRUE);/*update mean&variance&weight of the cluster*/
			printf("\t\tIteration %d with measure %f\n",j,measure/numPnts);fflush(stdout);
			/*printVector(grps[0]->mean,"m1");
			printVector(grps[1]->mean,"m1");*/
		}
		for(i=0;i<grpCount;i++){
			pGrps[i]->weight=grps[i]->weight;
			pGrps[i]->gConst=grps[i]->gConst;
			CopyVector(grps[i]->mean,pGrps[i]->mean);
			if(fc)
				CopyMatrix(grps[i]->inv,pGrps[i]->inv);
			else
				CopyVector(grps[i]->var,pGrps[i]->var);
		}
	}
	/*for(i=0;i<numPnts;i++){
		for(j=1;j<=vSize;j++){
			pnts[i]->value[j]=pnts[i]->value[j]/scale;
		}
	}
	for(i=0;i<numClust;i++){
		for(j=1;j<=vSize;j++){
			grps[i]->mean[j]=grps[i]->mean[j]/scale;
			if(fc){
				for(k=1;k<=j;k++){
					grps[i]->inv[j][k]=grps[i]->inv[j][k]*scale*scale;
				}
			}
			else
				grps[i]->var[j]=grps[i]->var[j]/(scale*scale);
		}
	}*/
	return grps;
}
/*min c'x, s.t. Ax=b, x>=0, assume the input x is feasible*/


void Interior_point_method(Vector c, Matrix A, Vector b, Vector x, char *ident){
	const int maxK=100;
	const float objThresh=1E-3;
	const float sigThresh=1E-3;
	const float conThresh=1E-3;
	const float alpThresh=1E-10;
	int m=NumRows(A);
	int n=NumCols(A);
	int M=2*n+m;
	int i,j,k;
	int sign;
	float mu,sigma,alpha,initObj,prevObj,curObj,temp;
	Matrix Jac=CreateMatrix(&gstack,M,M);
	Vector s=CreateVector(&gstack,n);
	Vector lmd=CreateVector(&gstack,m);
	Vector rb=CreateVector(&gstack,m);
	Vector residual=CreateVector(&gstack,M);
	IntVec perm=CreateIntVec(&gstack,M);
	Matrix tempX=CreateMatrix(&gstack,M,M);
	for(i=1;i<=n;i++){/*Initialization, >0*/
		s[i]=1;
	}
	for(i=1;i<=m;i++){/*Initialization*/
		lmd[i]=1;
	}
	for(i=1;i<=M;i++){/*Jacobian=[0 A' I; A 0 0; S^ 0 X^]*/
		for(j=1;j<=M;j++){
			Jac[i][j]=0;
			if(i<=n){/*first row*/
				if(j>n&&j<=n+m){/*second column*/
					Jac[i][j]=A[j-n][i];
				}
				if(j>n+m&&i==j-n-m){/*third column*/
					Jac[i][j]=1;
				}
			}
			if(i>n&&i<=n+m&&j<=n){/*second row, first column*/
				Jac[i][j]=A[i-n][j];
			}
		}
	}
	prevObj=0;
	for(i=1;i<=n;i++){
		prevObj+=c[i]*x[i];
	}
	initObj=prevObj;
	sigma=1;
	for(k=0;k<maxK;k++){/*Iterative update*/
		mu=0;
		for(i=1;i<=n;i++){/*update duality measure*/
			mu+=x[i]*s[i];
		}
		mu/=n;
		for(i=1;i<=n;i++){/*update Jacobian element S^ and X^*/
			Jac[i+n+m][i]=s[i];
		}
		for(i=n+m+1;i<=M;i++){
			Jac[i][i]=x[i-n-m];
		}
		for(i=1;i<=n;i++){/*rc=A'lmd+s-c*/
			temp=s[i]-c[i];
			for(j=1;j<=m;j++){
				temp+=A[j][i]*lmd[j];
			}
			residual[i]=-temp;
		}
		for(i=1;i<=m;i++){/*rb=Ax-b*/
			temp=-b[i];
			for(j=1;j<=n;j++){
				temp+=A[i][j]*x[j];
			}
			rb[i]=temp;
			residual[i+n]=-temp;
		}
		for(i=1;i<=n;i++){/*XSe*/
			residual[i+n+m]=-x[i]*s[i]+sigma*mu;
		}
		CopyMatrix(Jac,tempX);
		LUDecompose(tempX,perm,&sign);
		LinSolve(tempX,perm,residual);
		alpha=1;
		while(TRUE){
			for(i=1;i<=M;i++){
				if(i<=n){/*update x*/
					if(x[i]+alpha*residual[i]<=0)break;
				}
				if(i>n+m){/*update s*/
					if(s[i-n-m]+alpha*residual[i]<=0)break;
				}
			}
			if(i==M+1)break;
			alpha/=2;
		}
		for(i=1;i<=M;i++){
			if(i<=n){/*update x*/
				x[i]+=alpha*residual[i];
			}
			else if(i<=n+m){/*update lambda*/
				lmd[i-n]+=alpha*residual[i];
			}
			else{/*update s*/
				s[i-n-m]+=alpha*residual[i];
			}
		}
		curObj=0;
		for(i=1;i<=n;i++){
			curObj+=c[i]*x[i];
		}
		if(fabs(curObj-prevObj)<=objThresh&&sigma<=sigThresh&&rb[1]<=conThresh){
			break;
		}
		prevObj=curObj;
		sigma/=2;
	}
	if(ident)
		printf("%s, iter=%2d, initObj=%.3f, curObj=%.3f, alpha=%.3f, sigma=%.3f, constr=%.3f\n",ident,k,initObj,curObj,alpha,sigma,rb[1]);
	FreeMatrix(&gstack,Jac);
}



/* --------------------- Initialisation ---------------------- */

/* EXPORT->InitMath: initialise this module */
void InitMath(void)
{
   int i;
   double f;

   Register(hmath_version,hmath_vc_id);
   RandInit(-1);
   minLogExp = -log(-LZERO);
   numParm = GetConfig("HMATH", TRUE, cParm, MAXGLOBS);
   if (numParm>0){
      if (GetConfInt(cParm,numParm,"TRACE",&i)) trace = i;
      if (GetConfFlt(cParm,numParm,"MINEIGEN",&f)){
          min_eigen=f;
      }
      if (GetConfFlt(cParm,numParm,"NNMFTHRESHOLD",&f)){
          threshold=f;
      }
   }
}

/* ------------------------- End of HMath.c ------------------------- */
