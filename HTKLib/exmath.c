#include "HShell.h"        /* HTK Libraries */
#include "HMem.h"
#include "HMath.h"
#include "exmath.h"
/* cal mean of DMatrix type=1, mean per col, type=2: mean per row*/
DMatrix mean(DMatrix A, int type)
{
	DMatrix temp;
	int m=NumDRows(A), n=NumDCols(A),row,i,j;
	double sum;
	row = m>n?m:n;
	temp = CreateDMatrix(&gstack,row,1);
	for(i=1;i<=m;i++)
	{
		sum=0;
		for(j=1;j<=n;j++)
			{
			if(type==2)
				sum+=A[i][j];
			}

		temp[i][1]=(double)sum/n;
	}
	return temp;
}
DMatrix reshape(DMatrix A, int m, int n)
{
	DMatrix temp = CreateDMatrix(&gstack,m,n);
	int i,j;
	if(m*n==NumDRows(A)*NumDCols(A)&& NumDCols(A)==1 )
	{
		for(i=1;i<=m;i++)
			for(j=1;j<=n;j++)
				temp[i][j] = A[i+(j-1)*m][1];
	}

	return temp;
}
/* Info read *
 *
 */
char* stringcat(char *des, char * from )
{
	char* buffer;
	buffer = malloc(60*sizeof(char));
	char *filepath = buffer;

	filepath = stpcpy (filepath, des);
	filepath = stpcpy (filepath,from);
	return buffer;
}
char** TextRead(char *filepath, DVector info)
{
	FILE *fs;
	int i;
	char **arCleanInfo;
	fs = fopen (filepath, "r" ) ;
	filepath="";
	if ( fs == NULL )
	{
		puts ( "Cannot open source file" ) ;

	}

	for( i=1;i<3;i++ )
	{
		char buffer[10];
		fgets(buffer, 10,fs);
		info[i]=atof(buffer);
	}

	int n = 3*(int)info[1]+3;
	arCleanInfo= malloc(n * sizeof(char *));
	for(i=3;i<n;i++)
	{
		char* buffer = malloc(20*sizeof(char));
		fgets(buffer,20,fs);
		if(buffer==NULL)
			break;
		/*puts(buffer);*/
		/*arCleanInfo[i]=(char*) malloc(sizeof(char)*15);*/
		arCleanInfo[i]=buffer;
	}
   fclose ( fs ) ;
   return arCleanInfo;

}

/* convert string to double with e-00 format*/
double strtodouble(char *s)
{
	double temp;
	char * pEnd;
	temp=(double)strtod(s,&pEnd);
	return temp;
}
DMatrix toRow(DMatrix m)
{
	DMatrix temp;
	int i,col = NumDCols(m);
	int row = NumDRows(m);
	if(col==1)
	{
		temp = CreateDMatrix(&gstack,1,row);
		for(i=1;i<=row;i++)
			temp[1][i] = m[i][1];
		return temp;
	}else
	return m;

}
void ShowDMatrix1(char * title,DMatrix m)
{
	printf("nc=%d nr=%d",NumDCols(m),NumDRows(m));
	ShowDMatrix(title,m,NumDCols(m),NumDRows(m));

}
void initDMatrix(int start,int step, int end,DMatrix m)
{
	int i,j,nr,nc;

	   nr=NumDRows(m); nc=NumDCols(m);
	   ZeroDMatrix(m);
	   for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	      {
	    	  if(nr > nc && start+step*(i-1)<=end)
	    		  m[i][j]=start+step*(i-1);
	    	  else if(start+step*(j-1)<=end)
	    		  m[i][j]=start+step*(j-1);
	      }

}
/* Mul DMatrix with m1 dec numrowx1, m2 with dec 1xnumcol*/
DMatrix mulDMatrixRC(DMatrix m1,DMatrix m2)
{
	int i,j,nr,nc;
	DMatrix m;
	nr=NumDRows(m1); nc=NumDCols(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	   for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	      {
	    	  m[i][j]=m1[i][1]*m2[1][j];
	      }

	return m;

}
/* Mul Matrix m1 with m2, type=1: m1*m2 with nc1 = nc2*/
DMatrix mulMatrixOp(DMatrix m1,DMatrix m2)
{
	int i,j,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	if(nc1!=nc)
	 printf("cant use this func");
   for (i=1;i<=nr;i++)
   {
	   for (j=1;j<=nc;j++)
	  {

		  if(m1[i][j]!=0 && m2[j][j]!=0)
		  m[i][j]+=m1[i][j]*m2[j][j];
	  }
   }
	return m;

}
/*  this function only for trajMeanVect = trajVarMat * rq; % cq trajMeanVect = trajMeanVect'; % make it row DVector
	m1:trajVarMat(1)
	return trajMeanVect';
*/
DMatrix mulMatrixOp1(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
		DMatrix m;
		nr=NumDRows(m1);  nc=NumDCols(m2);
		nc1=NumDCols(m1);
		m= CreateDMatrix(&gstack,nc,nr);
		ZeroDMatrix(m);

	for (i=1;i<=26;i++)
	   {
		m[1][i]=m1[i][i]*m2[i][1]+m1[i][i+26]*m2[i+26][1] + m1[i][i+52]*m2[i+52][1];
	   }
	for(i=27;i<=52;i++)
		m[1][i]=m1[i][i]*m2[i][1]+m1[i][i+26]*m2[i+26][1] + m1[i][i-26]*m2[i-26][1];
	for(i=53;i<=78;i++)
			m[1][i]=m1[i][i]*m2[i][1]+m1[i][i-26]*m2[i-26][1] + m1[i][i-52]*m2[i-52][1];
		return m;

}


/* Mul Matrix m1 with m2, type=1: m1*m2*/
DMatrix mulMatrix(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);

   for (i=1;i<=nr;i++)
   {
	   for (j=1;j<=nc;j++)
	  {
		  for(k=1;k<=nc1;k++)
			  {
			  	  if(m1[i][k]!=0 && m2[k][j]!=0)
				  m[i][j]+=m1[i][k]*m2[k][j];
			  }

	  }
   }
	return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixTest(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
		DMatrix m;
		nr=NumDRows(m1);  nc=NumDCols(m2);
		nc1=NumDCols(m1);
		m= CreateDMatrix(&gstack,nr,nc);
		ZeroDMatrix(m);
		int r=0;
	   for (i=1;i<=nr;i++)
	   {
		   for (j=1;j<=nc;j++)
		  {
			  for(k=1;k<=nc1;k++)
				  {
				  	  if(m1[i][k]!=0 && m2[k][j]!=0)
				  	  {
				  		m[i][j]+=m1[i][k]*m2[k][j];
						r++;
						/*if(i>=27 && i<=52 && j>=14 && j<=26)*/
						printf("%d.m[%d][%d]+=m1[%d][%d]*m2[%d][%d] \n",r,i,j,i,k,k,j);

				  	  }

				  }

		  }
	   }
		return m;

}


/*  this function only for a tri-diagonal matrix m1 mul with m2 and the same size
*/
DMatrix mulMatrixOp3(DMatrix m1,DMatrix m2)
{
	int i,j,nr,nc,k;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);

/*Row1*/
  for (i=1;i<=26;i++)
  {
	  /*Col1*/
	  for (j=1;j<=13;j++)
		 for(k=1;k<=13;k++)
		/*  if(m1[i][k]!=0 && m2[k][j]!=0)*/
		  	m[i][j]+=m1[i][k]*m2[k][j];
	  /*Col2*/
	  for (j=27;j<=39;j++)
		 for(k=1;k<=13;k++)
		 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
		  	m[i][j]+=m1[i][k]*m2[k][j];
	  /*Col3*/
	  for (j=53;j<=65;j++)
		 for(k=1;k<=13;k++)
		  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
		  	m[i][j]+=m1[i][k]*m2[k][j];
  }

  /*Row2*/
    for (i=27;i<=52;i++)
    {
  	  /*Col1*/
  	  for (j=14;j<=26;j++)
  		 for(k=13;k<=26;k++)
  		 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
  		  	m[i][j]+=m1[i][k]*m2[k][j];
  	  /*Col2*/
  	  for (j=40;j<=52;j++)
  		 for(k=13;k<=26;k++)
  		  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
  		  	m[i][j]+=m1[i][k]*m2[k][j];
  	  /*Col3*/
  	  for (j=66;j<=78;j++)
  		m[i][j]=m[i-26][j-13];
    }
    /*Row3*/
      for (i=53;i<=78;i++)
      {
    	  /*Col1*/
    	  for (j=1;j<=13;j++)
    	  {
    		 m[i][j] = m[i-52][j+26];
    		 m[i][j+26] = m[i][j];
    		 m[i][j+52] = m[i][j];

    	  }
      }
  /*Row4*/
   for (i=79;i<=104;i++)
   {
	  /*Col1*/
	  for (j=14;j<=26;j++)
	  {
		 m[i][j] = m[i-26][j-13];
		 m[i][j+26] = m[i][j];
		 m[i][j+52] = m[i][j];

	  }
   }
  /*Row5*/
   for (i=105;i<=130;i++)
   {
	  /*Col1*/
	  for (j=1;j<=13;j++)
	  {
		 m[i][j] = m[i-104][j+52];

	  }
	  /*col2*/
	  for(j=27;j<=39;j++)
		  m[i][j] = m[i-104][j];
	  /*Col3*/
	  for(j=53;j<=65;j++)
		  m[i][j] = m[i-78][j-39];
   }
   /*Row6*/
    for (i=131;i<=156;i++)
    {
 	  /*Col1*/
 	  for (j=14;j<=26;j++)
 	  {
 		 m[i][j] = m[i-130][j+39];

 	  }
 	  /*col2*/
 	  for(j=40;j<=52;j++)
 		  m[i][j] = m[i-26][j-13];
 	  /*Col3*/
 	  for(j=66;j<=78;j++)
 		  m[i][j] = m[i-130][j-65];
    }

	return m;

}
/* Mul Matrix for dataMfcc = mulMatrixTest(dctDMatrix,dataLog);*/
DMatrix mulMatrixOp4(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);

/* row 1*/
   for (i=1;i<=13;i++)
   {
	   for (j=1;j<=26;j++)
	  {
		  for(k=1;k<=26;k++)
				  if(m1[i][k]!=0 && m2[k][j]!=0)
			  		  m[i][j]+=m1[i][k]*m2[k][j];

	  }
	   for (j=53;j<=78;j++)
		  {
			  for(k=1;k<=26;k++)
					  if(m1[i][k]!=0 && m2[k][j]!=0)
						  m[i][j]+=m1[i][k]*m2[k][j];

		  }
	   for (j=105;j<=78+52;j++)
		  {
			  for(k=1;k<=26;k++)
					  if(m1[i][k]!=0 && m2[k][j]!=0)
						  m[i][j]+=m1[i][k]*m2[k][j];

		  }
   }

   /* row 2*/
   for (i=14;i<=26;i++)
      {
   	   for (j=27;j<=52;j++)
		  {
			  for(k=27;k<=52;k++)
				  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

		  }

   	for (j=78;j<=104;j++)
   			  {
   				  for(k=27;k<=53;k++)
   						 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
   							  m[i][j]+=m1[i][k]*m2[k][j];

   			  }
   for (j=131;j<=156;j++)
	  {
		  for(k=27;k<=52;k++)
				 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

	  }

      }
   /* row 3*/
  for (i=27;i<=39;i++)
  {
   for (j=1;j<=26;j++)
   {
	  for(k=53;k<=78;k++)
			  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
				  m[i][j]+=m1[i][k]*m2[k][j];

   }
   for (j=53;j<=78;j++)
	  {
		  for(k=53;k<=78;k++)
				 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

	  }
   for (j=105;j<=78+52;j++)
	  {
		  for(k=53;k<=78;k++)
				  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

	  }
  }
  /* row 4*/
  for (i=40;i<=52;i++)
     {
  	   for (j=27;j<=52;j++)
		  {
			  for(k=79;k<=104;k++)
				 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

		  }

  	   for (j=78;j<=104;j++)
  		{
  		   m[i][j]=m[i-13][j-26];

  		}
  	   for (j=131;j<=156;j++)
  	   {
		 for(k=79;k<=104;k++)
				/*  if(m1[i][k]!=0 && m2[k][j]!=0)*/
					  m[i][j]+=m1[i][k]*m2[k][j];

  	   }

     }
  /* row 5*/
   for (i=53;i<=65;i++)
   {
    for (j=1;j<=26;j++)
    {
 	  for(k=105;k<=130;k++)
 			 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
 				  m[i][j]+=m1[i][k]*m2[k][j];

    }
    for (j=53;j<=78;j++)
 	  {
 		  for(k=105;k<=130;k++)
 				 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
 					  m[i][j]+=m1[i][k]*m2[k][j];

 	  }
    for (j=105;j<=78+52;j++)
 	  {
 		  for(k=105;k<=130;k++)
 				/*  if(m1[i][k]!=0 && m2[k][j]!=0)*/
 					  m[i][j]+=m1[i][k]*m2[k][j];

 	  }
   }

   /* row 6*/
   for (i=66;i<=78;i++)
      {
   	   for (j=27;j<=52;j++)
 		  {
 			  for(k=131;k<=156;k++)
 				/*  if(m1[i][k]!=0 && m2[k][j]!=0)*/
 					  m[i][j]+=m1[i][k]*m2[k][j];

 		  }

   	   for (j=78;j<=104;j++)
   		{
   		   for(k=131;k<=156;k++)
			 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
				  m[i][j]+=m1[i][k]*m2[k][j];

   		}
   	   for (j=131;j<=156;j++)
   	   {
   		   for(k=131;k<=156;k++)
   			/* if(m1[i][k]!=0 && m2[k][j]!=0)*/
   						  m[i][j]+=m1[i][k]*m2[k][j];

   	   }

      }
	return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixOp12(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
		DMatrix m;
		nr=NumDRows(m1);  nc=NumDCols(m2);
		nc1=NumDCols(m1);
		m= CreateDMatrix(&gstack,nr,nc);
		ZeroDMatrix(m);
		int r=0;

	   for (i=1;i<=nr;i++)
	   {

		  if(i<=13)
			   m[i][1]=m1[i][i+26]*m2[i+26][1];
		   if(i>=14 && i<=26)
			   m[i][1] = m1[i][i]*m2[i][1] + m1[i][i+52]*m2[i+52][1];
		  if(i>=27 && i<=52)
			   m[i][1] = m1[i][i-26]*m2[i-26][1] +m1[i][i]*m2[i][1] + m1[i][i+26]*m2[i+26][1];
		   if(i>=53 && i<=65)
			   m[i][1] = m1[i][i-26]*m2[i-26][1] + m1[i][i]*m2[i][1];
		   if(i>=66)
			   m[i][1]=m1[i][i]*m2[i][1];

	   }
		return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixOp11(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	int r=0,p=0;
/* row1*/
	nr = nc =13;
	nc1 = 39;
   for (i=1;i<=nr;i++)
   {
	  for (j=1;j<=nc;j++)
	  {

		  for(k=27;k<=nc1;k++)
			  {

			  	 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
			  	 	  m[i][j]+=m1[i][k]*m2[k][j];
			  }

		  m[i][j+52] = m[i][j];
	  }

	  for(j=27;j<=39;j++)
	  {
		  for(k=27;k<=nc1;k++)
			  {

				  if(m1[i][k]!=0 && m2[k][j]!=0)
					  m[i][j]+=m1[i][k]*m2[k][j];
			  }
	  }

   }

/* row2*/
   nr = nc =26;
   	nc1 = 52;
      for (i=14;i<=nr;i++)
      {
		  for (j=14;j<=nc;j++)
		  {

			  for(k=14;k<=nc1;k++)
				  {

					 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
						  m[i][j]+=m1[i][k]*m2[k][j];
				  }


		  }

		  for(j=40;j<=52;j++)
			  for(k=14;k<=nc1;k++)
				  {

					  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
						  m[i][j]+=m1[i][k]*m2[k][j];
				  }


		  for(j=66;j<=78;j++)
			  for(k=14;k<=nc1;k++)
				  {

					  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
						  m[i][j]+=m1[i][k]*m2[k][j];
				  }

      }
/* row3*/


		nc1 = 65;
		nr = 39;
		nc = 13;
		for (i=27;i<=nr;i++)
		{
		  for (j=1;j<=nc;j++)
		  {

			  for(k=1;k<=nc1;k++)
				  {

					  /*if(m1[i][k]!=0 && m2[k][j]!=0)*/
						  m[i][j]+=m1[i][k]*m2[k][j];
				  }

		  }

		  for(j=27;j<=39;j++)
			  for(k=1;k<=nc1;k++)
				  {

					 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
						  m[i][j]+=m1[i][k]*m2[k][j];
				  }

		  for(j=53;j<=65;j++)
		  			  for(k=1;k<=nc1;k++)
		  				  {

		  					 /* if(m1[i][k]!=0 && m2[k][j]!=0)*/
		  						  m[i][j]+=m1[i][k]*m2[k][j];
		  				  }
		}
/* row 4*/
	nr = 52;
	nc = 26;
	for(i=40;i<=nr;i++)
	{
		for(j=14;j<=nc;j++)
		{
			for(k=40;k<=52;k++)
				 m[i][j]+=m1[i][k]*m2[k][j];

			m[i][j+26] = m[i-39][j+13];
		}

		for(j=66;j<=78;j++)
			for(k=40;k<=52;k++)
				m[i][j]+=m1[i][k]*m2[k][j];
	}

/* row 5*/
	nr = 65;
	nc = 13;
	for(i=53;i<=nr;i++)
	{
		for(j=1;j<=nc;j++)
		{
			for(k=27;k<=65;k++)
					m[i][j]+=m1[i][k]*m2[k][j];

		}

		for(j=27;j<=39;j++)
			for(k=27;k<=65;k++)
				m[i][j]+=m1[i][k]*m2[k][j];

		for(j=53;j<=65;j++)
				for(k=27;k<=65;k++)
					m[i][j]+=m1[i][k]*m2[k][j];
	}
	/* row 5*/
	nr = 78;
	nc = 26;
	for(i=66;i<=nr;i++)
	{
		for(j=14;j<=nc;j++)
		{
			for(k=14;k<=78;k++)
				m[i][j]+=m1[i][k]*m2[k][j];
		}

		for(j=40;j<=52;j++)
			{
				for(k=14;k<=78;k++)
					m[i][j]+=m1[i][k]*m2[k][j];
			}

		for(j=66;j<=78;j++)
			{
				for(k=14;k<=78;k++)
					m[i][j]+=m1[i][k]*m2[k][j];
			}
	}

	return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixOp10(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	int r=0;
	nr = 65;
   for (i=1;i<=nr;i++)
   {
		  for(k=1;k<=nc1;k++)
			  {

	  		  m[i][1]+=m1[i][k]*m2[k][1];
			  }
		  if(i>52)
		  m[i+13][1]=m[i-52][1];

   }
	return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixOp9(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	int r=0;
	nc1 =26;
   for (i=1;i<=nr;i++)
   {

	   if(i<=13)
	   {
		   m[i][i] = m1[i][i+26]*m2[i+26][i];
		   m[i][i+26] = m1[i][i+26]*m2[i+26][i+26];
		   m[i][i+52] = m1[i][i+26]*m2[i+26][i+52];
	   }


	   if(i>13 && i<=26)
	   {
		   m[i][i] = m1[i][i]*m2[i][i] + m1[i][i+52]*m2[i+52][i];
		   m[i][i+52] = m1[i][i+52]*m2[i+52][i+52];

	   }


	   if(i>26 && i<=52)
		   m[i][i] = m1[i][i-26]*m2[i-26][i]+m1[i][i]*m2[i][i]+m1[i][i+26]*m2[i+26][i];
	   if(i>52 && i<=65)
		   m[i][i] = m1[i][i-26]*m2[i-26][i]+m1[i][i]*m2[i][i];
	   if(i>65)
		   m[i][i] = m1[i][i]*m2[i][i];

	   if(i>13 && i<=26)
	  	   {
	  		   m[i][i+26] =  m1[i][i]*m2[i][i+26] + m1[i][i+52]*m2[i+52][i+26];
	  	   }
	   if(i>26 && i<40)
	   {
		   m[i][i+26] =  m1[i][i]*m2[i][i+26] + m1[i][i+26]*m2[i+26][i+26];
		   m[i][i-26] =  m1[i][i]*m2[i][i-26];
	   }
	   if(i>=40 && i<53)
		   m[i][i+26] = m1[i][i+26]*m2[i+26][i+26];

	  if(i>=40 && i<=52)
		  m[i][i-26] = m1[i][i-26]*m2[i-26][i-26] + m1[i][i+26]*m2[i+26][i-26];
	  if(i>52 && i<=65)
		  m[i][i-26] = m1[i][i-26]*m2[i-26][i-26] + m1[i][i]*m2[i][i-26];

	  if(i>65)
		   m[i][i-26] = m1[i][i]*m2[i][i-26];


	   /* i = j+52*/
	   if(i>65)
		   m[i][i-52] = m1[i][i]*m2[i][i-52];

	   if(i<=65 && i>52)
		   m[i][i-52] = m1[i][i-26]*m2[i-26][i-52];


   }

	return m;

}
/* Mul Matrix test*/
DMatrix mulMatrixOp8(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	int r=0,p;
	j=1;
	nr=26;
	nc1 = 13;
   for (i=1;i<=nr;i++)
   {

		  for(k=1;k<=nc1;k++)
			  {
			  	  if(m1[i][k]!=0 && m2[k][j]!=0)
			  	  {

			  		  for(p=0;p<6;p++)
			  		  m[i+p*nr][j]+=m1[i][k]*m2[k+p*nc1][j];
			  		  /*m[i+13][j]+=m1[i][k]*m2[k+13][j];*/
			  		/*printf("%d.m[%d][%d]+=m1[%d][%d]*m2[%d][%d] \n",r,i,j,i,k,k,j);*/
			  	  }
			  }


   }
	return m;

}
/* Mul Matrix only for dataLog = mulMatrixTest(invDctDMatrix, dataMfcc);*/
DMatrix mulMatrixOp7(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);
	int r=0,p;
	j=1;
	nr=26;
   for (i=1;i<=nr;i++)
   {
	     for(k=1;k<=nc1;k++)
			  {


			  		  m[i][j]+=m1[i][k]*m2[k][j];
			  		/*printf("%d.m[%d][%d]+=m1[%d][%d]*m2[%d][%d] \n",r,i,j,i,k,k,j);*/

			  }
	     for(p=1;p<6;p++)
	     m[i+p*nr][j] = m[i][j];

   }
	return m;

}
/* Mul Matrix mulMatrixOp6*/
DMatrix mulMatrixOp6(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	nc1=NumDCols(m1);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);

	/*nr = nc =39;*/
   for (i=1;i<=nr;i++)
   {
	   for (j=1;j<=nc;j++)
	  {


		  m[i][j]=m2[i][j]!=0?m1[1][1]*m2[i][j]:m1[1][1]*m2[i+39][j];

	  }
   }
	return m;

}

 /*  this function only for a diagonal matrix m1 mul with m2 and the same size*/
DMatrix mulMatrixOp2(DMatrix m1,DMatrix m2)
{
	int i,j,nr,nc,k;
	DMatrix m;
	nr=NumDRows(m1);  nc=NumDCols(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	ZeroDMatrix(m);

/*int r=0;*/
   for (i=1;i<=nr;i++)
   {
	   for (j=i;j<=nc;j++)
	  {

		 if(m1[j][j]!=0 && m2[i][j]!=0)
			 m[i][j]=m1[j][j]*m2[i][j];
			  /*printf("%d.m[%d][%d]+=%f*%f \n",r,i,j,m1[i][j],m2[j][j]);*/

		/* if(m[i][j]!=0)*/
			 m[j][i] = m[i][j];


	  }
   }
	return m;

}

/*  this function only for //% cepstral liftering dataMfcc = cepDMatrix * dataMfcc;
	m1:dataMfcc(1)
	return dataMfcc;
*/
DMatrix mulMatrixOp5(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,nc1;
		DMatrix m;
		nr=NumDRows(m1);  nc=NumDCols(m2);
		nc1=NumDCols(m1);
		m= CreateDMatrix(&gstack,nr,nc);
		ZeroDMatrix(m);

		/*int r=0;*/
	   for (i=1;i<=nr;i++)
	 	   m[i][1]=m1[i][i]*m2[i][1];

		return m;

}
/* Mul DMatrix m1 with m2:m1'*m2*/
DMatrix mulMatrix2(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr1,nr,nc;
	DMatrix m;
	nr1=NumDRows(m1);
	nr = NumDCols(m1);
	nc = NumDCols(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	int r=0;
	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
		  for(k=1;k<=nr1;k++)
			  if(m1[k][i]!=0 &&m2[k][j]!=0)
			  {
				  r++;
			  m[i][j]+=m1[k][i]*m2[k][j];
			/*  printf("%d.m[%d][%d]+=m1[%d][%d]*m2[%d][%d]%f*%f=%f\n",r,i,j,k,i,k,j,m1[j][i],m2[j][j],m[i][j]);*/
			  }

	return m;

}
/* Mul DMatrix m1 with m2:m1'*m2*/
DMatrix mulMatrix2Test(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr1,nr,nc;
	DMatrix m;
	nr1=NumDRows(m1);
	nr = NumDCols(m1);
	nc = NumDCols(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	int r=0;
	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
		  for(k=1;k<=nr1;k++)
			  if(m1[k][i]!=0 &&m2[k][j]!=0)
			  {
				  r++;
			  m[i][j]+=m1[k][i]*m2[k][j];

			 /* printf("%d.m[%d][%d]+=m1[%d][%d]*m2[%d][%d]%f*%f=%f\n",r,i,j,k,i,k,j,m1[k][i],m2[k][j],m[i][j]);*/

			  }

	return m;

}
/* Mul DMatrix m1 with m2:m1'*m2 with nr1=nc2*/
DMatrix mulMatrix2Op(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr1,nr,nc;
	DMatrix m;
	nr1=NumDRows(m1);
	nr = NumDCols(m1);
	nc =nr;
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	int r=0;
	/* row1*/
	for(i=1;i<=13;i++)
	{
		m[i][i+26]= m1[i+26][i]*m2[i+26][1];

	}
	/* row2*/
	for(i=14;i<=26;i++)
	{
		m[i][i]= m1[i][i]*m2[i][1];
		m[i][i+52]= m1[i+52][i]*m2[i+52][1];

	}
	/* row3*/
	for(i=27;i<=52;i++)
	{
		m[i][i]= m1[i][i]*m2[i][1];
		m[i][i-26]= m1[i-26][i]*m2[i-26][1];
		m[i][i+26]= m1[i+26][i]*m2[i+26][1];
	}
	/* row4*/
	for(i=53;i<=65;i++)
	{
		m[i][i]= m1[i][i]*m2[i][1];
		m[i][i-26]= m1[i-26][i]*m2[i-26][1];

	}
	/* row4*/
	for(i=66;i<=78;i++)
	{
		m[i][i]= m1[i][i]*m2[i][1];


	}
	return m;

}

/* Mul DMatrix m1 with m2:m1'*m2 with nr1=nc2*/
DMatrix mulMatrix2Op_t(DMatrix m1,DMatrix m2)
{

	int i,j,k,nr1,nr,nc;
		DMatrix m;
		nr1=NumDRows(m1);
		nr = NumDCols(m1);
		nc =nr;
		m= CreateDMatrix(&gstack,nr,nc);

		ZeroDMatrix(m);
		int r=0;
		/* row1*/
		for(i=1;i<=13;i++)
		{
			m[i][i+26]= m1[i+26][i]*m2[i+26][1];

		}
		/* row2*/
		for(i=14;i<=26;i++)
		{
			m[i][i]= m1[i][i]*m2[i][1];
			m[i][i+52]= m1[i+52][i]*m2[i+52][1];

		}
		/* row3*/
		for(i=27;i<=52;i++)
		{
			m[i][i]= m1[i][i]*m2[i][1];
			m[i][i-26]= m1[i-26][i]*m2[i-26][1];
			m[i][i+26]= m1[i+26][i]*m2[i+26][1];
		}
		/* row4*/
		for(i=53;i<=65;i++)
		{
			m[i][i]= m1[i][i]*m2[i][1];
			m[i][i-26]= m1[i-26][i]*m2[i-26][1];

		}
		/* row4*/
		for(i=66;i<=78;i++)
		{
			m[i][i]= m1[i][i]*m2[i][1];


		}
		return m;

}

/* Mul DMatrix m1 with m2:m1'*m2 with m1=m2 size 1,156*/
DMatrix mulMatrix2Op1(DMatrix m1)
{
	int i,j,k,nr1,nr,nc;
	DMatrix m;
	nr = NumDCols(m1);
	nc = nr;
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	  for (j=i;j<=nc;j++)
	  {

		m[i][j]=m1[1][i]*m1[1][j];

		m[j][i] = m[i][j];

	  }


	return m;

}
/* Mul DMatrix m1 with m2:m1*m2'*/
DMatrix mulMatrix3(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
		  for(k=1;k<=nc1;k++)
			  if(m1[i][k]!=0 &&m2[j][k]!=0)
			  m[i][j]+=m1[i][k]*m2[j][k];

	return m;

}
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Test(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	int r=0;

	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
		  for(k=1;k<=nc1;k++)
			  if(m1[i][k]!=0 &&m2[j][k]!=0)
			  {
				  r++;
				  m[i][j]+=m1[i][k]*m2[j][k];
				  if(i>=105 && i<=130 &&j>=105 && j<=130 )
				  printf("%d.m[%d][%d]=m1[%d][%d]*m2[%d][%d] \n",r,i,j,i,k,j,k);
			  }


	return m;

}
/* Mul DMatrix m1 with m2:m1*m2' j=1, nc1 = nr1*/
DMatrix mulMatrix3Op2(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc,p;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	/*Row1*/
	for (i=1;i<=26;i++)
	{
	/* col1*/
	  for (j=1;j<=26;j++)
		  for(k=1;k<=13;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
			     m[i][j]+=m1[i][k]*m2[j][k];
	/* col2*/
	  for (j=53;j<=78;j++)
		  for(k=27;k<=39;k++)
			  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 m[i][j]+=m1[i][k]*m2[j][k];
	  /* col3*/
	  for (j=105;j<=130;j++)
		  for(k=53;k<=65;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 m[i][j]+=m1[i][k]*m2[j][k];
	}

	/*Row2*/
	for (i=27;i<=52;i++)
	{
	/* col1*/
	  for (j=27;j<=52;j++)
		  for(k=14;k<=26;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 m[i][j]+=m1[i][k]*m2[j][k];
	  /* col2*/
	  for (j=79;j<=104;j++)
		  m[i][j] = m[i-26][j-26];

	  /* col3*/
	  	  for (j=131;j<=156;j++)
	  		  m[i][j] = m[i-26][j-26];
	}
	/*Row3 & 4*/
	for (i=53;i<=78;i++)
	{
	/* col1*/
	  for (j=1;j<=26;j++)
	  {
		  m[i][j]=m[i-52][j+52];
		  m[i][j+52] = m[i][j];
		  m[i][j+104] = m[i][j];
		  m[i+26][j+26] = m[i][j];
		  m[i+26][j+78] = m[i][j];
		  m[i+26][j+130] = m[i][j];
		  m[i+52][j+52] = m[i][j];
	  }


	}
/* Row5*/
	for (i=105;i<=130;i++)
	{
		for(j=1;j<=26;j++)
		  for(k=1;k<=13;k++)
			  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
			     m[i][j]+=m1[i][k]*m2[j][k];

	  /* col3*/
	  for (j=105;j<=130;j++)
		  for(k=53;k<=65;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 m[i][j]+=m1[i][k]*m2[j][k];
	}

	/* Row6*/
	for (i=131;i<=156;i++)
	{
		 for (j=27;j<=52;j++)
		 {
			 m[i][j] = m[i-26][j-26];
			 m[i][j+52] = m[i-26][j+26];
			 m[i][j+104] = m[i-130][j-26];
		 }

	}

	return m;

}
/* Mul DMatrix3 only for dataMfcc = mulMatrix3(dataMfcc,dctDMatrix);*/
DMatrix mulMatrix3Op3(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
/*Row1*/
	for (i=1;i<=13;i++)
	{
		for (j=1;j<=13;j++)
		  for(k=1;k<=26;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
			 	  m[i][j]+=m1[i][k]*m2[j][k];
		/* col2*/
		for(j=27;j<=39;j++)
		  for(k=53;k<=78;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];
		/* col3*/
		for(j=53;j<=65;j++)
		  for(k=105;k<=130;k++)
			  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];

	}
/* Row2*/
	for (i=14;i<=26;i++)
	{
		/*Col1*/
		for (j=14;j<=26;j++)
		  for(k=26;k<=52;k++)
			/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
			 	  m[i][j]+=m1[i][k]*m2[j][k];
		/* col2*/
		for(j=40;j<=52;j++)
		  for(k=78;k<=104;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];
		/* col3*/
		for(j=66;j<=78;j++)
		  for(k=130;k<=156;k++)
			  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];

	}

	/*Row3*/
		for (i=27;i<=39;i++)
		{
			for (j=1;j<=13;j++)
			  for(k=1;k<=26;k++)
				/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 	  m[i][j]+=m1[i][k]*m2[j][k];
			/* col2*/
			for(j=27;j<=39;j++)
			  for(k=53;k<=78;k++)
				/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col3*/
			for(j=53;j<=65;j++)
			  for(k=105;k<=130;k++)
				/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];

		}
	/* Row4*/
		for (i=40;i<=52;i++)
		{
			/*Col1*/
			for (j=14;j<=26;j++)
			  for(k=26;k<=52;k++)
				 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col2*/
			for(j=40;j<=52;j++)
			  for(k=78;k<=104;k++)
				  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col3*/
			for(j=66;j<=78;j++)
			  for(k=130;k<=156;k++)
				 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];

		}
	/*Row5*/
		for (i=53;i<=65;i++)
		{
			for (j=1;j<=13;j++)
			  for(k=1;k<=26;k++)
				 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col2*/
			for(j=27;j<=39;j++)
			  for(k=53;k<=78;k++)
				 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col3*/
			for(j=53;j<=65;j++)
			  for(k=105;k<=130;k++)
				  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];

		}
	/* Row6*/
		for (i=66;i<=78;i++)
		{
			/*Col1*/
			for (j=14;j<=26;j++)
			  for(k=26;k<=52;k++)
				 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col2*/
			for(j=40;j<=52;j++)
			  for(k=78;k<=104;k++)
				/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
			/* col3*/
			for(j=66;j<=78;j++)
			  for(k=130;k<=156;k++)
				/*  if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];

		}
	return m;

}
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Op9(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	/*int r=0;*/
	nr=13;
	nc=39;
	nc1 = 39;
	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
	  {
		  for(k=27;k<=nc1;k++)
			 /* if(m1[i][k]!=0 &&m2[j][k]!=0)*/
			   m[i][j]+=m1[i][k]*m2[j][k];

		  if(j>=26)
			  m[i][j+26] = m[i][j];
	  }

	nr=26;
	nc=26;
	nc1 = 52;
	for (i=14;i<=nr;i++)
	{
		for (j=14;j<=nc;j++)
		{
		  for(k=14;k<=nc1;k++)
			  if(m1[i][k]!=0 &&m2[j][k]!=0)

				  m[i][j]+=m1[i][k]*m2[j][k];
		  m[i][j+26] = - m[i-13][j+13];
		}

		for (j=66;j<=78;j++)
			{
			  for(k=14;k<=nc1;k++)
				   m[i][j]+=m1[i][k]*m2[j][k];

			}

	}

	nr=39;
	nc=39;
	nc1 = 65;
	for (i=27;i<=nr;i++)
	{
		for (j=27;j<=nc;j++)
		{
		  for(k=1;k<=nc1;k++)
			  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];

		  m[i][j-26] = m[i-26][j];
		}

		for (j=53;j<=65;j++)
				{
				  for(k=27;k<=nc1;k++)
					  /*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
					  m[i][j]+=m1[i][k]*m2[j][k];
				}
	}
/* row4*/
	for(i=40;i<=52;i++)
		for(j=14;j<=26;j++)
		{
			m[i][j] = -m[i-13][j-13];
			m[i][j+26] = m[i-39][j-13];
			m[i][j+52] = m[i-13][j-13];
		}
	/* row5*/
	for(i=53;i<=65;i++)
	{
		for(j=1;j<=13;j++)
		{
			m[i][j]= m[i-26][j];
			m[i][j+26]= m[i-26][j+52];

		}
		for(j=53;j<=65;j++)
			for(k=27;k<=65;k++)
			{
				/*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				  m[i][j]+=m1[i][k]*m2[j][k];
			}



	}/* end for i row 5*/

	/* Row 6*/

	for(i=66;i<=78;i++)
	{
		for(j=14;j<=26;j++)
		{
			m[i][j] = m[i-52][j+52];
			m[i][j+26] = m[i-13][j-13];
		}
		for(j=66;j<=78;j++)
		{
			for(k=14;k<=78;k++)
				/*if(m1[i][k]!=0 &&m2[j][k]!=0)*/
				 m[i][j]+=m1[i][k]*m2[j][k];
		}
	}



	return m;

}
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Op8(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	  	m[i][1]+=m1[i][i]*m2[1][i];

	return m;

}
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Op7(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	int p=0;
	nc1 = 26;
	for (i=1;i<=nr;i++)
	{

	   for(k=1;k<=nc1;k++)
		  m[i][1]+=m1[i][k+p*nc1]*m2[1][k+p*nc1];

	   if(i%13==0)
		   p++;
	}
	return m;

}
/* Mul DMatrix3 Test*/
DMatrix mulMatrix3Op6(DMatrix m1,DMatrix m2)
{
	int i,j,k,nc1,nr,nc;
	DMatrix m;
	nc1=NumDCols(m1);
	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);
	nr  =39;
	ZeroDMatrix(m);

	for (i=1;i<=nr;i++)
	{
		for(k=1;k<=65;k++)
		  /*if(m1[i][k]!=0 &&m2[1][k]!=0)*/

			 m[i][1]+=m1[i][k]*m2[1][k];
		if(i<=13)
		m[i+nr][1] = m[i][1];

	}
	nr  =78;
	for (i=53;i<=nr;i++)
		{
			for(k=13;k<=nc1;k++)
			  /*if(m1[i][k]!=0 &&m2[1][k]!=0)*/

				 m[i][1]+=m1[i][k]*m2[1][k];
				 /* printf("%d.m[%d][%d]=m1[%d][%d]*m2[%d][%d] \n",r,i,1,i,k,1,k);*/


		}

	return m;

}
/* Mul DMatrix3 Only for dataMfcc = mulMatrix3Op5(invCepDMatrix,dataMfcc);*/
DMatrix mulMatrix3Op5(DMatrix m1,DMatrix m2)
{
	int nr,i;
	DMatrix m;
	nr = NumDRows(m1);
	m= CreateDMatrix(&gstack,nr,1);

	ZeroDMatrix(m);
	for (i=1;i<=nr;i++)
		m[i][1]+=m1[i][i]*m2[1][i];


	return m;

}
/* only for m2 is factmatrix*/
DMatrix mulMatrix3Op4(DMatrix m1,DMatrix m2)
{
	int i,j,nr,nc;
	DMatrix m;

	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	/*int r=0;*/

	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
	  m[i][j]=m1[i][j]+m1[i][j+39];




	return m;

}


/* Mul Matrix m1 with m2:m1*m2' with nc1 = nr2*/
DMatrix mulMatrix3Op(DMatrix m1,DMatrix m2)
{
	int i,j,nc1,nr,nc;
	DMatrix m;

	nr = NumDRows(m1);
	nc = NumDRows(m2);
	m= CreateDMatrix(&gstack,nr,nc);

	ZeroDMatrix(m);
	/*int r=0;*/
	for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
	  {
		  if(m1[i][j]!=0 &&m2[j][j]!=0)
		  {
			  /*r++;*/
			  m[i][j]+=m1[i][j]*m2[j][j];
			 /* printf("%d.m[%d][%d]=%f*%f \n",r,i,j,m1[i][j],m2[j][j]);*/
		  }

	  }
	return m;

}
/* Mul Matrix m1 with m2:m1*m2' with nc1 = 2*nr2*/
DMatrix mulMatrix3Op1(DMatrix m1,DMatrix m2)
{
	int i,j,k,nr,nc,p;
		DMatrix m;

		nr = NumDRows(m1);
		nc = NumDRows(m2);
		m= CreateDMatrix(&gstack,nr,nc);

		ZeroDMatrix(m);

		for (i=1;i<=26;i++)
		  for (j=1;j<=26;j++)
		  {
			  for(k=1;k<=26;k++)
				   m[i][j]+=m1[i][k]*m2[j][k];

			  for(p=1;p<6;p++)
				  m[i+p*26][j+p*26] = m[i][j];
		  }

		return m;

}

/* Mul DMatrix m1 with m2, type=1: m1*m2, type=2:m1'*m2, type=3:m1*m2'*/
DMatrix mulDMatrix(DMatrix m1,DMatrix m2,int type)
{
	int i,j,k,nr1,nc1,nr2,nc2,nr,nc;
	DMatrix m;
	nr1=NumDRows(m1); nc1=NumDCols(m1);
	nr2=NumDRows(m2); nc2=NumDCols(m2);
	if (type==1 && nc1==nr2)
	{
		nr=nr1; nc=nc2;
		m= CreateDMatrix(&gstack,nr1,nc2);
	}
	else if (type==2 && nr1==nr2){
		nr=nc1; nc=nc2;
		m= CreateDMatrix(&gstack,nr,nc);

	}else if(type==3)
	{
		nr=nr1; nc=nr2;
		m= CreateDMatrix(&gstack,nr,nc);

	}
	ZeroDMatrix(m);

   for (i=1;i<=nr;i++)
	  for (j=1;j<=nc;j++)
	  {
		  if(type==1)
			  for(k=1;k<=nc1;k++)
			  {
				  m[i][j]+=m1[i][k]*m2[k][j];
			  }
		  else if(type==2)
			  for(k=1;k<=nr1;k++)
			  {
			 		m[i][j]+=m1[k][i]*m2[k][j];
			  }
		  else if(type==3)
		  {

			  for(k=1;k<=nc1;k++)
			  {
					m[i][j]+=m1[i][k]*m2[j][k];
			  }
		  }
	  }

	return m;

}
/*cal DMatrix m1 with m2:m=m1+m2*/
DMatrix sumEDMatrix1(DMatrix m1,DMatrix m2)
{
	DMatrix temp;
	int nr1,nc1,nr2,nc2;
	int i,j;

	nr1 = NumDRows(m1); nc1 = NumDCols(m1);
	temp = CreateDMatrix(&gstack,nr1,nc1);

	for(i=1;i<=nr1;i++)
		for(j=1;j<=nc1;j++)
			temp[i][j] = m1[i][j]+m2[i][j];
return temp;
}

/* DMatrix m1 with m2:m=m1-m2*/
DMatrix minusEDMatrix(DMatrix m1,DMatrix m2)
{
	DMatrix temp;
	int nr1,nc1,nr2,nc2;
	int i,j;

	nr1 = NumDRows(m1); nc1 = NumDCols(m1);
	temp = CreateDMatrix(&gstack,nr1,nc1);

	for(i=1;i<=nr1;i++)
		for(j=1;j<=nc1;j++)
			temp[i][j] = m1[i][j]-m2[i][j];
return temp;
}
/* DMatrix m1 with m2:m=m1.*m2*/
DMatrix mulEDMatrix(DMatrix m1,DMatrix m2)
{
	DMatrix temp;
	int nr1,nc1,nr2,nc2;
	int i,j;

	nr1 = NumDRows(m1); nc1 = NumDCols(m1);
	temp = CreateDMatrix(&gstack,nr1,nc1);

	/*ZeroDMatrix(temp);*/
	for(i=1;i<=nr1;i++)
		for(j=1;j<=nc1;j++)
			temp[i][j] = m1[i][j]*m2[i][j];
return temp;
}
/* DMatrix m1 with m2:m=m1/m2*/
DMatrix divEDMatrix(DMatrix m1,DMatrix m2)
{
	DMatrix temp;
	int nr1,nc1,nr2,nc2;
	int i,j;

	nr1 = NumDRows(m1); nc1 = NumDCols(m1);
	temp = CreateDMatrix(&gstack,nr1,nc1);

	/*ZeroDMatrix(temp);*/
	for(i=1;i<=nr1;i++)
		for(j=1;j<=nc1;j++)
			temp[i][j] = m1[i][j]/m2[i][j];
return temp;
}
/* cal DMatrix m1 with m2, type=1: m1+m2, type =2:m1.*m2, type=3:m1./m2, type =4:m1-m2*/
DMatrix calDMatrix(DMatrix m1,DMatrix m2,int type)
{
	DMatrix temp;
	int nr1,nc1,nr2,nc2;
	int i,j;

	nr1 = NumDRows(m1); nc1 = NumDCols(m1); nr2 = NumDRows(m2); nc2 = NumDCols(m2);
	temp = CreateDMatrix(&gstack,nr1,nc1);
	if(nr1!=nr2 && nc1!=nc2)
	{
		puts("problem with size of two DMatrix");
		return temp;
	}
	/*ZeroDMatrix(temp);*/
	for(i=1;i<=nr1;i++)
		for(j=1;j<=nc1;j++)
			if(type==1)
				temp[i][j] = m1[i][j]+m2[i][j];
			else if (type==2)
				temp[i][j] = m1[i][j]*m2[i][j];
			else if (type==3)
				temp[i][j] = m1[i][j]/m2[i][j];
			else if (type==4)
				temp[i][j] = m1[i][j]-m2[i][j];
return temp;
}

/* Mul a Number with a DMatrix m*/
DMatrix mulNumberwDMatrix(double num, DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=num*m[i][j];
	return temp;
}
/* Sum a Number with a DMatrix m*/
DMatrix sumNumberwDMatrix(double num, DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
	/*ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	   	  temp[i][j]=num+m[i][j];
	return temp;
}

/* sum all element of DMatrix m*/
double sumEDMatrix(DMatrix m)
{
	double sum=0;
	int i,j,nr,nc;
	nr=NumDRows(m); nc=NumDCols(m);
	for (i=1;i<=nr;i++)
		 for (j=1;j<=nc;j++)
			 if(m[i][j]!=0)
		  	 sum+= m[i][j];
	return sum;

}
/* Divide a Number with a DMatrix m*/
DMatrix divNumberwDMatrix(double num, DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
/*	ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=num/m[i][j];
	return temp;
}


/* Cos a DMatrix*/
DMatrix cosDMatrix(DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
/*	ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=cos(m[i][j]);
	return temp;
}

/* Sin a DMatrix*/
DMatrix sinDMatrix(DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
	/*ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=sin(m[i][j]);
	return temp;
}
/* Exp a DMatrix*/
DMatrix expDMatrix(DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
/*	ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=exp(m[i][j]);
	return temp;
}
/* pow a DMatrix*/
DMatrix powDMatrix(double num,DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
	/*ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=pow(m[i][j],num);
	return temp;
}

/* log a DMatrix*/
DMatrix logDMatrix(DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m); nc=NumDCols(m);
	temp= CreateDMatrix(&gstack,nr,nc);
	/*ZeroDMatrix(temp);*/
	for (i=1;i<=nr;i++)
	      for (j=1;j<=nc;j++)
	    	  temp[i][j]=log(m[i][j]);
	return temp;
}
/* calc a number with element of DMatrix
 * type is type of cal s: sin, c:cos, d divide, t: sum, m: multiply,p: pow, e:exp
 * type ='l':log
 */
DMatrix calNumWDMatrix(char type,double num,DMatrix m)
{
	int i,j,nr,nc;
		DMatrix temp;
		nr=NumDRows(m); nc=NumDCols(m);
		temp= CreateDMatrix(&gstack,nr,nc);
		ZeroDMatrix(temp);
		for (i=1;i<=nr;i++)
		      for (j=1;j<=nc;j++)
		    	  switch(type)
		    	  {
		    	  case 's' :
		    		  /*if(m[i][j]!=0)*/
		    		  temp[i][j]=sin(m[i][j]);
		    		  break;
		    	  case 'c':
		    		  temp[i][j]=cos(m[i][j]);
		    		  break;
		    	  case 'd':
		    		  temp[i][j]=num/m[i][j];
		    		  break;
		    	  case 't':
		    		  temp[i][j]=num+m[i][j];
		    		  break;
		    	  case 'm':
		    		  if(m[i][j]!=0)
		    		  temp[i][j]=num*m[i][j];
		    		  break;
		    	  case 'p':
		    		  if(m[i][j]!=0)
		    		  temp[i][j]=pow(m[i][j],num);
		    		  break;
		    	  case 'e':
					  temp[i][j]=exp(m[i][j]);
					  break;
		    	  case 'l':
					  temp[i][j]=log(m[i][j]);
					  break;
		    	  }

		return temp;
}

/* calc a number with element of DVector
 * type is type of cal s: sin, c:cos, d divide, t: sum, m: multiply,p: pow
 */
DVector calNumWDVector(char type,double num,DVector v)
{
	int i,nc;
	DVector temp;
	nc=DVectorSize(v);
	temp= CreateDVector(&gstack,nc);
	ZeroDVector(temp);
	for (i=1;i<=nc;i++)
			  switch(type)
			  {
			  case 's' :
				  temp[i]=sin(v[i]);
				  break;
			  case 'c':
				  temp[i]=cos(v[i]);
				  break;
			  case 'd':
				  temp[i]=num/v[i];
				  break;
			  case 't':
				  temp[i]=num+v[i];
				  break;
			  case 'm':
				  temp[i]=num*v[i];
				  break;
			  case 'p':
				  temp[i]=pow(v[i],num);
				  break;
			  }

	return temp;
}
DMatrix eye(int nr, int nc)
{
	int i,j;
		DMatrix temp;
		temp= CreateDMatrix(&gstack,nr,nc);
		ZeroDMatrix(temp);
		for (i=1;i<=nr;i++)
		      for (j=1;j<=nc;j++)
		    	 if(i==j)
		    		 temp[i][j]=1;
		return temp;
}

/* Kronecker tensor product*/
DMatrix kron(DMatrix a, DMatrix b)
{
	int i,j,ii,jj,nc,nr,nra,nca,nrb,ncb;
	int k,q;
	nra=NumDRows(a); nrb=NumDRows(b);nca= NumDCols(a); ncb =NumDCols(b);
	nr=nra*nrb; nc=nca*ncb;
	DMatrix temp;
	temp= CreateDMatrix(&gstack,nr,nc);
	for (i=1;i<=nra;i++)
		  for (j=1;j<=nca;j++)
			  for (ii=1;ii<=nrb;ii++)
				  for (jj=1;jj<=ncb;jj++)
				  {
					  k = nrb*(i-1)+ii;
					  q = ncb*(j-1)+jj;
					  temp[k][q]=a[i][j]*b[ii][jj];
				  }
	return temp;
}

/* return  DMatrix col column form m DMatrix*/
DMatrix cutColDMatrix(DMatrix m, int col)
{
	int i,j,nr;
	DMatrix temp;
	nr=NumDRows(m);
	temp= CreateDMatrix(&gstack,nr,col);
	for (i=1;i<=nr;i++)
			  for (j=1;j<=col;j++)
				  temp[i][j]=m[i][j];
	return temp;
}
/* diag of DMatrix m*/
DMatrix diag(DMatrix m)
{
	int i,j,nr,nc;
	DMatrix temp;
	nr=NumDRows(m);
	nc=NumDCols(m);

	if(nc!=1 && nr!=1)
	{
		temp= CreateDMatrix(&gstack,nr,1);
		ZeroDMatrix(temp);
		for(i=1;i<=nc;i++)
			for(j=1;j<=nr;j++)
				if(i==j)
				temp[i][1] = m[i][j];

		return temp;
	}
	if(nr==1)
		nr=nc;
	if(nc==1)
		nc=nr;
	temp= CreateDMatrix(&gstack,nr,nc);
 int nc1 = NumDCols(m);
 int nr1 = NumDRows(m);
	for (i=1;i<=nr;i++)
			  for (j=1;j<=nr;j++)
				  if(i==j && nc1==1 )
				  temp[i][j]=m[i][1];
				  else if (i==j && nr1==1 )
					  temp[i][j]=m[1][j];

	return temp;
}
/* Convert ones matlab function to C*/
DMatrix ones(int nr, int nc)
{
	int i,j;
			DMatrix temp;
			temp= CreateDMatrix(&gstack,nr,nc);
			for (i=1;i<=nr;i++)
			      for (j=1;j<=nc;j++)

			    		 temp[i][j]=1;
			return temp;

}
DMatrix* initializeMatrices(int numObsDVectors, int numDctDVectors,int numChans,int numceps,int cepLifter,int M,int D,int windowL,char *option, char* wOpt)
{


	DMatrix dctDMatrix,cepDMatrix,invDctDMatrix,invCepDMatrix,W,VS,W0;
	DMatrix *initMatrices;
	DMatrix *dctMatrices;
	/*DMatrix Matrices[8];*/
	DMatrix* Matrices=(DMatrix*)New(&gstack,sizeof(DMatrix)*8);
	/*%% Get DCT matrices %%*/

	dctMatrices = initializeDctMatrices(numDctDVectors, numChans, numceps, cepLifter);

	 dctDMatrix = dctMatrices[1];
	 cepDMatrix = dctMatrices[2];
	 invDctDMatrix = dctMatrices[3];
	 invCepDMatrix = dctMatrices[4];

	/*%% Selection DMatrix %%
		Sq = kron(ones(numObsDVectors, 1), eye(D*M)); % used for all trajectory based methods
	*/
	/*Sq = kron(ones(numObsDVectors, 1),eye(D*M,D*M));
	/*
	%% Trajectory weight DMatrix %%
		W = [];
		if ~isempty(regexp(option, 'traj'))
			fprintf(1,'van cu vao');
			[W, W0] = getDynamicFeatureDMatrix(M, D, windowL, numObsDVectors, wOpt);
			[NumDRows, NumDCols] = size(W0);
			if (NumDRows ~= D*numObsDVectors) || (NumDCols ~= numDctDVectors)
				fprintf(1, '(NumDRows %d ~= D*numObsDVectors %d) || (NumDCols %d ~= numDctDVectors %d)\n', NumDRows, D*numObsDVectors, NumDCols, numDctDVectors);
				return;
			end
		end
*/

	 if(strncmp(option, "traj", 4)==0)
		 {

		 initMatrices = getDynamicFeatureDMatrix( M,  D,  windowL,  numObsDVectors,  wOpt);

		 W = initMatrices[1];
		 W0 = initMatrices[2];


			 if(NumDRows(W0)!=D*numObsDVectors ||NumDCols(W0)!=numDctDVectors)
			 {
				 printf("(NumDRows %d != D*numObsDVectors %d) || (NumDCols %d != numDctDVectors %d)\n", NumDRows(W0), D*numObsDVectors, NumDCols(W0), numDctDVectors);
				 return NULL;
			 }
			/* %% VS DMatrix for ML optimization %%*/
			 /* Thang's matlab code is not working with trajAvg-tML Option*/
			 if(strcmp(option, "trajAvg-tML")==0)
			 		VS = getVSMat(M, D, numDctDVectors, W, W0);



			/* ShowDMatrix1("DMatrix",W0);*/

		 }

	 Matrices[1]=dctDMatrix;
	 Matrices[2]=cepDMatrix;
	 Matrices[3]=invDctDMatrix;
	 Matrices[4]=invCepDMatrix;
	 Matrices[6]=W;
	 Matrices[7]=VS;

	 return Matrices;

}
/* VS DMatrix for ML optimization*/
DMatrix getVSMat(int M,int D,int numDctDVectors,DMatrix W,DMatrix W0)
{
	DMatrix VS;

	return VS;

}
/*get Dynamic Feature DMatrix*/
/*M: num of static parameters, D: order of dynamic parameters*/
/*Typical setup: M=13, D=3, windowL=1, T=2*/
DMatrix* getDynamicFeatureDMatrix(int M, int D, int windowL, int T, char* wOpt)
{
	int rangeLength = 2;
	DMatrix delCoeffs,accCoeffs,coefficients;
	DMatrix W, W0;
	/*DMatrix FeatureDMatrix[3];*/
	DMatrix* FeatureDMatrix=(DMatrix*)New(&gstack,sizeof(DMatrix)*3);
	int i,k,h,j,ii,jj;
	int startRow,endRow,startCol,endCol;
	delCoeffs = CreateDMatrix(&gstack,1,windowL*2+1);
	accCoeffs = CreateDMatrix(&gstack,1,2*rangeLength*windowL + 1);
	coefficients = CreateDMatrix(&gstack,3,2*rangeLength*windowL+1);
	W = CreateDMatrix(&gstack,T*D, T + 2*rangeLength*windowL);
	ZeroDMatrix(accCoeffs);
	/*% Assume D = 3*/
	/*if(D!=3) FeatureDMatrix;*/
	initDMatrix(-windowL,1,windowL,delCoeffs);
	delCoeffs = calNumWDMatrix('m',1/sumEDMatrix(calNumWDMatrix('p',2,delCoeffs)),delCoeffs);

	for(k=1;k<=2*windowL+1;k++)
		for(h=1;h<=2*windowL+1;h++)
		{
			i=k+h-1;
			accCoeffs[1][i]+= delCoeffs[1][k]*delCoeffs[1][h];
		}

	coefficients[1][2*windowL+1]=1;
	for(i=windowL+1,k=1;i<=3*windowL+1;i++,k++)
		coefficients[2][i] = delCoeffs[1][k];

	for(i=1,k=1;i<=2*rangeLength*windowL+1;i++,k++)
		coefficients[3][i] = accCoeffs[1][k];
	/*DMatrix W1 = coefficients;*/
	/*% construct full W DMatrix*/
	for (k=1;k<=T;k++)
	{
		startRow = (k-1)*D + 1;
	   	endRow = k*D;
	   	startCol = k;
	   	endCol = k + 2*rangeLength*windowL;
	   	for(i=startRow,ii=1;i<=endRow;i++,ii++)
	   		for(j=startCol,jj=1;j<=endCol;j++,jj++)
	   			W[i][j] = coefficients[ii][jj];
	}
	/*ShowDVector("DVector",coefficients[1],6);*/
	W0=W;
	W = kron(W0, eye(M,M));
	FeatureDMatrix[1] = W;
	FeatureDMatrix[2] = W0;
	return FeatureDMatrix;
	/*ShowDMatrix1("DMatrix",W0);*/



}
/*function [dctDMatrix, cepDMatrix, invDctDMatrix, invCepDMatrix]*/
DMatrix*  initializeDctMatrices(int numDctDVectors, int numChans, int numceps, int cepLifter)
{

	DMatrix temp,temp1,temp2;
	DMatrix dctDMatrix,cepDMatrix,invDctDMatrix,invCepDMatrix;
	/*DMatrix dctMatrices[5];*/
	DMatrix* dctMatrices=(DMatrix*)New(&gstack,sizeof(DMatrix)*5);
	/*% DCT DMatrix*/
  /*dctDMatrix = sqrt(2.0/numChans) * cos ( (PI/(2*numChans)) * (0:1:numceps)' * (1:2:(2*numChans-1)) );*/
  /*dctDMatrix = kron(eye(numDctDVectors), dctDMatrix);/* % for trajectory model if numDctDVectors > 1*/
	temp=CreateDMatrix(&gstack,numceps+1,1);
	temp1=CreateDMatrix(&gstack,1,(int)ceil((2*numChans-1)/2)+1);
	initDMatrix(0,1,numceps+1,temp);
	initDMatrix(1,2,2*numChans-1,temp1);
	temp2 = mulDMatrixRC(temp,temp1);
	temp2 = mulNumberwDMatrix(PI/(2*numChans),temp2);
	temp2 = cosDMatrix(temp2);
	dctDMatrix = mulNumberwDMatrix(sqrt(2.0/numChans),temp2);
	dctDMatrix = kron(eye(numDctDVectors,numDctDVectors),dctDMatrix);

	/*% inv DCT DMatrix
	  invDctDMatrix = inv(sqrt(2.0/numChans) * cos ( (pi/(2*numChans)) * (0:1:(numChans-1))' * (1:2:(2*numChans-1)) ));
	  invDctDMatrix = invDctDMatrix(:, 1:(numceps+1)); % cut off unnecessary cols
	  invDctDMatrix = kron(eye(numDctDVectors), invDctDMatrix); % for trajectory model if numDctDVectors > 1
	 */
	temp=CreateDMatrix(&gstack,numChans,1);
	initDMatrix(0,1,numChans,temp);
	temp2 = mulDMatrixRC(temp,temp1);
	temp2 = mulNumberwDMatrix(PI/(2*numChans),temp2);
	temp2 = cosDMatrix(temp2);
	invDctDMatrix = mulNumberwDMatrix(sqrt(2.0/numChans),temp2);
	DMatInvert(&gstack, invDctDMatrix,invDctDMatrix);
	invDctDMatrix= cutColDMatrix(invDctDMatrix, numceps+1);/*% cut off unnecessary cols*/
	invDctDMatrix = kron(eye(numDctDVectors,numDctDVectors), invDctDMatrix); /*% for trajectory model if numDctDVectors > 1
	/*% cepstral liftering coefficient
		cepDMatrix = diag(1+cepLifter*0.5*sin((0:numceps)'*pi/cepLifter));
		cepDMatrix = kron(eye(numDctDVectors), cepDMatrix); % for trajectory model if numDctDVectors > 1
	*/
	temp = CreateDMatrix(&gstack,numceps+1,1);
	initDMatrix(0,1,numceps+1,temp);
	temp = mulNumberwDMatrix(PI/cepLifter,temp);
	temp = mulNumberwDMatrix(cepLifter*0.5,sinDMatrix(temp));
	temp2 = sumNumberwDMatrix(1,temp);
	cepDMatrix = diag(temp2);
	cepDMatrix = kron(eye(numDctDVectors,numDctDVectors), cepDMatrix);/* % for trajectory model if numDctDVectors > 1
	/*% inv cepstral liftering coefficient
	  invCepDMatrix = diag(1./(1+cepLifter*0.5*sin((0:numceps)*pi/cepLifter)));
	  invCepDMatrix = kron(eye(numDctDVectors), invCepDMatrix);  % for trajectory model if numDctDVectors > 1
	*/
	temp = CreateDMatrix(&gstack,1,numceps+1);
	initDMatrix(0,1,numceps+1,temp);
	temp = mulNumberwDMatrix(PI/cepLifter,temp);
	temp = mulNumberwDMatrix(cepLifter*0.5,sinDMatrix(temp));
	temp2 = sumNumberwDMatrix(1,temp);
	temp2 = divNumberwDMatrix(1,temp2);
	invCepDMatrix = diag(temp2);
	invCepDMatrix = kron(eye(numDctDVectors,numDctDVectors), invCepDMatrix); /* % for trajectory model if numDctDVectors > 1*/

	dctMatrices[1]=dctDMatrix;
	dctMatrices[2]=cepDMatrix;
	dctMatrices[3]=invDctDMatrix;
	dctMatrices[4]=invCepDMatrix;

	return dctMatrices;

}
/* DLUDecompose: perform LU decomposition on Matrix a, the permutation
       of the rows is returned in perm and sign is returned as +/-1
       depending on whether there was an even/odd number of row
       interchanges */
Boolean DLUDecompose1(DMatrix a, int *perm, int *sign)
{
   int i,imax,j,k,n;
   double scale,sum,xx,yy;
   DVector vv,tmp;

   n = NumDRows(a); imax = 0;
   vv = CreateDVector(&gstack,n);
   *sign = 1;
   for (i=1; i<=n; i++) {
      scale = 0.0;
      for (j=1; j<=n; j++)
         if ((xx = fabs(a[i][j])) > scale )
            scale = xx;
      if (scale == 0.0) {
         HError(-1,"LUDecompose: Matrix is Singular");
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
   FreeDVector(&gstack,vv);
   return(TRUE);
}

/* DLinSolve: solve the set of linear equations Ax = b, returning
        the result x in  b */
void DLinSolve1(DMatrix a, int *perm, double *b)
{
	  int i,ii=0,ip,j,n;
	   double sum;
	   int r=0;
	   n=NumDRows(a);
	   for (i=1;i<=n;i++) {
	      ip=perm[i]; sum=b[ip]; b[ip]=b[i];
	      if (ii)
	         for (j=ii;j<=i-1;j++)
	         {
	        	 if(a[i][j]!=0 && b[j]!=0)
	        	 {

	        		 sum -=a[i][j]*b[j];

	        	 }

	         }

	      else
	         if (sum) ii=i;
	      b[i]=sum;
	   }
	   for (i=n; i>=1; i--) {
	      sum=b[i];
	      for (j=i+1; j<=n; j++)
	    	  if(a[i][j]!=0 && b[j]!=0)
	         sum -=a[i][j]*b[j];
	      b[i]=sum/a[i][i];
	   }
}

/* Inverting a double matrix */
void DMatrixInvert(DMatrix c)
{

	DMatrix a;
   double col[100];
   int sign;
   int n,i,j,perm[100];

   n=NumDRows(c);
   a=CreateDMatrix(&gstack,n,n);
   CopyDMatrix(c,a);           /* Make a copy of c */
   DLUDecompose1(a,perm,&sign);      /* Do LU Decomposition */
   for (j=1; j<=n; j++) {     /* Invert matrix */
	  for (i=1; i<=n; i++)
		 col[i]=0.0;
	  col[j]=1.0;
	  DLinSolve1(a,perm,col);
	  for (i=1; i<=n; i++)
		 c[i][j] = col[i];
   }

   FreeDMatrix(&gstack,a);


}
