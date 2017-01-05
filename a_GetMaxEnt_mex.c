#include <math.h>
#include "mex.h"
#include "matrix.h"

#define LOWER pow(10,-15.0)
#define PI (3.141592653589793)

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *lw;				/* 1x3 input */
    double *mu, *var;		/* 3*Nx1 input matrix */
	size_t m;				/* matrix dimensions */
	size_t nrows;
    double *B;				/* output Nx1 matrix */

	mwSize i,j,k,n;
	int indi, indj;
    double varc, wa, wb, B1, B2, B3;	
	double nlpdf[3][3];
	double lzv[3];

    /* create a pointer to the real data in the input matrix  */
	lw = mxGetPr(prhs[0]);
    mu = mxGetPr(prhs[1]);
    var = mxGetPr(prhs[2]);

    /* get rows of the input matrix */
	m = mxGetM(prhs[1]);
	nrows = m/3;

	/* create the matrix */
    plhs[0] = mxCreateDoubleMatrix(nrows,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    B = mxGetPr(plhs[0]);

    /* call the computational routine */
	n=nrows;
	for (k=0; k<n; k++)
	{

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				indi = k+nrows*i;
				indj = k+nrows*j;
				varc = var[indi]+var[indj];
				nlpdf[i][j] = log(varc*2*PI)/2.0 + (mu[indi]-mu[indj])*(mu[indi]-mu[indj])/2.0/varc;
			}
		}

		B1=0;
		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				lzv[j] = lw[j] - nlpdf[i][j];
			}
			B1 -= exp(lw[i])*log( exp(lzv[0])+exp(lzv[1])+exp(lzv[2]) );
		}
    
		B[k]=B1;
	}

}

