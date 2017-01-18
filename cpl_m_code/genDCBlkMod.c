
#include <math.h> /* Needed for the ceil() prototype */
#include "mex.h"


/*#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif*/

/* Generates block model 
   first input is "c", second input is "P", third input "nzmax" */

void mexFunction(
		 int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]
		 )
{
    /* Declare variable */
    mwSize n, K;
    mwSize nzmax;
    mwIndex *Air,*Ajc,j,k;
    int cmplx, HadToResize; /*,isfull; */
    double *A, *P, *c, *theta;
    double test;
    
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 4) {
        mexErrMsgTxt("Four input arguments required.");
    } 
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0]))) {
        mexErrMsgTxt("Input argument must be of type double.");
    }	
    
    if (mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgTxt("'c' must be  2-dimensional (a vector) \n");
    }

    if (mxGetNumberOfDimensions(prhs[1]) != 2) {
        mexErrMsgTxt("'P' must be  2-dimensional (a matrix) \n");
    }

    /* Get the size and pointers to input data */
    
    /* first input is "c", second input is "P", third input "nzmax", fourth input "theta" */
    K = mxGetM(prhs[1]);
    n = mxGetM(prhs[0]);
    
    c = mxGetPr(prhs[0]);
    P = mxGetPr(prhs[1]);
    theta = mxGetPr(prhs[3]);
    
    
    /* Allocate space for the sparse adjacency matrix 
       nzmax is our initial estimate of the number of nonzeros */
    nzmax = (mwSize) mxGetPr(prhs[2])[0];  
    
    cmplx = 0;
    plhs[0] = mxCreateSparse(n,n,nzmax,cmplx);
    
    A  = mxGetPr(plhs[0]);
    
    Air = mxGetIr(plhs[0]);
    Ajc = mxGetJc(plhs[0]);
    
    
    k = 0; 
    
    HadToResize = 0;
    
    for (j = 0; j < n; j++) {
      
      mwSize i;
      Ajc[j] = k;
      
      for (i = 0; i < j; i++) {

          /* Check to see if non-zero element will fit in 
           * allocated output array.  If not, increase percent_sparse 
           * by 10%, recalculate nzmax, and augment the sparse array
          */
          if (k>=nzmax) {
            mwSize oldnzmax = nzmax;
            nzmax = (mwSize) ceil(((double) nzmax)*1.1);
            
            HadToResize++;
            
            /* make sure nzmax increases atleast by 1 */
            if (oldnzmax == nzmax) 
              nzmax++;

            mxSetNzmax(plhs[0], nzmax); 
            mxSetPr(plhs[0], mxRealloc(A, nzmax*sizeof(double)));
           
            mxSetIr(plhs[0], mxRealloc(Air, nzmax*sizeof(mwIndex)));

            A  = mxGetPr(plhs[0]);
            Air = mxGetIr(plhs[0]);
          } 

            
          /* The acutal work of generating the model is the following if statement */
          test =  theta[i] * theta[j] * P[ ((mwIndex) c[j]-1)*K + ((mwIndex) c[i])-1 ];
          if ( ((double) rand())/RAND_MAX < test) { 

              A[k] = 1; 
              Air[k] = i;
              k++;
          }
      }
    }
    Ajc[n] = k;
    
    if (HadToResize > 0) {
        mexPrintf("Had to reallocate %d times.\n",HadToResize);
    }
}
 
