/*				UndoRLE_fast.c  by Anthony Giljum
 *   This code undoes the run-length encoding used by fSTO_fast_RLE.c.
 */

#include <math.h>
#include <stdlib.h>
#include "mex.h"


/*********MEX Interface********/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *u; /* Input array to be de-RLE'd. */
    double *count_a; /* For moving values between u and u_full_scanner */
    double *u_full_scanner; /* Output array pointer. */
    double *end; /* Final address in record_ds. */
    int *dsRat_sqrt; /* 2D-downsampling ratio */
    int *res; /* The total number of pixels to be output. */
    bool *record_ds; /* Record of where was downsampled. */
    bool *ind; /* For scanning over record_ds */
    int count_b; /* For spreading out value in LR pixel */
    int dsRat;
    int m,n;
    int rlePx; /* Number of pixels in RLE-vector. */
    
    /* Error checking */
    if(nlhs != 1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "One output required.");
    }
    if(nrhs != 4) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                      "Four inputs required.");
    }
    
    u = mxGetPr(prhs[0]);
    record_ds = mxGetPr(prhs[1]);
    dsRat_sqrt = mxGetPr(prhs[2]);
    res = mxGetPr(prhs[3]);
    
    dsRat = (*dsRat_sqrt)*(*dsRat_sqrt);
    
    if(!mxIsClass(prhs[0],"double")){
        mexErrMsgTxt("Data must be of type double!!!\n");
    }
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if(n>1 && m>1){
        mexErrMsgTxt("Only vectors allowed.  No Matrices.\n");
    }
    
    if(n==1)
        n=m;
    
    rlePx = mxGetM(prhs[1]);
    if(rlePx==1)
        rlePx = mxGetN(prhs[1]);
    
    end = record_ds+rlePx;
    
    /* Create output array */
    plhs[0] = mxCreateDoubleMatrix(*res,1,mxREAL);
    /* Get pointer to output array values */
    u_full_scanner = mxGetDoubles(plhs[0]);
    /* Undo RLE for Input Vector */
    count_a = u;
    for(ind=record_ds;ind<end;ind+=1){
        if(*ind){
            /* If we did downsample, spread that value over the entire low-
             * resolution pixel in the new array. */
            for(count_b=0;count_b<dsRat;count_b+=1){
                // *u_full_scanner = *count_a / ((*dsRat_sqrt));
                *u_full_scanner = *count_a;
                u_full_scanner+=1;
            }
            count_a+=1;
        } else{
            /* Otherwise, we only set one high-resolution pixel equal to
             * that value in the new array. */
            *u_full_scanner = *count_a;
            count_a+=1;
            u_full_scanner+=1;
        }
    }
    
    return;
}


    

