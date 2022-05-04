/*				fSTO_fast_RLE.c  by Anthony Giljum
 *   This code implements the fast foveated Sum-To-One Transform
 * 
 *   Based on STO_fast.c by Tom Goldstein
 *   Updated for the Fast Foveated STOne Transform by Anthony Giljum
 *
 *   The idea behind the algorithm is to build the full transform at the
 *   low resolution, and then go over the high-resolution bits separately.
 *   It is important to keep in mind the correct scaling when doing the 
 *   high resolution regions of the transform.
 *   Once the full transform is performed, foveation is implemented by 
 *   downsampling the low-resolution region, representing each pixel as a
 *   single value independent of effective pixel size.
 *
 */

#include <math.h>
#include <stdlib.h>
#include "mex.h"


/*********MEX Interface********/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* u;
    double *ind; /* For looping over full transform */
    double* HR_ind; /* For recording high-res values into RLE vector. */
    double* u_rle; /* Run-Length Encoded Foveated Vector */
    double* end;
    double *first, *second, *third, *fourth;
    double *start, *stop;
    double *sumInd; /* For looping over LR px */
    double sumVal; /* Sum of HR px within LR px. */
    int *highres, *lowres; /* High- and low-resolutions */
    int *LRreg; /* Low-resolution region indices */
    int dsRat; /* Downsample ratio */
    int dsRat_sqrt; /* Square root of downsample ratio */
    int paramLen; /* Length of the foveation-parameters vector */
    int gap;
    int m,n;
    double a,b,c,d;
    bool DoFoveate; /* Perform the foveation. */
    
    u = mxGetPr(prhs[0]);
    
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
    
    end = u+n;
    
    /* This funky block checks that the length is a power of 4 */
    if ( n & (n - 1)  
     || (n & 0xAAAAAAAB ) ){
        mexErrMsgTxt("Vector length must be power of 4.\n");
    }
    
    
    /* Get length of foveation-parameters vector */
    if(mxGetM(prhs[1])==1){
        paramLen = mxGetN(prhs[1]);
    } else {
        paramLen = mxGetM(prhs[1]);
    }
    /* Check whether we are foveating. */
    if(paramLen>2){ /* Are foveating */
        highres = mxGetPr(prhs[1]); /* Added */
        lowres = highres+1; /* Added */
        LRreg = highres+2; /* Added */
        DoFoveate = 1;
        dsRat_sqrt = *highres / *lowres; /* Compute downsampling ratio */
        dsRat = dsRat_sqrt*dsRat_sqrt;
    } else if(paramLen==2){ /* Not foveating */
        highres = mxGetPr(prhs[1]); /* Added */
        lowres = highres+1; /* Added */
        DoFoveate = 0;
    } else if(paramLen<2){ /* Error */
        mexErrMsgTxt("Foveation parameters not in correct form.\n");
    }
    
    
    /*Do the First level transform.  This is done in its own loop for efficiency*/
    first = u;
    second = u+1;
    third = u+2;
    fourth = u+3;
    for(;first<end; first+=4, second+=4, third+=4, fourth+=4){
        a = (*first)+(*second)+(*third)-(*fourth);
        b = (*first)+(*second)-(*third)+(*fourth);
        c = (*first)-(*second)+(*third)+(*fourth);
        d = -(*first)+(*second)+(*third)+(*fourth);
        *first = a*.5;
        *second = b*.5;
        *third = c*.5;
        *fourth = d*.5;
    }
    
    /* Do every other level of the transform */
    for(gap = 4; gap<n ;gap*=4){
        first = u;
        second = u+gap;
        third = u+2*gap;
        fourth = u+3*gap;
        for(start = u, stop=u+gap; start<end; start+=gap*4, stop+=gap*4){
            for(;first<stop; first++, second++, third++, fourth++){
                a = (*first)+(*second)+(*third)-(*fourth);
                b = (*first)+(*second)-(*third)+(*fourth);
                c = (*first)-(*second)+(*third)+(*fourth);
                d = -(*first)+(*second)+(*third)+(*fourth);
                *first = a*.5;
                *second = b*.5;
                *third = c*.5;
                *fourth = d*.5;
            }
            first+=3*gap;
            second+=3*gap;
            third+=3*gap;
            fourth+=3*gap;
        }
    }
    
    /* Foveation Section
     * By this point, we have a full-resolution transform
     * of the input. So, since our transform has this nice
     * sum-to-one property, we downsample by summing together
     * the high-res pixels that make up the low-res region.
     */
    if(DoFoveate){
        u_rle = u;
        /* Loop over every pixel in the vector */
        for(ind = u; ind<end ;ind+=dsRat){
            if(ind==u+*LRreg-1){ /* If we are not in the foveated region */
                /* Sum together all high res px contained by the
                 * low res px. */
                sumVal = 0;
                for(sumInd = ind; sumInd<ind+dsRat; sumInd+=1){
                    sumVal = sumVal + *sumInd;
                }
                /* Set values within low-res pixel */
                // *u_rle = sumVal/(dsRat_sqrt);
                *u_rle = sumVal;
                u_rle+=1;
                /* Move to next LR pixel */
                LRreg += dsRat;
            } else{ /* If in HR-region, record values for RLE vector */
                for (HR_ind = ind; HR_ind<ind+dsRat; HR_ind+=1){
                    *u_rle = *HR_ind;
                    u_rle+=1;
                }
            }
        }
    }
    return;
}


    

