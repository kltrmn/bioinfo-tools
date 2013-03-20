/*=================================================================
 *  motifind.cpp
 *
 *  ind = motifind(seq1,seq2,pct_ident)
 *
 *  returns indicies in seq1 where seq2 has >= pct_ident 
 *  percentage of characters in common excluding overlapping words
 * 
 *  
 *  Brian Kolterman 8/2012
 *=================================================================*/


#include <stdio.h>
#include <string.h> /* strlen */
#include "mex.h"


#define PID     prhs[2]
#define OUT     plhs[0]

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    char    *str1, *str2;
    int     lSt1, lSt2, iPos, iCh, iLast, nHits;
    double  score, pct_ident, *ind;
    int     *indTemp;
    
    
    // Check for correct number of arguments     
    
    if (nrhs != 3) 
    {
        mexErrMsgTxt("Usage: ind = motifind(seq1,seq2,pct_ident)\n");
    } 
    if (nlhs > 1)
    {
        mexErrMsgTxt("Usage: ind = motifind(seq1,seq2,pct_ident)\n");
    }
    
    // Check to be sure inputs are correct
    
    if (!(mxIsChar(prhs[0])) && !(mxIsChar(prhs[1])))
    {
        mexErrMsgTxt("seq1 and seq2 must be of type string.\n.");
    }
    
    str1=mxArrayToString(prhs[0]);
    str2=mxArrayToString(prhs[1]);
    
    lSt1 = (int)strlen(str1);
    lSt2 = (int)strlen(str2);
    
    if (lSt1 < lSt2)
    {
        mxFree(str1);
        mxFree(str2);    
        mexErrMsgTxt("Length of str1 must be longer than or equal to length of str2.\n");
    }
    
     if (mxGetM(PID) != 1 && mxGetN(PID) != 1)
    {
        mexErrMsgTxt("pct_ident must be a scalar.\n.");
    }
    
    pct_ident = mxGetScalar(PID);
    
    
    // Set up temproary storage for hit indicies
    
    indTemp = (int*)mxCalloc(lSt1,sizeof(int));
    
    iLast = (lSt1 - lSt2 + 1);
    
    nHits = 0;
    
    
    // Do the comparison
    
    for (iPos = 0; iPos < iLast; iPos++)
    {
        score = 0.0;
        
        for (iCh = 0; iCh < lSt2; iCh++)
        {
            if (str1[iCh+iPos] == str2[iCh])
            {
                score = score + 1.0;
            }
        }
        
        score = score/(double)lSt2;
        
        if (score >= pct_ident)
        {
            indTemp[nHits] = iPos+1;
            nHits++;
            iPos += (lSt2-1);
        }
    }
     
    // Remove extra zeros 
    
    OUT = mxCreateDoubleMatrix(1, nHits, mxREAL);
    ind = mxGetPr(OUT);
    
    for (iPos = 0; iPos < nHits; iPos++)
    {
        ind[iPos] = indTemp[iPos];
    }
    
    mxFree(str1);
    mxFree(str2);
    
    return;
}
