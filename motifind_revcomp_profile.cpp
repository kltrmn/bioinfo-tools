/*=================================================================
 *  motifind_revcomp_profile.cpp
 *
 *  ind = motifind_revcomp_profile(seq1,motif_profile,pct_ident)
 *
 *  returns indicies in seq1 where motif_profile has >= pct_ident 
 *  percentage of characters in common counting reverse-compliments 
 *  and excluding overlaping words
 * 
 *  motif_profile is a 4 x N matrix of nucleotide counts with 
 *      N = motif length and nucleotides order A C G T  
 *
 *  Brian Kolterman 8/2011 : upgraded to profile inputs 8/2012
 *
 *=================================================================*/


#include <stdio.h>
#include <string.h> /* strlen */
#include "mex.h"


#define PID     prhs[2]
#define OUT     plhs[0]
#define MAX      30
#define MAXSEQ   10000


void RevComp(double *substr, double *substrR, mwSize smotif);
void seqToInt(char *seqstr, int *seq, mwSize seqlen);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    char    *str1;
    int     lSt1, iPos, iCh, iLast, nHits;
    double  *motif_profile, *motif_profileR, score, scoreR, pct_ident, *ind;
    int     *indTemp, i, *sequence;
    mwSize  lSt2;
    
    // Check for correct number of arguments     
    
    if (nrhs != 3) 
    {
        mexErrMsgTxt("Usage: ind = motifind_revcomp_profile(seq1,motif_profile,pct_ident)\n");
    } 
    if (nlhs > 1)
    {
        mexErrMsgTxt("Usage: ind = motifind_revcomp_profile(seq1,motif_profile,pct_ident)\n");
    }
    
    // Check to be sure inputs are correct
    
    if (!(mxIsChar(prhs[0])))
    {
        mexErrMsgTxt("seq1 must be of type string.\n.");
    }
    
    if (!(mxIsDouble(prhs[1])) || !(mxGetM(prhs[1]) == 4))
    {
        mexErrMsgTxt("motif_profile must be 4 X motif_length double matrix (order ACGT).\n.");
        
    }
    
    
    
    
    str1=mxArrayToString(prhs[0]);
    
    lSt1 = (int)strlen(str1);
    lSt2 = (mwSize)mxGetN(prhs[1]);
    
    if (lSt1 > MAXSEQ)
    {
        mexErrMsgTxt("seq1 must be <= 10000 bp.\n.");
    }
    
    if (lSt2 > MAX)
    {
        mexErrMsgTxt("motif_length must be <= 30.\n.");
    }
    
    
    
    
    // convert main nuc. sequence into an integer array
    
    sequence = (int*)mxCalloc(lSt1,sizeof(int));
    
    seqToInt(str1,sequence,lSt1);
   
    
    
    motif_profile = (double*)mxCalloc(4*lSt2,sizeof(double));
    motif_profileR = (double*)mxCalloc(4*lSt2,sizeof(double));
   
    motif_profile=mxGetPr(prhs[1]);
   
    
     
    // normalize profile
    
    for (iPos = 0; iPos < lSt2; iPos++) 
    {
        i = iPos*4;
        score = motif_profile[i] + motif_profile[i+1] + motif_profile[i+2] + motif_profile[i+3];
        
        motif_profile[i] /= score;
        motif_profile[i+1] /= score;
        motif_profile[i+2] /= score;
        motif_profile[i+3] /= score;
    }
    
   
    RevComp(motif_profile,motif_profileR,lSt2);
    
    
    if (lSt1 < lSt2)
    {
        mxFree(str1);   
        mexErrMsgTxt("Length of str1 must be longer than or equal to length of motif.\n");
    }
    
     if (mxGetM(PID) != 1 && mxGetN(PID) != 1)
    {
        mexErrMsgTxt("pct_ident must be a scalar 0< pct_ident <=1 .\n.");
    }
    
    pct_ident = mxGetScalar(PID);
    
    if ((pct_ident <= 0) || (pct_ident > 1))
    {
        mexErrMsgTxt("pct_ident must be a scalar 0< pct_ident <=1 .\n.");
    }
    
  
    
    
    
    // Set up temproary storage for hit indicies
    
    indTemp = (int*)mxCalloc(lSt1,sizeof(int));
    
    iLast = (lSt1 - lSt2 + 1);
    
    nHits = 0;
    
    
    // Do the comparison
    
    for (iPos = 0; iPos < iLast; iPos++)
    {
        score = 0.0;
        scoreR = 0.0;
        
        
        for (iCh = 0; iCh < lSt2; iCh++)
        {
            
            i = sequence[iCh+iPos] + iCh*4;
            
            score += motif_profile[i];
            scoreR += motif_profileR[i];
        
        }
        
        
        score = score/(double)lSt2;
        scoreR = scoreR/(double)lSt2;
        
        if (score >= pct_ident || scoreR >= pct_ident)
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
   
    return;
}


// void RevComp takes substr profile and returns reverse compliment 
// profile in substrR 

void RevComp(double *substr, double *substrR, mwSize smotif)
{
    
    int i,j,k;
    
    for (i = 0; i < smotif; i++)
    {
        j = i*4;
        k = (smotif-i-1)*4;
        substrR[k] = substr[j+3];
        substrR[k+1] = substr[j+2];
        substrR[k+2] = substr[j+1];
        substrR[k+3] = substr[j];
    }
}


void seqToInt(char *seqstr, int *seq, mwSize seqlen) {
    
    int i;
    
    for (i = 0; i < seqlen; i++) 
    {
        
        switch  (seqstr[i]) 
        {
            case 'A':
                seq[i] = 0;
                break;
                
            case 'C':
                seq[i] = 1;
                break;
                
            case 'G':
                seq[i] = 2;
                break;
                
            case 'T':
                seq[i] = 3;
                break;
        }
        
    }
    
}
