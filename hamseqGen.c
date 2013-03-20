/*=================================================================
 *  motifcount.c
 *
 *  str = hamseqGen(motif_size)
 *
 *  returns string containing all possible motifs 
 *  of length (motif_size) 
 *  
 *  Brian Kolterman 9/2012
 *=================================================================*/


#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"


#define MS       prhs[0]
#define STR      plhs[0]
#define MAX      13

const char bases[4] = "ATCG";

void GetMotif(const int iMot, const mwSize smotif, char *substr);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	char    *motifs, *substr;
	int     iPos, iCh, count;
	mwSize *dims, ndim, smotif, nmotif, cmotif, nsubseq;

	ndim = 2;

	/* Check for correct number of arguments */     

	if (nrhs != 1) 
	{
		mexErrMsgTxt("Usage: ind = hamseqGen(motif_size)\n");
	} 
	if (nlhs > 1)
	{
		mexErrMsgTxt("Usage: ind = hamseqGen(motif_size)\n");
	}

	/* Check to be sure inputs are correct */

	
	if (mxGetM(MS) != 1 && mxGetN(MS) != 1)
	{
		mexErrMsgTxt("motif_size must be a scalar.\n.");
	}

	smotif = (mwSize)mxGetScalar(MS);

	if (smotif > 13) 
	{
		mexErrMsgTxt("motif_size must <= 13.\n.");
	}

    nmotif = (mwSize)pow(4.0,(double)smotif);
     
	/* Set up temproary storage for motifs, indicies and counts */

    substr = (char*)mxCalloc(smotif,sizeof(char*));
	motifs = (char*)mxCalloc(smotif*nmotif,sizeof(char*));
    
    
    /* create matlab array with all motifs */

        
    count = 0;
    for (iPos = 0; iPos < nmotif; iPos++)
    {
        GetMotif(iPos,smotif,substr);
    
        for (iCh = 0; iCh < smotif; iCh++)
        {
            motifs[count] = substr[iCh];
            count++;
        }
    }
   
	 
    
    STR = mxCreateString(motifs);
    
 
    mxFree(substr);
    
    return;
}

void GetMotif(const int iMot, const mwSize smotif, char *substr)
{
    int i, p, rem;
 
    rem = iMot;
    
    for (i = 0; i < smotif; i++)
    {
        p = (int)floor((double)rem/pow(4.0,(double)(smotif-i-1)));
        rem = rem % (int)pow(4.0,(double)(smotif-i-1));
        substr[i] = bases[p];
    }
    
}


