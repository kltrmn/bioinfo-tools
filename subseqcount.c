/*=================================================================
 *  subseqcount.c
 *
 *  ind = subseqcount(seq1,seq2,motif_size,pct_ident)
 *
 *  returns cell array containing subseq and # of repeats found in seq1 for each subsequence 
 *  of length (motif_size) in 
 *  seq2 having >= pct_ident 
 *  percentage of characters in common
 * 
 *  
 *  Brian Kolterman 8/2012
 *=================================================================*/


#include <stdio.h>
#include <string.h> 
#include "mex.h"
#include "matrix.h"

#define MS      prhs[2]
#define PID     prhs[3]
#define OUT     plhs[0]

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	char    *str1, *str2, *substr;
	int     lSt1, lSt2, iPos, iCh, iLast, iSub;
	double  score, pct_ident, nHits;
	mwSize *dims, ndim, smotif, nmotif;

	ndim = 2;

	/* Check for correct number of arguments */     

	if (nrhs != 4) 
	{
		mexErrMsgTxt("Usage: ind = subseqcount(seq1,seq2,motif_size,pct_ident)\n");
	} 
	if (nlhs > 1)
	{
		mexErrMsgTxt("Usage: ind = subseqcount(seq1,seq2,motif_size,pct_ident)\n");
	}

	/* Check to be sure inputs are correct */

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


	if (mxGetM(MS) != 1 && mxGetN(MS) != 1)
	{
		mexErrMsgTxt("motif_size must be a scalar.\n.");
	}

	smotif = (mwSize)mxGetScalar(MS);

	if (smotif > 15) 
	{
		mexErrMsgTxt("motif_size must <= 15.\n.");
	}

	nmotif = lSt2 - smotif + 1;

	/* Set up temproary storage and output cell array */

	substr = (char*)mxCalloc(smotif,sizeof(char*));

	iLast = (lSt1 - smotif + 1);

	dims = (mwSize*)mxCalloc(2,sizeof(mwSize*)); 
	dims[0] = nmotif;
	dims[1] = 2;


	OUT = mxCreateCellArray(ndim, dims);
	/*fout = (mxArray*)mxCreateCellMatrix(100, 2);*/


	/* Start iterating through subsequences */

	for (iSub = 0; iSub < nmotif; iSub++)
	{

		for (iPos = 0; iPos < smotif; iPos++)
		{
			substr[iPos] = str2[iSub+iPos];
		}

		nHits = 0.0;


        

	/* Do the comparison */
        
        for (iPos = 0; iPos < iLast; iPos++) 
        {
            score = 0.0;
            
             for (iCh = 0; iCh < smotif; iCh++) 
             {
                        
                if (str1[iCh+iPos] == substr[iCh]) 
                {
               
                    score += 1.0;
                
                }
            }
            
            score = score/(double)smotif;
            
            if (score >= pct_ident) 
            {
                nHits += 1.0;
                iPos += (smotif - 1); /* avoid overlaps */
            }
        }
        
        mxSetCell(OUT,iSub,mxCreateString(substr));
        
        mxSetCell(OUT,iSub+nmotif,mxCreateDoubleScalar(nHits));
       
    }
     
    
    
    
    mxFree(str1);
    mxFree(str2);
    mxFree(substr);
    
    return;
}
