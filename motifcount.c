/*=================================================================
 *  motifcount.c
 *
 *  ind = motifcount(seq,motif_size)
 *
 *  returns cell array containing motif and # of repeats found in seq for each subsequence 
 *  of length (motif_size) including reverse compliment 
 *  overlaps not included
 *  
 *  
 *  Brian Kolterman 9/2011
 *=================================================================*/


#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

#define SEQ      prhs[0]
#define MS       prhs[1]
#define OUT      plhs[0]
#define MAX      13

const char bases[4] = "ATCG";

void GetMotif(const int iMot, const mwSize smotif, char *substr);
void GetIndex(const char *substr, int *iMot);
void RevComp(char *substr);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	char    *str, *substr;
	int     iPos, iCh, iSub, iMot1, iMot2, *iNext, *moCount, lSeq;
	mwSize *dims, ndim, smotif, nmotif, cmotif, nsubseq;

	ndim = 2;

	/* Check for correct number of arguments */     

	if (nrhs != 2) 
	{
		mexErrMsgTxt("Usage: ind = subseqcount(seq,motif_size)\n");
	} 
	if (nlhs > 1)
	{
		mexErrMsgTxt("Usage: ind = subseqcount(seq,motif_size)\n");
	}

	/* Check to be sure inputs are correct */

	if (!(mxIsChar(SEQ)))
	{
		mexErrMsgTxt("seq must be of type string.\n.");
	}

	str=mxArrayToString(SEQ);
    
	lSeq = (int)strlen(str);

    for (iPos = 0; iPos < lSeq; iPos++)
    {
        if (str[iPos] != 'A' && str[iPos] != 'G' && str[iPos] != 'C' && str[iPos] != 'T')
        {
            mexErrMsgTxt("invalid sequence.\n.");
        }
    }
    
	if (mxGetM(MS) != 1 && mxGetN(MS) != 1)
	{
		mexErrMsgTxt("motif_size must be a scalar.\n.");
	}

	smotif = (mwSize)mxGetScalar(MS);

	if (smotif > 13) 
	{
		mexErrMsgTxt("motif_size must <= 13.\n.");
	}

	nsubseq = lSeq - smotif + 1;
    
    nmotif = (mwSize)pow(4.0,(double)smotif);
     
	/* Set up temproary storage for motifs, indicies and counts */

	substr = (char*)mxCalloc(smotif,sizeof(char*));
    
    moCount =  (int*)mxCalloc(nmotif,sizeof(int*));
	iNext = (int*)mxCalloc(nmotif,sizeof(int*));
   
    /*initialize moCount to remove reverse comp in results*/
    
    for (iPos = 0; iPos < nmotif; iPos++)
    {
        GetMotif(iPos,smotif,substr);
        RevComp(substr);
        GetIndex(substr,&iMot1);
        if (iMot1 < iPos)
        {
            moCount[iPos] = -1;
        }
    }

	/* iterate through subsequences */

	for (iSub = 0; iSub < nsubseq; iSub++)
	{

		for (iCh = 0; iCh < smotif; iCh++)
		{
			substr[iCh] = str[iSub+iCh];
		}

		/* Get motif index and increment count (include reverse compliment) */
         
            
            GetIndex(substr,&iMot1);
            RevComp(substr);
            GetIndex(substr,&iMot2);
                        
            if (iMot2 < iMot1)
            {
                iMot1 = iMot2;
            }
            
            /* ignore overlaps */
            
            if (iSub < iNext[iMot1]) continue; 
            
            moCount[iMot1]++;
            iNext[iMot1] = iSub + smotif; 
        
    }
    
    
    /* create matlab cell array with motif counts */

    cmotif = 0;
        
    for (iPos = 0; iPos < nmotif; iPos++)
    {
        if (moCount[iPos] != -1) cmotif++;
    }
   
	dims = (mwSize*)mxCalloc(2,sizeof(mwSize*)); 
	dims[0] = cmotif;
	dims[1] = 2;

    OUT = mxCreateCellArray(ndim, dims);
    
    iCh = 0;
    for (iPos = 0; iPos < nmotif; iPos++)
    {
        
        if (moCount[iPos] != -1)
        {
            GetMotif(iPos,smotif,substr);
            mxSetCell(OUT,iCh,mxCreateString(substr));
            mxSetCell(OUT,iCh+cmotif,mxCreateDoubleScalar(moCount[iPos]));
            iCh++;
        }
    }
     
    
    
    
    mxFree(str);
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

void GetIndex(const char *substr, int *iMot)
{
    int i, num, smotif;
    i = 0;
    *iMot = 0;
    smotif = strlen(substr);
    
    while (substr[i] != 0)
    {
       
        switch (substr[i])
        {
    
        case 'A':
        num = 0;
        break;
        
        case 'T':
        num = 1;
        break;
       
        case 'C':
        num = 2;
        break;
        
        case 'G':
        num = 3;
        break;
        }    
        
        *iMot += num*(int)pow(4.0,(double)(smotif-i-1));
        
        i++;
    }
    
}

void RevComp(char *substr)
{
    char temp[MAX];
    int i, smotif;
    smotif = strlen(substr);
    
    for (i = 0; i < smotif; i++)
    {
        temp[i] = substr[i];
    }
    
    i = 0;
    
    while (smotif >= 0)
    {
        smotif--;
        
        switch  (temp[smotif])
        {
            case 'A':
            substr[i] = 'T';
            break;
            
            case 'T':
            substr[i] = 'A';
            break;

            case 'C':
            substr[i] = 'G';
            break;
            
            case 'G':
            substr[i] = 'C';
            break;
        }
            
        i++;
    }
    
}
        
   

 
