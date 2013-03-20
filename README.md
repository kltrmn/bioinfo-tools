Here are a couple of simple DNA motif searching tools written as mex files for use in MATLAB. The input and output formats are based on those found in the Bioinformatics toolbox for easier integration into existing MATLAB scripts.
Compile in MATLAB using: mex \<filename>

#### hamseqGen.c
  returns string containing all possible nucleotide ('A','T','C','G') words of a given length
  
#### motifcount.c 
  returns cell array containing the number of repeats found in a given DNA sequence for all possible words of a given 
  size. reverse compliments are counted together and overlaping words are counted as 1. 
  
#### motifind.cpp
  returns indicies in DNA sequence where a given motif has >= pct_ident
  percentage of characters in common excluding overlapping words. only forward strand is matched. 
  overlapping words are counted as 1.
  
#### motifind_revcomp.cpp 
  returns indicies in DNA sequence where a given input motif has >= pct_ident
  percentage of characters in common. reverse compliments are counted and 
  overlapping words are counted as 1.
  
#### motifind_revcomp_profile.cpp 
  same as above except input motif is represented as a position weight matrix of nucleotide frequencies.
  
#### subseqcount.c
  returns cell array containing all subsequences and # of repeats found in seq1 for each subsequence
  of length (motif_size) in seq2 having >= pct_ident percentage of characters in common.
