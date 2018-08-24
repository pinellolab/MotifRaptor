# Fast matching of Transcription Factors & SNPs
### Author: Paolo Ferragina 
### August 2018 @ PinelloLab

import cython
cimport cython
import numpy as np
cimport numpy as np
import pandas as pd
import time

ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.float_t DTYPE_FLOAT



################################################################################################################
# Parameters:
#
# int    n : length of the input sequence
# int    num_blocks : number of blocks forming the input sequence
# int    block_size : size of the blocks forming the input sequence (61)
# char[] text : array of n characters denoting the sequence to be processed
# int[]  sa : array of n integers, starting positions of text suffixes text[sa[i]: ] sorted lexicographically 
# int[]  lcp : array of n integers, denoting the lcp between suffixes T[sa[i-1]: ] and T[sa[i]: ]
# char[] motif_file_name : name of the file where the motif is stored
#################################################################################################################
@cython.boundscheck(False)
@cython.nonecheck(False)
def motif_matching(int n, int num_blocks, int block_size, np.ndarray[np.uint_t, ndim=1] lcp, np.ndarray[np.uint_t, ndim=1] sa, str text, np.ndarray[np.float_t, ndim=2] motif_pfm_matrix):

    # We keep an array of scores, one per prefix of the motif
    # 
    # At every step we compare the MOTIF against the substring s = T[ sa[i] : sa[i] + motif_length]
    #   but the comparison takes advantage of the comparison with the "previous" substring
    #   s' = T[ sa[i-1] : sa[i-1] + motif_length] in the sorted lexicographic order given by 
    #   the suffix array SA and of the scores stored in SCORE_ARRAY for that substring s'.
    #   In fact, the code skips the recomputation of the scores for the prefix of s shared with s', that is,
    #   it skips the first lcp[i-1] bases, because they are shared by the two substrings and thus their scores
    #   are the same.


    cdef int motif_length
    #cdef np.ndarray[DTYPE_FLOAT, ndim=2] motif_pfm_matrix
    
    # Loading the motif to match
    #print("Loading the motif : " + motif_file_name)
    #pseudo_count=0.001
    #motif_pfm_matrix=generate_pfm_matrix(motif_file_name,pseudo_count)
    motif_length = motif_pfm_matrix.shape[1]
    #print(str(motif_length) + " bases")
    #print("done!")


    #print("Starting scoring computation... ")
    start_time = time.time()

    # Create and initializes the two matrices that contain the scores
    #   ref_scores = contains the scores for the reference sequence
    #   alt_scores = contains the scores for the alternative sequence
    #
    # Keep attention that actually an offset larger than 30, we count from 0, gives 
    #   a not meaningful score, so the corresponding entry is zeroed to do not account for it.
    # We need anyway to compute the scores because of the scanning procedure operated over 
    #   the suffix array

    cdef np.ndarray[DTYPE_FLOAT, ndim=2] ref_scores
    cdef np.ndarray[DTYPE_FLOAT, ndim=2] alt_scores

    ref_scores = np.zeros((num_blocks,block_size), dtype = np.float)
    alt_scores = np.zeros((num_blocks,block_size), dtype = np.float)


    # The first SCORE_ARRAY is computed from scratch
    cdef np.ndarray[DTYPE_FLOAT, ndim=1] score_array
    score_array = np.zeros((motif_length), dtype = np.float)

    # Additional infos for debugging
    # print("\n   index in SA = ", 0) 
    # print("   substring = ", text[sa[0]:sa[0]+motif_length].upper()) # very slow but for check
    # print("   lcp with previous substring in the lexicographic order = ", 0)
    # print("   entire score array = ", score_array)
   
    
    # The other SCORE_ARRAYs are computed by "difference"
    #   by recomputing only starting from the lcp[i-1]-th base (if any)

    cdef int i,j,start_range
    
    
    for i in range(0,n):

        # The first suffix in the sorted order has not predecessor and thus
        #  the computation of the score_array has to start from the beginning
        if (i == 0):
            start_range = 0
        else:
            start_range = min(lcp[i-1], motif_length)
    
        # Notice that the computation occurs only when lcp[i-1] < motif_length
        #  thus, only when there are characters changed in the matched part of the motif
        for j in range(start_range, min(motif_length, n - sa[i])):

            char = text[sa[i]+j].upper()
            if (char not in ['A','C','G','T', 'N']):  # unknown char in the sequence
                raise Exception("\n\nChar not allowed in sequence file\n")

            if (char == 'N'): # wildcard char
                score_array[j] = 0
            elif (char == 'A'): 
                score_array[j] = motif_pfm_matrix[0,j]
            elif (char == 'C'): 
                score_array[j] = motif_pfm_matrix[1,j]
            elif (char == 'G'): 
                score_array[j] = motif_pfm_matrix[2,j]
            else: 
                score_array[j] = motif_pfm_matrix[3,j]

            if (j>0):
                score_array[j] = score_array[j] + score_array[j-1] 
   
        # 0 for debugging, instead of -1 (which is impossible, of course)
        #if ((i>0) and (i % 10000000 == 0)):
        #    t = time.time() - start_time
        #    remaining_time = t * ((float(n) / float(i)) - 1)
        #    s0 = time.strftime("%H:%M:%S", time.gmtime(t)) 
        #    s1 = time.strftime("%H:%M:%S", time.gmtime(remaining_time))
        #    print(" parsed ", end='')
        #    print(f"{i:,d}: ", end='')
        #    print(" elapsed time %s, estimated completion time %s. " % (s0, s1))

           
            
        # The following IF excludes the last motif_length positions in the file
        if (sa[i] <= n - motif_length):

            # For debugging purposes:
            #
            # if (lcp[i-1] < motif_length):
            #    print("\n\n position in sa = ", i, "position of the match = ", sa[i], ", score = ", score_array[motif_length-1])
            # else:
            #    print("\n\n [free] position in sa = ", i, "position of the match = ", sa[i], ", score = ", score_array[motif_length-1])
            #
            # Additional infos for debugging
            # print("\n   index in SA = ", i) 
            # print("   substring = ", text[sa[i]:sa[i]+motif_length].upper())  # very slow, only for check
            # print("   lcp with previous substring in the lexicographic order = ",lcp[i-1])
            # print("   entire score array = ", score_array) # very slow, only for check

            # Store the score in the appropriate place, either ref_scores or alt_score, depending on the
            #   block where the suffix text[sa[i]: ] starts.
            snp_block = ( sa[i] // block_size )
            snp = ( sa[i] // (2*block_size) )            
            block_offset = sa[i] % block_size
            if (block_offset < 31): # We initialize only whether it is before the SNP position (30)
                if (snp_block % 2 == 0):  # it is a SNP in the reference sequence
                    ref_scores[snp][block_offset] = score_array[motif_length-1]
                else:                 # it is a SNP in the alternative sequence
                    alt_scores[snp][block_offset] = score_array[motif_length-1]

        

    #print(" done in " +str( (time.time() - start_time) )+" seconds " + "for motif: " + motif_file_name)

    return motif_length, ref_scores, alt_scores
