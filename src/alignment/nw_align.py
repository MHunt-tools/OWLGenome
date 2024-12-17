import numpy as np
from numba import njit

@njit(cache=True)
def banded_nw_align_numba(seq1_array, seq2_array, n, m, half_bandwidth, gap_penalty, match_score, mismatch_score):
    """
    Performs a banded needleman-wunsch (nw) global alignment using numpy arrays and numba enabled just in time machine code compilation
    banded_nw_align_numba constructs the dynamic programming (dp) table as a numpy array for use in a seperate traceback function

    INPUTS:
        seq1_array - numpy array of DNA characters converted into ascii values
        seq2_array - numpy array of DNA characters converted into ascii values
        n - integer length of seq1_array
        m - integer length of seq2_array
        half_bandwidth - integer distance from central band of dynamic programming table to edge of band, whole bandwidth is 2 * half_bandwidth + 1
        gap_penalty - integer gap penalty value for nw dp table filling
        match_score - integer match score value for nw dp table filling
        mismatch_score - integer mismatch_score value for nw dp table filling

    OUTPUTS:
        max_val_dp - numpy 2D array that is the filled dp table for the nw global alignment between seq1_array, seq2_array 
        traceback_dp -numpy 2D array that contains traceback pointers for every cell in max_val_dp for use in a seperate traceback function
        max_score - integer value representing the maximum alignment score found for the nw global alignment between seq1_array, seq2_array
    """
    # Initialize dp tables as 2D numpy arrays of dimensions n + 1 by 2 * half_bandwidth + 1
    max_val_dp = np.zeros((n + 1, 2 * half_bandwidth + 1), dtype=np.int32) # Set these values to be int32 to account for a safe range of integers
    traceback_dp = np.zeros((n + 1, 2 * half_bandwidth + 1), dtype=np.int8) # Set these values to be int8 as we know exactly what the pointer integers can be

    # Initialize pointer integers for the traceback dp table
    END = 0
    MATCH_MISMATCH = 1
    GAP_IN_SEQ1 = 2  
    GAP_IN_SEQ2 = 3  

    '''
    We define j <= m to be the column index for a non-banded 2D array
    j_band is the column index to access the band within a banded 2D array
    This is necessary as in a banded array, the rows are of length 2 * half_bandwidth + 1 and not of length m
    The filled band begins in the first row (i = 0) with the middle cell of the band over the j = 0 cell
    To convert j to j_band, we use j_band = half_bandwidth - i + j
    '''

    # Fill in the base case values for the first column in the dp table based off nw rules
    for i in range(1, n + 1):
        j_band = half_bandwidth - i + 0  # Convert j to j_band for the column indexed by j = 0
        if 0 <= j_band <= 2 * half_bandwidth: # Restrict the bound of j_band to not go out of bounds of the band
            max_val_dp[i, j_band] = gap_penalty * i # Fill in the first column according to the gap_penalty * i base case rule of nw
            traceback_dp[i, j_band] = GAP_IN_SEQ2  # Fill traceback array with the pointer for GAP_IN_SEQ2 which corresponds to i = (1:n), j = 0
        else:
            break # Break loop when bandwidth bound broken

    # Fill in the base case values for the first row in the dp table based off nw rules
    for j in range(1, m + 1):
        j_band = half_bandwidth - 0 + j # Convert j to j_band for the column indexed by i = 0
        if 0 <= j_band <= 2 * half_bandwidth:  # Restrict the bound of j_band to not go out of bounds of the band
            max_val_dp[0, j_band] = gap_penalty * j  # Fill in the first row according to the gap_penalty * i base case rule of nw
            traceback_dp[0, j_band] = GAP_IN_SEQ1  # Fill traceback array with the pointer for GAP_IN_SEQ1 which corresponds to i = 0, j = (1:m)
        else:
            break # Break loop when bandwidth bound broken

    # Fill in i, j = 0, 0 value and traceback pointer
    max_val_dp[0, half_bandwidth] = 0
    traceback_dp[0, half_bandwidth] = END
    max_score = 0

    # Fill the DP table
    for i in range(1, n + 1):
        start_j = max(1, i - half_bandwidth) # Get j index for start of band
        end_j = min(m, i + half_bandwidth) # Get j index for end of band

        for j in range(start_j, end_j + 1):
            j_band = half_bandwidth - i + j # Convert j indices into j_band indices

            # Calculate NW score for match or mismtach
            if 0 <= j_band <= 2 * half_bandwidth:
                if seq1_array[i - 1] == seq2_array[j - 1]: # match
                    score = match_score
                else: # mismatch
                    score = mismatch_score
                match_mismatch = max_val_dp[i - 1, j_band] + score # NW score
            else:
                match_mismatch = float('-inf') # Set scores to -inf if an out of bounds j_band is accessed

            # Calculate NW score for a gap in seq1
            if 0 <= j_band - 1 <= 2 * half_bandwidth:
                seq1_gap = max_val_dp[i, j_band - 1] + gap_penalty # NW score
            else:
                seq1_gap = float('-inf') # Set scores to -inf if an out of bounds j_band is accessed

            # Calculate NW score for a gap in seq2
            if 0 <= j_band + 1 <= 2 * half_bandwidth:
                seq2_gap = max_val_dp[i - 1, j_band + 1] + gap_penalty # NW score, moving up in i results in left shift of j_band value, so to stay at the same j we have to access j_band + 1
            else:
                seq2_gap = float('-inf') # Set scores to -inf if an out of bounds j_band is accessed

            # Find which movement gives the maximum NW score and set that as the value for the current cell of the NW dp table
            max_val = max(match_mismatch, seq1_gap, seq2_gap)
            max_val_dp[i, j_band] = max_val

            # Set the max_score to be the bottom right cell in the NW dp table
            if i == n and j_band == end_j:
                max_score = max_val_dp[i, j_band]


            # Set traceback pointers for the traceback table based off which move resulted in the maximum NW score
            if max_val == match_mismatch:
                traceback_dp[i, j_band] = MATCH_MISMATCH
            elif max_val == seq1_gap:
                traceback_dp[i, j_band] = GAP_IN_SEQ1
            elif max_val == seq2_gap:
                traceback_dp[i, j_band] = GAP_IN_SEQ2

    return max_val_dp, traceback_dp, max_score # Return the NW dp table, the traceback table, and the max score

def banded_nw_traceback(seq1, seq2, traceback_dp, n, m, half_bandwidth):
    """
    Using a pre-filled banded needleman-wunsch (nw) dynamic programming (dp) traceback table, construct the global alignment between seq1 and seq2 as alignment strings

    INPUTS:
        seq1 - string of DNA to be aligned
        seq2 - string of DNA to be aligned
        traceback_dp - numpy 2D array that contains traceback pointers for reconstruction of the global alignment between seq1 and seq2
        n - integer length of seq1
        m - integer length of seq2
        half_bandwidth - integer distance from central band of dynamic programming table to edge of band, whole bandwidth is 2 * half_bandwidth + 1

    OUTPUTS:
        aligned_seq1 - string of DNA for the seq1 global alignment
        alinged_seq2 - string of DNA for the seq2 global alignment
    """
  
    # Initialize pointer integers for the traceback dp table
    MATCH_MISMATCH = 1
    GAP_IN_SEQ1 = 2 
    GAP_IN_SEQ2 = 3  

    # Get indices for the bottom right cell in the traceback table
    i = n
    j_band = m - i + half_bandwidth

    # Initialize arrays to hold the alignment string characters
    aligned_seq1 = []
    aligned_seq2 = []

    # Iterate through the traceback until reaching the 0,0 in the traceback table
    # Alignment will proceed until 0,0 as no move that goes out of bounds was allowed when constructing the NW DP table
    while i > 0 and 0 <= j_band < 2 * half_bandwidth:
        traceback_code = traceback_dp[i, j_band] # Extract the traceback code from the array

        # Convert j_band to j for accessing the sequence strings
        j = j_band - half_bandwidth + i

        # Based off of the alignmend code, append the DNA string characters, or a gap, to the alignment string character arrays
        if traceback_code == MATCH_MISMATCH:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1 # A diagonal move preserves j_band, so only decrement i
        elif traceback_code == GAP_IN_SEQ1:
            aligned_seq1.append('_')
            aligned_seq2.append(seq2[j - 1])
            j_band -= 1  # Move left in the band
        elif traceback_code == GAP_IN_SEQ2:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('_')
            i -= 1
            j_band += 1  # Moving up in i shifts j_band to the left, so to preserve the j position we add 1 to j_band
        else:
            break

    # Reverse the characters array and join into a string to construct the final alignment strings
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return aligned_seq1, aligned_seq2 # Return the alignment strings


def nw_align(seq1, seq2, gap_penalty, match_score, mismatch_score):
    '''
    Performs a needleman-wunsch (nw) dynamic programming (dp) based global alignment and uses traceback to construct the alignment strings

    INPUTS:
        seq1 - string of DNA to be aligned
        seq2 - string of DNA to be aligned
        gap_penalty - integer gap penalty value for nw dp table filling
        match_score - integer match score value for nw dp table filling
        mismatch_score - integer mismatch_score value for nw dp table filling

    OUTPUTS:
        seq1_align - string of DNA for the seq1 global alignment
        seq2_align - string of DNA for the seq2 global alignment
    '''

  
    # DP tables are initialized to hold NW score and traceback information
    traceback_dp = [["gap"] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] # This table will keep backtracking pointers for reconstructing the sequence
    max_val_dp = [[gap_penalty] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] # This table holds the NW score values
    

    
    def score(a,b):
        '''
        Simple scoring function for character match or mismatches

        INPUTS:
            a - string character to be compared
            b - string character to be comapred

        OUTPUTS:
            match_score - integer score for if the character are the same, fed from nw_align
            mismatch_score - integer score for if the characters are different, fed from nw_align
        '''

        if a == b:
            return match_score

        else:
            return mismatch_score
        

    def traceback(traceback_dp, i, j): 
        """
        Using a banded needleman-wunsch (nw) dynamic programming (dp) traceback table, construct the global alignment between seq1 and seq2 as alignment strings

        INPUTS:
            traceback_dp - 2D array that contains traceback pointers for reconstruction of the global alignment between seq1 and seq2
            i - seq1 index
            j - seq2 index

        OUTPUTS:
            seq1_align - string of DNA for the seq1 global alignment
            seq2_align - string of DNA for the seq2 global alignment
        """

        seq_1_align = '' # Initialize strings
        seq_2_align = ''
        while i > 0 or j > 0: # Set up while loop break condition to end when eithe index is zero
            if traceback_dp[i][j]=="seq1_gap": # if a gap in seq1 is found, then append a gap to aligned seq1 and append the curret seq2 char to aligned seq2
                seq_1_align = "_" + seq_1_align 
                seq_2_align = seq2[j-1] + seq_2_align
                j -= 1 # Move to j-1
                continue
            elif traceback_dp[i][j]=="match_mismatch": # if a match/mismatch is found, append the current char for seq1 and seq2 to their aligned sequences
                seq_1_align = seq1[i-1] + seq_1_align
                seq_2_align = seq2[j-1] + seq_2_align
                i -= 1 # Move back diagonally
                j -= 1
                continue
            elif traceback_dp[i][j]=="seq2_gap": # if a gap in seq2 is found, then append a gap to aligned seq 2 and append the curret seq1 char to aligned seq1
                seq_1_align = seq1[i-1] + seq_1_align
                seq_2_align = "_" + seq_2_align 
                i -= 1 # Move to i - 1

        return seq_1_align, seq_2_align # Return aligned sequences
    
    
    # Dynamic table creation, iterate through both strings and fill in DP table
    for i in range(0,len(seq1)+1): 
        for j in range(0,len(seq2)+1):
            # Calculate values according to NW algorithm rules
            match_mismatch = max_val_dp[i-1][j-1] + score(seq1[i-1],seq2[j-1])
            seq1_gap = max_val_dp[i][j-1] + gap_penalty
            seq2_gap = max_val_dp[i-1][j] + gap_penalty
            max_val = max(match_mismatch,seq1_gap,seq2_gap)
            
            # Set DP table values, store direction for traceback
            if max_val == match_mismatch: # Set pointer stored in traceback depending on which condition gave the max
                traceback_dp[i][j] = "match_mismatch"
            elif max_val == seq1_gap:
                traceback_dp[i][j] = "seq1_gap"
            elif max_val == seq2_gap:
                traceback_dp[i][j] = "seq2_gap"
    
    # Call traceback function with traceback_dp to construct alignment strings
    # Set starting index for traceback to be the bottom right cell in the traceback dp table
    seq1_align, seq2_align = traceback(traceback_dp, len(seq1), len(seq2))
    
    return seq1_align, seq2_align, max_val_dp[len(seq1)][len(seq2)] # Return optimal alignments, and optimal alignment score