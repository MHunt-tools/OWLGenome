import numpy as np
from numba import njit, prange

@njit(cache=True)
def banded_sw_align_numba(seq1_array, seq2_array, n, m, half_bandwidth, gap_penalty, match_score, mismatch_score):
    """
    Performs a banded smith-waterman (sw) local alignment using numpy arrays and numba enabled just in time machine code compilation
    banded_sw_align_numba constructs the dynamic programming (dp) table as a numpy array for use in a seperate traceback function
    Function is only compiled on first call, and is thereafter cached

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
        max_score - integer value representing the maximum alignment score found for the sw global alignment between seq1_array, seq2_array
        max_pos - 2 len tuple (i, j) containing the index within traceback_dp that max_score was found in the sw dp table
        traceback_dp - numpy 2D array that contains traceback pointers for every cell in max_val_dp for use in a seperate traceback function
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

    # Fill in the base case values for the first column in the dp table based off sw rules
    for i in range(1, n + 1):
        j_band = half_bandwidth - i + 0  # Convert j to j_band for the column indexed by j = 0
        if 0 <= j_band <= 2 * half_bandwidth: # Restrict the bound of j_band to not go out of bounds of the band
            max_val_dp[i, j_band] = 0 # All j = 0 values are 0 per sw rules
            traceback_dp[i, j_band] = END # Fill traceback array with the pointer for ending the alignment 
        else:
            break # Break loop when bandwidth bound broken

    # Fill in the base case values for the first row in the dp table based off sw rules
    for j in range(1, m + 1):
        j_band = half_bandwidth - 0 + j # Convert j to j_band for the column indexed by i = 0
        if 0 <= j_band <= 2 * half_bandwidth: # Restrict the bound of j_band to not go out of bounds of the band
            max_val_dp[0, j_band] = 0 # Fill in the first row according to the gap_penalty * i base case rule of nw
            traceback_dp[0, j_band] = END # Fill traceback array with the pointer for ending the alignment
        else:
            break # Break loop when bandwidth bound broken

    # Sets i, j = 0, 0 value and traceback pointer
    max_val_dp[0, half_bandwidth] = 0 
    traceback_dp[0, half_bandwidth] = END

    # Initialize values to hold max_score and max_pos, which will update as the table is filled
    max_score = 0
    max_pos = (0,0)

    # Fill the DP table
    for i in range(1, n + 1):
        start_j = max(1, i - half_bandwidth) # Get j index for start of band
        end_j = min(m, i + half_bandwidth) # Get j index for end of band

        for j in range(start_j, end_j + 1):
            j_band = half_bandwidth - i + j # Convert j indices into j_band indices

            # Calculate SW score for match or mismtach
            if 0 <= j_band <= 2 * half_bandwidth:
                if seq1_array[i - 1] == seq2_array[j - 1]:
                    score = match_score
                else:
                    score = mismatch_score
                match_mismatch = max_val_dp[i - 1, j_band] + score
            else:
                match_mismatch = float('-inf') # Set scores to -inf if an out of bounds j_band is accessed

            # Calculate SW score for a gap in seq1
            if 0 <= j_band - 1 <= 2 * half_bandwidth:
                seq1_gap = max_val_dp[i, j_band - 1] + gap_penalty
            else:
                seq1_gap = float('-inf')

            # Calculate SW score for a gap in seq2
            if 0 <= j_band + 1 <= 2 * half_bandwidth:
                seq2_gap = max_val_dp[i - 1, j_band + 1] + gap_penalty
            else:
                seq2_gap = float('-inf')

            # Find which movement gives the maximum SW score with a floor of 0 and set that as the value for the current cell of the SW dp table
            max_val = max(match_mismatch, seq1_gap, seq2_gap, 0)
            max_val_dp[i, j_band] = max_val

            # Set traceback pointers for the traceback table based off which move resulted in the maximum NW score
            if max_val == 0:
                traceback_dp[i, j_band] = END
            elif max_val == match_mismatch:
                traceback_dp[i, j_band] = MATCH_MISMATCH
            elif max_val == seq1_gap:
                traceback_dp[i, j_band] = GAP_IN_SEQ1
            elif max_val == seq2_gap:
                traceback_dp[i, j_band] = GAP_IN_SEQ2

            # Update the max_val if a greater max_val is found
            if max_val > max_score:
                max_score = max_val
                max_pos = (i ,j) # Update the position as well

    return max_score, max_pos, traceback_dp # Return max_score, max_pos, and the traceback table

@njit(cache=True)
def banded_traceback_numba(seq1_array, seq2_array, traceback_dp, half_bandwidth, max_pos):
    """
    Using a pre-filled banded smith-waterman (sw) dynamic programming (dp) traceback table, construct the local alignment between seq1_array and seq2_array as alignment arrays containing DNA chracters as ascii values
    traceback_numba utilizes numpy arrays and numba enabled just in time machine code compilation for efficient computation
    Function is only compiled on first call, and is thereafter cached

    INPUTS:
        seq1_array - ascii DNA character numpy array of DNA sequence to be aligned
        seq2_array - ascii DNA character numpy array of DNA sequence to be aligned
        traceback_dp - numpy 2D array that contains traceback pointers for reconstruction of the local alignment between seq1_array and seq2_array
        half_bandwidth - integer distance from central band of dynamic programming table to edge of band, whole bandwidth is 2 * half_bandwidth + 1
        max_pos - 2 len tuple (i, j) containing location of maximum value found in the sw local alignment, which corresponds to the same cell in traceback_dp
        
    OUTPUTS:
        aligned_seq1 - ascii DNA character array for the seq1 global alignment
        alinged_seq2 - ascii DNA character array for the seq2 global alignment
    """

    # Extract max position indices
    i, j = max_pos

    # Set traceback pointer values
    END = 0
    MATCH_MISMATCH = 1
    GAP_IN_SEQ1 = 2
    GAP_IN_SEQ2 = 3 

    # Initialize numpy arrays to hold the alignment ascii character values
    max_pos_len = i + j # Maximum possible alignment character array length
    aligned_seq1 = np.empty(max_pos_len, dtype=np.uint8)
    aligned_seq2 = np.empty(max_pos_len, dtype=np.uint8)
    index = max_pos_len - 1 # Initialize an index for the last value in the alignment arrays

    # Iterate through the traceback table and construct the alignment string
    while i > 0 and j > 0:
        j_band = j - i + half_bandwidth # Convert j to j_band

        if j_band < 0 or j_band >= 2 * half_bandwidth + 1: # End alignment if traceback takes us outside of bandwidth
            break

        traceback_code = traceback_dp[i, j_band] # Extract traceback code 

        # Based off of the alignment code, append the DNA string characters, or a gap, to the alignment string character arrays
        if traceback_code == END: # End alignment if end pointer is found
            break
        elif traceback_code == MATCH_MISMATCH:
            aligned_seq1[index] = seq1_array[i - 1]
            aligned_seq2[index] = seq2_array[j - 1]
            i -= 1
            j -= 1
        elif traceback_code == GAP_IN_SEQ1:
            aligned_seq1[index] = ord('_')
            aligned_seq2[index] = seq2_array[j - 1]
            j -= 1
        elif traceback_code == GAP_IN_SEQ2:
            aligned_seq1[index] = seq1_array[i -1]
            aligned_seq2[index] = ord('_')
            i -= 1
        else:
            break
        index -= 1 # Decrement alignment character arry index for every mvoe through traceback table

    rev_ind = index + 1 # Slice based off last aligned position
    # Reverse and slice arrays
    aligned_seq1 = aligned_seq1[rev_ind::]
    aligned_seq2 = aligned_seq2[rev_ind::]

    return aligned_seq1, aligned_seq2 # Return alignments as ascii character value arrays

@njit(cache = True)
def sw_align_numba(seq1_array, seq2_array, gap_penalty, match_score, mismatch_score):
    '''
    Function sw_align_numba uses the numba package to implement just-in-time compilation for creating the smith-waterman dynamic programming table
    Compilation occurs at runtime in the first function call and is cached for all subsequent calls (cache = True)

        INPUTS:
            seq1_array: numpy array of DNA characters converted into ascii values
            seq2_array: numpy array of DNA characters converted into ascii values
            gap_penalty - integer gap penalty value for nw dp table filling
            match_score - integer match score value for nw dp table filling
            mismatch_score - integer mismatch_score value for nw dp table filling


        OUTPUTS:
            max_score - integer value representing the maximum alignment score found for the sw global alignment between seq1_array, seq2_array
            max_pos - 2 len tuple (i, j) containing the index within traceback_dp that max_score was found in the sw dp table
            traceback_dp - numpy 2D array that contains traceback pointers for every cell in max_val_dp for use in a seperate traceback function

    '''

    # Get lengths of sequence arrays
    n = len(seq1_array)
    m = len(seq2_array)

    # Initialize numpy arrays for use as dp and traceback tables
    max_val_dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    traceback_dp = np.zeros((n + 1, m + 1), dtype=np.int8)

    # Set pointer values for traceback
    MATCH_MISTMATCH = 1
    GAP_IN_SEQ1 = 2
    GAP_IN_SEQ2 = 3
    END = 0

    # Initialize values to hold max_score and max_pos, which will update as the table is filled
    max_score = 0
    max_pos = (0,0)

    # Fill the DP table
    for i in range(1, n+1):
        seq1_i = seq1_array[i-1] # get current char for seq1
        for j in range(1, m + 1):
            seq2_j = seq2_array[j-1] # get current char for seq2

            # Calculate match/mismatch score
            if seq1_i == seq2_j:
                score = match_score
            else:
                score = mismatch_score

            # Calculate SW scores
            match_mismatch = max_val_dp[i - 1, j - 1] + score
            seq1_gap = max_val_dp[i, j - 1] + gap_penalty
            seq2_gap = max_val_dp[i - 1, j] + gap_penalty

            # Find which movement gives the max score with a floor of zero
            max_val = max(match_mismatch, seq1_gap, seq2_gap, 0)

            # Set dp table val to be the maximum score
            max_val_dp[i, j] = max_val

            # Set traceback pointers
            if max_val == 0:
                traceback_dp[i, j] = END
            elif max_val == match_mismatch:
                traceback_dp[i, j] = MATCH_MISTMATCH
            elif max_val == seq1_gap:
                traceback_dp[i, j] = GAP_IN_SEQ1
            elif max_val == seq2_gap:
                traceback_dp[i, j] = GAP_IN_SEQ2

            # Update max score and position
            if max_val > max_score:
                max_score = max_val
                max_pos = (i ,j)

    
    return max_score, max_pos, traceback_dp # Return max score, position, and traceback_dp
        

@njit(cache=True)
def traceback_numba(seq1_array, seq2_array, traceback_dp, max_pos):
    """
    Perform traceback for the full (non-banded) Smith-Waterman alignment using numpy arrays and numba just in time compilation to machine code
    
    INPUTS:
        seq1_array: numpy array of DNA characters converted into ascii values
        seq2_array: numpy array of DNA characters converted into ascii values
        traceback_dp: numpy 2D array that contains traceback pointers for the optimal local alignment between seq1 and seq2
        max_pos (tuple): (i,j) position of maximum score in DP/traceback table.

    OUTPUTS:
        aligned_seq1 - ascii DNA character array for the seq1 global alignment
        alinged_seq2 - ascii DNA character array for the seq2 global alignment
    """
    # Extract max pos indices
    i, j = max_pos

    # Set traceback codes
    END = 0
    MATCH_MISMATCH = 1
    GAP_IN_SEQ1 = 2
    GAP_IN_SEQ2 = 3

    # Get maximum possible alignment length
    max_pos_len = i + j

    # Initialize numpy arrays for aligned sequences
    aligned_seq1 = np.empty(max_pos_len, dtype=np.uint8)
    aligned_seq2 = np.empty(max_pos_len, dtype=np.uint8)
    index = max_pos_len - 1 # Set index to last entry in array

    # Iterate through traceback table
    while i > 0 and j > 0:
        traceback_code = traceback_dp[i, j] # Extract traceback code

        if traceback_code == END: # End alignment if end pointer found
            break
        elif traceback_code == MATCH_MISMATCH: # match/mismatch
            aligned_seq1[index] = seq1_array[i - 1]
            aligned_seq2[index] = seq2_array[j - 1]
            i -= 1
            j -= 1
        elif traceback_code == GAP_IN_SEQ1:
            aligned_seq1[index] = ord('_')  # Gap in seq1
            aligned_seq2[index] = seq2_array[j - 1]
            j -= 1
        elif traceback_code == GAP_IN_SEQ2:
            aligned_seq1[index] = seq1_array[i - 1]
            aligned_seq2[index] = ord('_')  # Gap in seq2
            i -= 1
        else:
            break

        index -= 1 # Decrement the alignment array index every iteration

    # Slice and reverse the alignment arrays to include only the aligned portion in the correct orientation
    rev_ind = index + 1
    aligned_seq1 = aligned_seq1[rev_ind:]
    aligned_seq2 = aligned_seq2[rev_ind:]

    return aligned_seq1, aligned_seq2 # Return alignment arrays 

def traceback(seq1_array, seq2_array, traceback_dp, max_pos):
    """
    Perform traceback for the full (non-banded) Smith-Waterman alignment run with numba njit
    
    INPUTS:
        seq1_array: numpy array of DNA characters converted into ascii values
        seq2_array: numpy array of DNA characters converted into ascii values
        traceback_dp: numpy 2D array that contains traceback pointers for the optimal local alignment between seq1 and seq2
        max_pos (tuple): (i,j) position of maximum score in DP/traceback table.

    OUTPUTS:
        aligned_seq1 - DNA string for the seq1 global alignment
        alinged_seq2 - DNA string for the seq2 global alignment
    """
    # Initialize arrays to hold alignment characters
    aligned_seq1 = []
    aligned_seq2 = []

    # Extract max pos indices
    i, j = max_pos

    # Set traceback pointers
    MATCH_MISTMATCH = 1
    GAP_IN_SEQ1 = 2
    GAP_IN_SEQ2 = 3
    END = 0

    # Iterate through the traceback table and construct the alignment string
    while i > 0 and j > 0:
        traceback_code = traceback_dp[i, j]  # Extract traceback code 
        if traceback_code == END:  # End alignment if end pointer is found
            break
        # Append gap or character for other pointers, non gap characters have to be converted back into characters from ascii values
        elif traceback_code == MATCH_MISTMATCH:
            aligned_seq1.append(chr(seq1_array[i - 1]))
            aligned_seq2.append(chr(seq2_array[j - 1]))
            i -= 1
            j -= 1
        elif traceback_code == GAP_IN_SEQ1:
            aligned_seq1.append('_')
            aligned_seq2.append(chr(seq2_array[j - 1]))
            j -= 1
        elif traceback_code == GAP_IN_SEQ2:
            aligned_seq1.append(chr(seq1_array[i - 1]))
            aligned_seq2.append('_')
            i -= 1

    # Join and reverse character arrays to form aligned sequence strings
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return aligned_seq1, aligned_seq2 # Return aligned sequence strings

def banded_sw_align(seq1, seq2, half_bandwidth, gap_penalty=-1, match_score=1, mismatch_score=-1):
    '''
    Function banded_sw_align performs a banded smith waterman local alignment of two input DNA strings

        INPUT
            seq1: DNA string to be aligned
            seq2: DNA string to be aligned
            gap_penalty: integer value for the gap penalty
            match_score: integer value for a match score
            mismatch_score: integer value for a mismatch score
            half_bandwidth: integer value specifying how far out from the main diagonal to fill in the DP table. Total bandwidth is 2 * half_bandwidth + 1
        OUTPUT
            seqs[0]: The aligned DNA string for seq1
            seqs[1]: The aligned DNA string for seq2
            max_pointer
    '''
    
    n, m = len(seq1), len(seq2) # Extract the lengths of seq1 and seq2

    # Create traceback_dp and match_mistmatch_dp tables containing n+1 rows as in standard SW
    # However, banded tables will only contain 2 * half_bandwidth + 1 rows
    # This will only include the diagonal region of width half_bandwidth around the central diagonal
    # 2 * half_bandwidth to include both side of the diagonal
    # + 1 to include the main diagonal cell itself

    traceback_dp = [["end"] * (2 * half_bandwidth + 1) for _ in range(n+1)] # Create a DP table with n+1 rows and 2 * half_bandwidth * 1 columns
    max_val_dp = [[0] * (2 * half_bandwidth + 1) for _ in range(n+1)] # Create a DP table with n+1 rows and 2 * half_bandwidth * 1 columns


    def score(a, b): # Simple scoring function
        '''
        Simple scoring function for character match or mismatches

        INPUTS:
            a - string character to be compared
            b - string character to be comapred

        OUTPUTS:
            match_score - integer score for if the character are the same, fed from banded_sw_align
            mismatch_score - integer score for if the characters are different, fed from banded_sw_align
        '''
    
        return match_score if a == b else mismatch_score
    

    def traceback(traceback_dp, i, j_band):
        """
        Using a banded sw dp traceback table, construct the local alignment between seq1 and seq2 as alignment strings

        INPUTS:
            traceback_dp - 2D array that contains traceback pointers for reconstruction of the local alignment between seq1 and seq2
            i - seq1 max pos index
            j_band - seq2 max pos index within the band

        OUTPUTS:
            seq1_align - string of DNA for the seq1 local alignment
            seq2_align - string of DNA for the seq2 local alignment
        """
        seq_1_align = '' # Initialize strings
        seq_2_align = ''
        while i > 0 and 0 <= j_band <= 2 * half_bandwidth: # Set up while loop break condition to end when i hits zero or the j_band goes out of the banded region
            
            if traceback_dp[i][j_band] == "end": # End local alignment when an end pointer is hit
                break

            j = j_band - half_bandwidth + i # Convert j_band to j to allow for indexing into seq2
            
            if traceback_dp[i][j_band] == "seq1_gap": # if a gap in seq1 is found, then append a gap to aligned seq1 and append the curret seq2 char to aligned seq2
                seq_1_align = "_" + seq_1_align
                seq_2_align = seq2[j - 1] + seq_2_align
                j_band -= 1 # Move to j_band - 1
            elif traceback_dp[i][j_band] == "match_mismatch": # if a match/mismatch is found, append the current char for seq1 and seq2 to their aligned sequences
                seq_1_align = seq1[i - 1] + seq_1_align
                seq_2_align = seq2[j - 1] + seq_2_align
                i -= 1 # Move back diagonally
                # Diagonal move cancels out the shift to j_band from moving up a row, so j_band is not updated
            elif traceback_dp[i][j_band] == "seq2_gap":  # if a gap in seq2 is found, then append a gap to aligned seq 2 and append the curret seq1 char to aligned seq1
                seq_1_align = seq1[i - 1] + seq_1_align
                seq_2_align = "_" + seq_2_align
                i -= 1 # Move to i - 1
                j_band += 1 # Moving up a row results in a right shift of the j_band index

        return seq_1_align, seq_2_align # Returned the constructed alignments
    
    
    max_pointer = float('-inf') # Initiate value to hold largest value found in max_val_dp
    max_index = (0,0) # Initiate tuple for holding where the maximum index is found in max_val_dp

    for i in range(1, n +1): # Iterates over rows, or positions in seq1
        start_j = max(1, i - half_bandwidth) # start of j range set to leftmost position within band at row i, min is one for first row
        end_j = min(m, i + half_bandwidth) # end of j range set to rightmost position within band at row i, max at m so index does not exceed seq2 length
        for j in range(start_j, end_j + 1): # For the current row, iterate across possible j positions
            j_band = j-i+half_bandwidth # Convert the j value to j_band, which will actually allow access within the row as each row is in range 0, 2 * half_bandwidth

            if 0 <= j_band <= 2 * half_bandwidth: # Check to see if previous j position is within the bandwidth
                # Previous j_band value will just be the same j_band for a diagonal move
                match_mismatch = max_val_dp[i - 1][j_band] + score(seq1[i - 1], seq2[j - 1]) # If so, calculate match_mismatch value
            else: 
                match_mismatch = float('-inf')  # If the previous j position is out of the bandwidth, set to -inf

            if 0 <= j_band - 1 <= 2 * half_bandwidth: # Check to see if previous j position is within the andwidth
                # A move to the left column is just j_band - 1
                seq1_gap = max_val_dp[i][j_band - 1] + gap_penalty # If so, calculate the seq1_gap value
            else:
                seq1_gap = float('-inf') # If the previous j position is out of the bandwidth, set to -inf

            if 0 <= j_band + 1 <= 2 * half_bandwidth: # Check to see if previous j position is within the bandwidth, no movement in j so tested value is just j_band
                # A move up a row results in j_band + 1
                seq2_gap = max_val_dp[i - 1][j_band+1] + gap_penalty # If so, calculate the seq2_gap value
            else:
                seq2_gap = float('-inf') # If the previous j position is out of the bandwidth, set to -inf

            max_val = max(match_mismatch, seq1_gap, seq2_gap, 0) # Find the move that gives the maximum score with a floor of zero
            max_val_dp[i][j_band] = max_val # Set the current position in max_val_dp to be that maximum score

            if max_val > max_pointer: # Update the maximum value holder across the loop to keep track of what and where it is
                max_pointer = max_val
                max_index = (i, j_band)
            
            if max_val == match_mismatch: # Set pointer stored in traceback depending on which condition gave the max
                traceback_dp[i][j_band] = "match_mismatch"
            elif max_val == seq1_gap:
                traceback_dp[i][j_band] = "seq1_gap"
            elif max_val == seq2_gap:
                traceback_dp[i][j_band] = "seq2_gap"
            else:
                traceback_dp[i][j_band] = "end"

    seqs = traceback(traceback_dp, max_index[0], max_index[1]) # Call traceback to construct alignment sequences
            
    return seqs[0], seqs[1], max_pointer # Returned seq1 alignment, seq2 alignment, and the optimal alignment score

def sw_align(seq1, seq2, gap_penalty, match_score, mismatch_score):
    '''
    Function sw_align performs a smith waterman local alignment of two input DNA strings

        INPUT:
            seq1: DNA string to be aligned
            seq2: DNA string to be aligned
            gap_penalty: integer value for the gap penalty
            match_score: integer value for a match score
            mismatch_score: integer value for a mismatch score

        OUTPUT:
            seqs[0]: The aligned DNA string for seq1
            seqs[1]: The aligned DNA string for seq2
            max_pointer
    '''

    # DP tables are initialized to hold score and traceback information
    traceback_dp = [["end"] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] # This table will keep backtracking pointers for reconstructing the sequence
    max_val_dp = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)] # This table will be for the best alignment ending with a match/mismatch
    
    
    def score(a,b, match_score=match_score, mismatch_score=mismatch_score):
        '''
        Simple scoring function for character match or mismatches

        INPUTS:
            a - string character to be compared
            b - string character to be comapred

        OUTPUTS:
            match_score - integer score for if the character are the same, fed from banded_sw_align
            mismatch_score - integer score for if the characters are different, fed from banded_sw_align
        '''
        
        if a == b:
            return match_score

        else:
            return mismatch_score


    def traceback(traceback_dp, i, j):
        """
        Using a sw dp traceback table, construct the local alignment between seq1 and seq2 as alignment strings

        INPUTS:
            traceback_dp - 2D array that contains traceback pointers for reconstruction of the local alignment between seq1 and seq2
            i - seq1 max pos index
            j - seq2 max pos index

        OUTPUTS:
            seq1_align - string of DNA for the seq1 local alignment
            seq2_align - string of DNA for the seq2 local alignment
        """
        seq_1_align = '' # Initialize strings
        seq_2_align = ''
        while i > 0 or j > 0: # Set up while loop break condition to end when eithe index is zero
            if traceback_dp[i][j] == "end": # Base case, end local alignment when an end pointer is hit
                i = 0 # Update indices to break loop
                j = 0
                continue
            else:
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
    
    
    max_pointer = float('-inf') # Initiate value to hold largest value found in table
    max_index = (0,0) # Initiate tuple for holding where the maximum index is found

    # Set dimension of tables to be 1 larger in length than the string for both strings
    # Iterate through every character in both input strings
    for i in range(1,len(seq1)+1): 
        for j in range(1,len(seq2)+1): # Range starting at 1 as array construction pre initializes all base cases, i=0 j=0, to 0
        
            # Calculate score for a match/mismatch
            match_mismatch = max_val_dp[i-1][j-1] + score(seq1[i-1],seq2[j-1])

            # Calculate score for gap in seq1
            seq1_gap = max_val_dp[i][j-1] + gap_penalty
            
            # Calculate score for gap in seq2
            seq2_gap = max_val_dp[i-1][j] + gap_penalty
          
            # Fill the match/mismatch table with the max value of the current index across all three possibilies with a floor of zero
            max_val = max(match_mismatch,seq1_gap,seq2_gap,0)
            max_val_dp[i][j] = max_val

            if max_val > max_pointer: # Update the maximum value holder across the loop to keep track of what and where it is
                max_pointer = max_val
                max_index = (i,j)

            if max_val == match_mismatch: # Set pointer stored in traceback depending on which condition gave the max
                traceback_dp[i][j] = "match_mismatch"
            elif max_val == seq1_gap:
                traceback_dp[i][j] = "seq1_gap"
            elif max_val == seq2_gap:
                traceback_dp[i][j] = "seq2_gap"
            elif max_val == 0:
                traceback_dp[i][j] = "end"

 
    # Generate optimal SW local alignment from DP table and format function outputs
    seqs = traceback(traceback_dp, max_index[0], max_index[1])

    return seqs[0],seqs[1],max_pointer # Return seq1 alignment, seq2 alignment, maximum value