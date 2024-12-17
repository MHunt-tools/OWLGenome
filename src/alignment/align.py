## Import packages
import os, sys

def CalcHammingDist(seq1, seq2):
    '''
    Calculates the Hamming Dist. btwn two seqs
    '''
    #mismatches = [i for i, bp in enumerate(seq1) if bp != seq2[i]]     # (Time) complexity O(n) but (Space) complexity is O(k); k stored mismatch list
    if len (seq1) != len(seq2):
        raise ValueError('Sequence lengths must be equal to compute hamming Dist.')
    mismatches = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))    # (Time) complexity O(n) but (Space) is O(1) from O(k); generation only, no mismatch list stored
    return mismatches

def is_match(seq, genome_seq, start, end, threshold=18):
    '''
    Chechs if the query sequence matches the genome sequence w/in allowed mismatch threshold.

    Inputs:
        seq or the read
        genome_seq or full genome sequence
        start position of alignment in genome 
        threshold allowed mismatches

    Returns:
        True if HD b/w seq & genome_seq[start:end] <= threshold
        Otherwise False
    '''
    genome_subseq = genome_seq[start:end]

    # Check for equal lengths
    if len(seq) != len(genome_subseq):
        return False
    # Calc. HD
    hamming_dist = CalcHammingDist(seq, genome_subseq)

    return hamming_dist <= threshold   # True if mismatches w/in threshold
