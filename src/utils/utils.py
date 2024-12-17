# src/utils/utils.py

import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import random
import multiprocessing
import logging 
import concurrent.futures
import time
from tqdm import tqdm
import pandas as pd
from collections import deque
import numba
from numba import njit
import pysam

# Adding root dir. of the project (OWLGENOME) to sys.path so Python can find the src module
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Import the required modules
from src.alignment.nw_align import banded_nw_align_numba, banded_nw_traceback  # Import the actual function from nw_align.py
from src.alignment.sw_align import sw_align  # Import the actual function from sw_align.py
from src.alignment.sw_align import banded_sw_align # Import the banded_sw_align function from sw_align.py
from src.alignment.sw_align import sw_align_numba, traceback, traceback_numba, banded_traceback_numba, banded_sw_align_numba


def RevCom(seq):
    '''
    Returns the reverse complement of a dna sequence
    '''
    com_dict = {'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G',
                'N': 'N'}
    revcom = ''.join([com_dict.get(base, 'N') for base in reversed(seq)])   # (memory) seq returns an iterator to return sequence for long reads; creating a reverse string uses additional mem
    return revcom

def corrupt_seq(seq, error_rate=0.15, indel_rate=0.05):
    '''
    Corrupt a DNA sequence based on given error and indel rates
    '''
    corrupted_seq = seq
    for pos in range(len(seq)):
        if random.random() < error_rate:
            corrupted_seq = corrupted_seq[:pos] + np.random.choice(['A', 'G', 'T', 'C']) + corrupted_seq[pos+1:]

    indices = range(int(len(seq) * indel_rate))
    indel_idx = [int(random.random() * len(seq)) for i in indices]
    for pos in indel_idx:
        insert = random.random() < 0.5
        if insert:
            corrupted_seq = corrupted_seq[:pos+1] + np.random.choice(['A', 'G', 'T', 'C']) + corrupted_seq[pos+1:]
        else:
            corrupted_seq = corrupted_seq[:pos] + corrupted_seq[pos+1:]
    return corrupted_seq


def seq_to_kmers(seq, k_size, sort=False, include_revcom=True, k_stride=5):    # (Time) complexity O(n) from O(n*k), (Space) O(k); k = # k-mers gen.
    '''
    Splits the sequence into kmers of length k_size along a sliding window

    If include_revcom: compute reverse complement for each seq and add to list
    '''
    indices = range(0, len(seq) - k_size + 1, k_stride)
    if not include_revcom:
        kmers = [seq[i:i + k_size] for i in indices]
    else: 
        kmers = []
        for i in indices:
            kmers += [seq[i:i + k_size]]
            kmers += [RevCom(seq[i:i + k_size])]
    if sort: 
        kmers.sort()
    return kmers


def seq_to_idx_kmers(seq, k_size, k_stride=1, include_revcom=False):
    '''
    Splits the sequence into an list of indexed kmers [(pos, seq)] of length k_size along a sliding window

    If include_revcom: compute reverse complement for each seq and add to list with same indexing as forward seq
    '''
    indices = range(0, len(seq) - k_size + 1, k_stride) 
    if not include_revcom:
        idx_kmers = [(i, seq[i:i + k_size], 1) for i in indices]
    else: 
        idx_kmers = []
        for i in indices:
            kmer = seq[i:i + k_size]
            idx_kmers += [(i, kmer, 1)]
            idx_kmers += [(i, RevCom(kmer), -1)]
    return idx_kmers


def file_to_idx_kmers(file_path, file_fmt, k_size, k_stride=1, include_revcom=False):
    '''
    Reads a fasta file and splits sequences into kmers of k_size

    Returns: 
        A list of tuples of (record id, indexed k-mers) for each record in the fasta file
        Indexed kmers are tuples of (pos, seq) for each position and kmer sequence in the record
    '''
    record_kmers = []
    for ri, record in enumerate(SeqIO.parse(file_path, format=file_fmt)):
        idx_kmers = seq_to_idx_kmers(str(record.seq.upper()), k_size=k_size, k_stride=k_stride, include_revcom=include_revcom)
        record_kmers.append((record.id, idx_kmers))
    return record_kmers

def coarse_read_map(seq, k_size, hash_table, k_stride=5, include_revcom=False, dropout_rate=0, chain_method='lis'):
    '''
    Coarse mapping of the read sequence against the genome hash table.
    
    Inputs:
        seq: read to map
        k_size: size of the kmers
        hash_table: genome hash table of kmers
        k_stride: stride for sliding window over the read
        include_revcom: if True, includes kmers from the reverse complement strand
        dropout_rate: proportion of kmers to randomly drop (for sampling)
        
    Returns:
        List of possible genome location ranges where the read might align.
    '''
    # Generate indexed kmers from the read sequence
    idx_kmers = seq_to_idx_kmers(seq, k_size, k_stride=k_stride, include_revcom=include_revcom)

    ## Retreive, filter, and flatten kmers
    kmer_locs = [(idx, genome_pos, genome_strand) for read_pos, (idx, kmer, strand) in enumerate(idx_kmers) if kmer in hash_table for genome_pos, genome_strand in hash_table[kmer]]

    if not kmer_locs:   # Return empty list if no kmers are found in the hash table
        return [], [], []

    # Get possible chains of genome positions where the kmers might align
    chains = get_possible_chains(kmer_locs, k_size, max_gap=len(seq) * 1.1, min_length=len(seq) * 0.5, method=chain_method)

    # Exit if no possible chains are found
    if not chains:
        return [], [], []

    # Calculate possible genome location ranges from the chains
    genome_locations = []
    for chain in chains:    # (Time) complexity O(k); k-mers in chain
        # Find min & max genome positions in chain
        min_genome_pos = min(chain[0][1], chain[-1][1])
        max_genome_pos = max(chain[0][1], chain[-1][1])
        min_position = min(chain[0][0], chain[-1][0])
        max_position = max(chain[0][0], chain[-1][0])
        seq_revcom = True if min_position != chain[0][0] else False

        # Calculate the offset between read and genome positions
        if not seq_revcom:
            start = min_genome_pos - min_position
            end = max_genome_pos + (len(seq) - max_position)
            genome_locations.append((start, end))
        else:
            offset = len(seq) - max_position - k_size
            start = min_genome_pos - offset
            end = max_genome_pos + (len(seq) - min_position) + k_size
            end = start + len(seq)
            genome_locations.append((end, start))

    return genome_locations, chains, idx_kmers

@njit
def get_lis_numba(locs):
    '''
    Numba implementation of the get_list function
    Get longest increaseing subsequence from mapped locations
    '''
    ## Reshape array
    n = len(locs)  
    arr = np.empty(n, dtype=np.int64)
    for i in range(n):
        arr[i] = locs[i][0]

    ## Initialize dp
    dp = np.ones(n, dtype=np.int64)
    seq = np.arange(n, dtype=np.int64)

    ## Fill dp table
    for i in range(n):
        for prev in range(i):
            if arr[prev] < arr[i] and dp[prev] + 1 > dp[i]:
                dp[i] = dp[prev] + 1
                seq[i] = prev

    ## Find the best result
    ans = dp[0]
    ans_ind = 0
    for i in range(1, n):
        if dp[i] > ans:
            ans = dp[i]
            ans_ind = i

    ## Construct the result sequence
    length = dp[ans_ind]
    res = np.empty(length, dtype=np.int64)
    idx = ans_ind
    for i in range(length - 1, -1, -1):
        res[i] = arr[idx]
        if seq[idx] == idx:
            break
        idx = seq[idx]

    return res

def get_lis(locs):
    '''
    Get longest increaseing subsequence from mapped locations
    Adapted from https://www.geeksforgeeks.org/construction-of-longest-increasing-subsequence-using-dynamic-programming/
    '''
    ## Reshape array
    arr = [loc[0] for loc in locs]
    n = len(arr)

    # Initialize dp array with 1
    dp = [1] * n

    # Initialize seq array with index values (to store previous indices in LIS)
    seq = list(range(n))

    # Compute dp and seq arrays
    for i in range(n):
        for prev in range(i):
            if arr[prev] < arr[i] and 1 + dp[prev] > dp[i]:
                dp[i] = 1 + dp[prev]
                seq[i] = prev

    # Find the index of the last element in the LIS
    ans = -1
    ans_ind = -1
    for i in range(n):
        if dp[i] > ans:
            ans = dp[i]
            ans_ind = i

    # Construct the result sequence using seq array
    res = []
    res.append(arr[ans_ind])
    while seq[ans_ind] != ans_ind:
        ans_ind = seq[ans_ind]
        res.append(arr[ans_ind])

    # Reverse the result to get the correct order
    res.reverse()
    return res


def get_possible_chains(kmer_locs, k_size, max_gap=10000, min_length=100, method='hough'):   
    '''
    Builds possible chains of kmer positions based on proximity in the genome.
    
    Inputs:
        kmer_locs: List of lists containing genome positions for each kmer.
                   Example: [[101, 543, 956], [156, 557, 958], ...]
        max_gap: Maximum allowed gap (in genome positions) between consecutive kmers in a chain.
                 Determines how close k-mers must be to form a valid chain.
        
    Returns:
        A list of chains, each chain is a list of tuples (read_pos, genome_pos, strand).
        Example: [
            [(0, 101, 1), (1, 154, -1), (2, 203, 1)],
            [(0, 505, -1), (1, 555, 1), (2, 656, 1)],
            ...
        ]
    '''
    # Sort locs by genome position
    f_kmer_locs = []
    r_kmer_locs = []
    for loc in kmer_locs:
        if loc[2] == 1:
            f_kmer_locs.append(loc)
        else:
            r_kmer_locs.append(loc)

    ## Set longest number of hits as true
    kmer_locs = f_kmer_locs if len(f_kmer_locs) > len(r_kmer_locs) else r_kmer_locs
    kmer_locs.sort(key=lambda x: x[1])  # Sort by the second element (genome position)

    if abs(kmer_locs[-1][1]-kmer_locs[0][1]) < max_gap: ## return whole thing if within max gap
        return [kmer_locs]

    if method == 'lis':
        ## Get longest common subsequence
        if len(f_kmer_locs) > len(r_kmer_locs):
            kmer_lcs = get_lis_numba(kmer_locs)
        else:
            kmer_lcs = get_lis_numba(kmer_locs[::-1])
        et = time.time()
        return [[loc for loc in kmer_locs if loc[0] in kmer_lcs]]
    elif method == 'hough':
        ## Get base positions
        if len(f_kmer_locs) > len(r_kmer_locs):
            ## Find "start" position
            pos_arr = np.array([i[1] - i[0] for i in kmer_locs])
        else:
            ## Find "end" position
            pos_arr = np.array([i[1] + i[0] for i in kmer_locs])

        ## Find most common position
        positions, counts = np.unique(pos_arr, return_counts=True)
        best_pos = positions[np.argmax(counts)]
        retreive_kmers = abs(pos_arr - best_pos) < 100
        
        return [[loc for li, loc in enumerate(kmer_locs) if retreive_kmers[li]]]

def generate_dna(l):
    '''
    Generate random DNA string of lenght l
    '''
    return ''.join(random.choice(['A', 'G', 'T', 'C']) for i in range(l))


def adjust_start(input_seq, genome_seq, output_pos):
    '''
    Run small local alignment to adjust the start position of the reads
    Inputs: 
        input_seq: query sequence
        genome_seq: genome sequence
        output_pos: output from coarse grained read mapping
    '''
    ## Genome back-trace
    genome_neg_pos = 30

    ## Get sequences
    start = min(output_pos)
    if output_pos[0] == start: ## Not reverse complement
        seq1 = input_seq[:genome_neg_pos]
        seq2 = genome_seq[output_pos[0]-genome_neg_pos:output_pos[0]+(2*genome_neg_pos)]
    else: ### If reverse complement
        seq1 = RevCom(input_seq)[:genome_neg_pos]
        seq2 = genome_seq[start-genome_neg_pos:start+(2*genome_neg_pos)]

    ## Align sequences
    seq1_array = np.frombuffer(seq1.encode('ascii'), dtype=np.uint8)
    seq2_array = np.frombuffer(seq2.encode('ascii'), dtype=np.uint8)
    score, max_pos, traceback_dp = sw_align_numba(seq1_array, seq2_array, -1, 1, -1)
    aln_1_ascii, aln_2_ascii = traceback_numba(seq1_array, seq2_array, traceback_dp, max_pos)
    aln1 = ''.join(map(chr, aln_1_ascii.tolist()))
    aln2 = ''.join(map(chr, aln_2_ascii.tolist()))
    
    ## Get seq locs
    seq1_loc = seq1.find(''.join([i for i in aln1 if i != '_']))
    seq2_loc = seq2.find(''.join([i for i in aln2 if i != '_']))

    ## Adjust starting position
    adjusted_start = start - genome_neg_pos + (seq2_loc - seq1_loc)
    return adjusted_start, score

def adjust_end(input_seq, genome_seq, output_pos):
    '''
    Run small local alignment to adjust the end position of the reads
    Inputs: 
        input_seq: query sequence
        genome_seq: genome sequence
        output_pos: output from coarse grained read mapping
    '''
    ## Genome back-trace
    genome_neg_pos = 150
    input_neg_pos = 30

    ## Get sequences
    start = min(output_pos)
    end = start + len(input_seq)
    if output_pos[0] == start: ## Not reverse complement
        seq1 = input_seq[-input_neg_pos:]
        seq2 = genome_seq[end-(2*genome_neg_pos):end+genome_neg_pos]
    else: ### If reverse complement
        seq1 = RevCom(input_seq)[-input_neg_pos:]
        seq2 = genome_seq[end-(2*genome_neg_pos):end+genome_neg_pos]

    ## Align sequences
    seq1_array = np.frombuffer(seq1.encode('ascii'), dtype=np.uint8)
    seq2_array = np.frombuffer(seq2.encode('ascii'), dtype=np.uint8)
    score, max_pos, traceback_dp = sw_align_numba(seq1_array, seq2_array, -1, 1, -1)
    aln_1_ascii, aln_2_ascii = traceback_numba(seq1_array, seq2_array, traceback_dp, max_pos)
    aln1 = ''.join(map(chr, aln_1_ascii.tolist()))
    aln2 = ''.join(map(chr, aln_2_ascii.tolist()))

    ## Get seq locs
    seq1_loc = seq1.find(''.join([i for i in aln1 if i != '_']))
    seq2_loc = seq2.find(''.join([i for i in aln2 if i != '_']))

    ## Adjust starting position
    adjusted_end = end + (seq2_loc - seq1_loc) + input_neg_pos - 2*genome_neg_pos
    return adjusted_end, score

def full_alignment(genome_seq, seq, _adjusted_start, _adjusted_end):
    '''
    Function for aligning a full read to a genome slice
    '''
    buffer = 5 # Set a buffer value to ensure genome slice covers read start and end coordinates
    genome_alignment_region = genome_seq[_adjusted_start-buffer:_adjusted_end+buffer] # Extract genome slice
    seq1_array = np.frombuffer(genome_alignment_region.encode('ascii'), dtype=np.uint8) # Convert genome slice to an ascii array
    seq2_array = np.frombuffer(seq.encode('ascii'), dtype=np.uint8) # Convert read to an ascii array

    ## Get lengths of genome slice and read
    n = len(seq1_array)
    m = len(seq2_array)

    half_bandwidth = 4 # Set a half bandwidth value for the banded sw alignment

    # Run banded sw alignment and traceback
    max_score, max_pos, traceback_dp = banded_sw_align_numba(seq1_array, seq2_array,n, m, half_bandwidth, -2, 4, -1)
    aln_1_ascii, aln_2_ascii = banded_traceback_numba(seq1_array, seq2_array, traceback_dp, half_bandwidth, max_pos)

    # Convert outputs back into characters, and join into alignment string
    genome_aln = ''.join(map(chr, aln_1_ascii.tolist()))
    query_aln = ''.join(map(chr, aln_2_ascii.tolist()))

    # Set the alignment score
    score = max_score

    return query_aln, genome_aln, score # Return alignments and score

def process_record(record, genome_seq, k_size, hash_table, k_stride, seq_start=None, dropout_rate=0):
    '''
    Function to process single fasta/fastq record for read mapping
    '''
    ## Process seq
    seq = str(record.seq.upper())
    ## Get coarse location
    output, chains, idx_kmers = coarse_read_map(seq=seq, 
        k_size=k_size, 
        hash_table=hash_table, 
        k_stride=k_stride, 
        dropout_rate=dropout_rate,
    )
    ## Adjust locations based on start seq alignment
    if len(output) == 0: return []
    full_output = []
    
    for output_pos in output: 

        ## Find adjusted start and end positions
        _adjusted_start, score = adjust_start(seq, genome_seq, output_pos)
        _adjusted_end, score = adjust_end(seq, genome_seq, output_pos)

        ## check if read is reversed
        start = min(output_pos)
        if output_pos[0] != start: seq = RevCom(seq)
        
        # Run full alignment of read to genome slice
        query_aln, genome_aln, score = full_alignment(genome_seq, seq, _adjusted_start, _adjusted_end)
        
        ## Append to list
        full_output.append((_adjusted_start, _adjusted_end, query_aln, genome_aln, score))

    return full_output

def process_records_file(file_path, file_fmt, genome_seq, k_size, hash_table, k_stride, dropout_rate=0, use_multiprocessing=True, verbose=True, max_processes=None):
    '''
    Function to process fasta/fastq file for read mapping
    '''
    ## Identify locations of each query sequence
    results = []
    if not use_multiprocessing:
        ## Process each read separately
        for record in tqdm(SeqIO.parse(file_path, file_fmt)): 
            results.append((record.id, process_record(record, genome_seq, k_size, hash_table, k_stride, dropout_rate))) 
    else:
        # Combine the inputs into a list of tuples
        records = [i for i in SeqIO.parse(file_path, file_fmt)] 
        combined_inputs = [(record, genome_seq, k_size, hash_table, k_stride, dropout_rate) for record in records]
        ## Get processes
        process_count = max(min(len(records) // 200, multiprocessing.cpu_count()), 1)
        if max_processes:
            process_count = min(process_count, max_processes)
        ## Process reads in batches
        with multiprocessing.Pool(processes=process_count) as pool:  ## verbose call does not print when in monitor_function
            results = pool.starmap(process_record, tqdm(combined_inputs, total=len(combined_inputs), disable=~verbose), chunksize=len(combined_inputs)//process_count)
        ## Output to results
        results = [(records[ri].id, result) for ri, result in enumerate(results)]
    return results 

def process_ground_truth(filepath):
    '''
    Function to load groundtruth data from filepath
    '''
    df = pd.read_csv(filepath, sep='\t', header=None)
    df.columns = ['id', 'start', 'end']
    df.set_index('id', inplace=True)
    df.start = df.start.astype(float)
    df.end = df.end.astype(float)
    return df

def reshape_process_output_df(output):
    '''
    Function to process positional data from simulated record data
    '''
    record_ids = []
    start_list = []
    end_list = []
    for i in output:
        record_ids.append(i[0])
        if len(i[1]) > 0:
            start_list.append(min(i[1][0][0:2]))
            end_list.append(max(i[1][0][0:2]))
        else:
            start_list.append(None)
            end_list.append(None)

    df = pd.DataFrame([record_ids, start_list, end_list]).transpose()
    df.columns = ['id', 'pred_start', 'pred_end']
    df.set_index('id', inplace=True)
    return df

def true_locs_to_df(locs):
    '''
    Function to process read location data and output df
    '''
    df = pd.DataFrame(locs)
    df.columns = ['id', 'start', 'end']
    df.set_index('id', inplace=True)
    return df

def calc_precision_recall(matched_df, distance_cutoff=5):
    '''
    Function to calculate precision and recall from matched read data
    '''
    ## Check matched_df for different kinds of hits
    matched_df['start_diff'] = matched_df['pred_start'] - matched_df['start']
    matched_df['end_diff'] = matched_df['pred_end'] - matched_df['end']
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0
    for row in matched_df.iterrows():
        calculated_start = row[1]['pred_start']
        true_start =  row[1]['start']
        calculated_end = row[1]['pred_end']
        true_end =  row[1]['end']
        if np.isnan(true_start) and not calculated_start: ## unmappable read is unmapped
            true_negative += 1
        elif np.isnan(true_start) and not not calculated_start: ## unmappable read is mapped
            false_positive += 1
        elif not calculated_start: ## Mappable read is unmapped
            false_negative += 1
        elif (abs(calculated_start-true_start) < distance_cutoff) \
                   and (abs(calculated_end-true_end) < distance_cutoff): ## mappable read is correctly mapped
            true_positive += 1
        else: ## Mappable read is incorrectly mapped
            false_positive += 1

    # Avoid division by zero
    recall = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0
    precision = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0

    ## Write to matched_df.csv
    matched_df.to_csv('matched_df.csv')
    return precision, recall

def filepath_to_filetype(filepath):
    '''
    Function to convert file paths to file types based on file suffixes
    '''
    filetype = filepath.split('.')[-1]
    filefmt = None
    if (filetype == 'fa') or (filetype == 'fasta'):
        filefmt = 'fasta'
    elif (filetype == 'fastq'):
        filefmt = 'fastq'
    return filefmt

def aln_to_cig(query_aln, ref_aln):
    '''
    Function to convert an alignment to a cigar string
    '''
    if not query_aln or not ref_aln:
        return None 
    
    ## Generate expanded cigar string
    cig = []
    for i, nuc in enumerate(query_aln):
        if nuc != '-':
            if ref_aln[i] != '-':
                if ref_aln[i] == nuc:
                    cig.append('M')
                else:
                    cig.append('X')
            else:
                cig.append('I')
        elif ref_aln[i] != '-':
            cig.append('D')

    ## Compress cigar string
    curr_let = cig[0]
    curr_let_count = 0
    compressed_cig = ''
    for i in cig: 
        if i == curr_let:
            curr_let_count += 1
        else:
            compressed_cig += str(curr_let_count) + curr_let
            curr_let = i
            curr_let_count = 1

    return compressed_cig

def convert_to_sam(results, output_path, reference_name="reference_genome", reference_length=1000):
    '''
    Function to convert alignment results to sam file
    Input: 
        results: alignment results
        output_path: path to output sam file
        reference_name: name of reference for SAM file
        reference_length: length of reference for SAM file

    '''
    # Define the SAM file header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': reference_length, 'SN': reference_name}]
    }

    # Open the SAM file for writing
    with pysam.AlignmentFile(output_path, "w", header=header) as outfile:
        for record_id, mappings in results:

            # Skip records with no mappings
            if not mappings:
                # Create a new AlignedSegment object
                read = pysam.AlignedSegment()
                read.query_name = record_id
                read.flag = 4  # Unmapped
                read.reference_id = outfile.get_tid(reference_name)
                read.reference_start = 0  # Convert to 0-based position
                read.mapping_quality = 0  # Default mapping quality
                read.cigarstring = None  # Simplified CIGAR string
                read.query_sequence = None
                read.query_qualities = None

                # Write the read to the SAM file
                outfile.write(read)
                continue

            # Process each alignment in the mappings
            for mapping in mappings:
                start_pos, end_pos, query_aln, genome_aln, score = mapping

                # Create a new AlignedSegment object
                read = pysam.AlignedSegment()
                read.query_name = record_id
                read.flag = 0  # Single-end alignment
                read.reference_id = outfile.get_tid(reference_name)
                read.reference_start = start_pos - 1  # Convert to 0-based position
                read.mapping_quality = 99  # Default mapping quality
                read.cigarstring = aln_to_cig(query_aln, genome_aln)  # Simplified CIGAR string
                read.query_sequence = query_aln
                read.query_qualities = pysam.qualitystring_to_array("I" * len(query_aln))

                # Write the read to the SAM file
                outfile.write(read)