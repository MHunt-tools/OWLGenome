## Import packages
import os
import sys
import numpy as np
import hashlib
from Bio import SeqIO
from collections import defaultdict

## Add path 
current_dir = os.path.dirname(os.path.abspath(__file__))    # set current directory as absolute path of location
sys.path.insert(1, os.path.abspath(os.path.join(current_dir, '../utils'))) # Adding 'src/utils' directory to sys path

## Import functions from utils
from utils import file_to_idx_kmers


def base_hash_table(file_path, file_fmt, k_size, k_stride, include_revcom):
    '''
    Build a hash table mapping k-mers to respective positions and strand in genome.

    Inputs:
        file_path: Path to the genome fasta file.
        file_fmt: Format of the file ('fasta').
        k_size: Size of k-mers.
        k_stride: Stride for k-mer generation.
        include_revcom: If True, include reverse complement k-mers.

    Returns:
        hash_table: Dictionary where keys are k-mers, and values are lists of tuples (position, strand).
    '''
    # Extract k-mers from the genome file (including strand information)
    genome_record_kmers = file_to_idx_kmers(file_path, file_fmt, k_size, k_stride, include_revcom)

    # Initialize the hash table
    hash_table = defaultdict(list)

    # Populate hash table with k-mer as key and (location, strand) as value
    for loc, kmer, strand in genome_record_kmers[0][1]:
        hash_table[kmer].append((loc, strand))  # Append location and strand

    return hash_table


def spacedSeeds_window(file_path, file_fmt, k_size, k_stride, w, num_seeds, include_revcom):
    '''
    Building a minimized hash table mapping k-mers to respective positions in genome,
    selecting minimizers using the lowest hash value in each sliding window.
    
    This version applies spaced seeds to hash values.

    Inputs:
        file_path: Path to genome fasta
        file_fmt: Format of file ('fasta')
        k_size: Size of k-mers
        k_stride: Stride for k-mer parsing
        include_revcom: If True, includes reverse complement k-mers
        w: Window size for minimizer selection
        num_seeds: Number of spaced seeds to use for k-mer minimization

    Returns:
        Dictionary of minimized k-mers and their positions in the genome.
    '''
    # Function to compute the hash of the k-mer using a spaced seed pattern
    def compute_spaced_seed_hash(kmer, spaced_seed_pattern):
        '''
        Compute a simple hash of the k-mer using the spaced seed pattern.

        Arguments:
        kmer -- The k-mer string (e.g., 'AGCTAG')
        spaced_seed_pattern -- Binary pattern (e.g., '101010') indicating the positions to check

        Returns:
        int -- Simple hash value of the k-mer based on the spaced seed
        '''
        # # Select positions in the k-mer based on the spaced seed pattern
        # selected_kmer = ''.join([kmer[i] for i in range(len(kmer)) if spaced_seed_pattern[i] == '1'])
        
        # # Simple hash: sum of ASCII values (ordinals) of the selected characters
        # return sum(ord(char) for char in selected_kmer)
        return ''.join([kmer[i] for i in range(len(kmer)) if spaced_seed_pattern[i] == '1'])
    
    # Function to generate spaced seed patterns
    def generate_spaced_seeds(k_size, num_seeds):
        '''
        Generate a set of spaced seed patterns.
        Each seed pattern is a binary string with 1s marking positions in the k-mer that will be checked.

        Arguments:
        k_size -- The size of the k-mer
        num_seeds -- Number of spaced seed patterns to generate

        Returns:
        List of spaced seed patterns
        '''
        seeds = []
        for _ in range(num_seeds):
            # Generate a (random/fixed) spaced seed pattern
            pattern = ('11100' * ((k_size + 1) // 2))[:k_size]  # Repeats "10" pattern and slices to k_size
            print(pattern)
            seeds.append(pattern)
        return seeds
    
    # Load k-mers from genome function `file_to_idx_kmers`
    genome_record_kmers = file_to_idx_kmers(file_path, file_fmt, k_size, k_stride, include_revcom)
    
    # Generate spaced seed patterns
    spaced_seeds = generate_spaced_seeds(k_size, num_seeds)

    # Initialize the hash table
    hash_table = defaultdict(list)

    # Process each sequence record
    for record_id, kmer_positions in genome_record_kmers:
        # Collect k-mers and hashes for the sequence using spaced seeds
        kmer_hashes = []
        for loc, kmer, strand in kmer_positions:
            # Apply spaced seed hash calculation for each spaced seed pattern
            hashes = [compute_spaced_seed_hash(kmer, seed) for seed in spaced_seeds]
            kmer_hashes.append((min(hashes), loc, kmer, strand))  # Use the minimum hash from the spaced seeds

        # Slide a window over the k-mer hashes and select minimizers    
        for i in range(len(kmer_hashes) - w + 1):
            window = kmer_hashes[i:i + w]
            min_kmer = min(window, key=lambda x: (x[0], x[1]))  # Select by (hash, position) for smallest & leftmost minimizer in ties
            # Store the minimizer in the hash table
            pos_list = hash_table.get(min_kmer[2], [])
            pos_list.append((min_kmer[1], min_kmer[3]))
            pos_list = list(set(pos_list))
            hash_table[min_kmer[2]] = pos_list   # (position, strand)

    return hash_table


def spacedSeeds(file_path, file_fmt, k_size, k_stride, num_seeds, include_revcom):
    '''
    Building a minimized hash table mapping k-mers to respective positions in genome,
    selecting minimizers using the lowest hash value from spaced seeds for each k-mer.
    
    Inputs:
        file_path: Path to genome fasta
        file_fmt: Format of file ('fasta')
        k_size: Size of k-mers
        k_stride: Stride for k-mer parsing
        include_revcom: If True, includes reverse complement k-mers
        num_seeds: Number of spaced seeds to use for k-mer minimization

    Returns:
        Dictionary of minimized k-mers and their positions in the genome.
    '''
    # Function to compute the hash of the k-mer using a spaced seed pattern
    def compute_spaced_seed_hash(kmer, spaced_seed_pattern):
        '''
        Compute a simple hash of the k-mer using the spaced seed pattern.

        Arguments:
        kmer -- The k-mer string (e.g., 'AGCTAG')
        spaced_seed_pattern -- Binary pattern (e.g., '101010') indicating the positions to check

        Returns:
        int -- Simple hash value of the k-mer based on the spaced seed
        '''
        return ''.join([kmer[i] for i in range(len(kmer)) if spaced_seed_pattern[i] == '1'])
    
    # Function to generate spaced seed patterns
    def generate_spaced_seeds(k_size, num_seeds):
        '''
        Generate a set of spaced seed patterns.
        Each seed pattern is a binary string with 1s marking positions in the k-mer that will be checked.

        Arguments:
        k_size -- The size of the k-mer
        num_seeds -- Number of spaced seed patterns to generate

        Returns:
        List of spaced seed patterns
        '''
        seeds = []
        for _ in range(num_seeds):
            pattern = ('01010' * ((k_size + 1) // 2))[:k_size]  # Repeats "10" pattern and slices to k_size
            seeds.append(pattern)
        return seeds
    
    # Loading k-mers from genome (function `file_to_idx_kmers'
    genome_record_kmers = file_to_idx_kmers(file_path, file_fmt, k_size, k_stride, include_revcom)
    
    # Generate spaced seed patterns
    spaced_seeds = generate_spaced_seeds(k_size, num_seeds)

    # Initialize the hash table
    hash_table = defaultdict(list)

    # Processing each sequence record
    for record_id, kmer_positions in genome_record_kmers:
        # Collect minimized k-mers for the sequence using spaced seeds
        for loc, kmer, strand in kmer_positions:
            # Apply spaced seed hash calculation for each spaced seed pattern
            hashes = [compute_spaced_seed_hash(kmer, seed) for seed in spaced_seeds]
            minimized_hash = min(hashes)  # Use the minimum hash across spaced seeds for this k-mer
            # Store the minimized k-mer in the hash table
            hash_table[kmer].append((loc, strand))  # (position, strand)

    return hash_table


def windowMin(file_path, file_fmt, k_size, k_stride, w, include_revcom):
    '''
    Building a minimized hash table mapping k-mers to respective positions in genome,
    selecting minimizers using the lowest hash value in each sliding window.

    Inputs:
        file_path: Path to genome fasta
        file_fmt: Format of file ('fasta')
        k_size: Size of k-mers
        k_stride: Stride for k-mer parsing
        include_revcom: If True, includes reverse complement k-mers
        w: Window size for minimizer selection

    Returns:
        Dictionary of minimized k-mers and their positions in the genome.

    Note:
	•	Hash: Minimizers are selected based on the hash value, with the smallest hash chosen to represent the window.
	•	Position: 
        For minimizers with identical hash values, the selection prioritizes the leftmost (earliest) position in the sequence.
	•	Strand: 
        If both the hash and position are the same (hash collisions minimal in similar cases), or the forward and reverse complements are locally adjacent, the strand becomes the final differentiating factor. 
        Since Python’s tuple sorting is lexicographic, it treats the strand information as a third “layer” in ordering.
    '''
    genome_record_kmers = file_to_idx_kmers(file_path, file_fmt, k_size, k_stride, include_revcom)
    
    # Initialize the hash table
    hash_table = defaultdict(list)

    # Process each sequence record
    for record_id, kmer_hashes in genome_record_kmers:
        # Slide a window over the k-mer hashes and select minimizers
        for i in range(len(kmer_hashes) - w + 1):
            window = kmer_hashes[i:i + w]
            min_kmer = min(window, key=lambda x: ((x[1]), x[0]))  # Select by (lexicographically smallest hash representing 'w', position)
            # Store the minimizer in the hash table with just loc (position)
            pos_list = hash_table.get(min_kmer[1], [])
            pos_list.append((min_kmer[0], min_kmer[2]))
            pos_list = list(set(pos_list))
            hash_table[min_kmer[1]] = pos_list  # Store loc with None as a placeholder
        # print(hash_table)
    return hash_table


def minimize_hash(hash_table, method='single'):
    '''
    Remove elements from hash table that appear too frequently.
    Methods: 'single' removes all non-singlicate entries.
             'std' removes all entries with > 2 standard deviations above the mean.

    Args:
        hash_table (dict): Original hash table of k-mers with their positions.
        method (str): Minimization method.

    Returns:
        dict: Minimized hash table.
    '''
    # Calculate mean and standard deviation of occurrence counts
    num_std = 1
    item_lengths = [len(positions) for positions in hash_table.values()]
    max_freq = np.round(np.mean(item_lengths) + (num_std * np.std(item_lengths)), 0)
    if method == 'single':
        max_freq = 1
    
    # Remove elements above the frequency threshold
    for key in list(hash_table.keys()):
        if len(hash_table[key]) > max_freq:
            del hash_table[key]

    return hash_table
