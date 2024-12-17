# strobmer.py
import psutil
import time
from threading import Thread
from memory_profiler import memory_usage

def monitor_process(interval=1):
    """
    Monitors the current process's CPU and memory usage at specified intervals.
    """
    process = psutil.Process()
    while True:
        cpu = process.cpu_percent(interval=0.1)
        mem = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MiB
        print(f"[Monitor] CPU Usage: {cpu}% | Memory Usage: {mem:.2f} MiB")
        time.sleep(interval)

# Starting the monitoring in a separate thread
monitor_thread = Thread(target=monitor_process, args=(5,), daemon=True)
monitor_thread.start()

import hashlib
from Bio import SeqIO
from collections import defaultdict
from typing import List, Tuple, Dict

def parse_sequences(file_path: str, file_format: str = "fasta") -> List[str]:
    """
    Parse sequences from a FASTA or FASTQ file.

    Args:
        file_path (str): Path to the input file.
        file_format (str): Format of the input file ('fasta' or 'fastq').

    Returns:
        List[str]: List of DNA sequences extracted from the file.
    """
    sequences = []  # Initialize an empty list to store sequences
    for record in SeqIO.parse(file_path, file_format):
        sequences.append(str(record.seq).upper())  # Convert sequences to uppercase strings
    return sequences  # Return the list of sequences

def generate_kmers(sequence: str, k: int) -> List[Tuple[str, int]]:
    """
    Generate k-mers and their frequencies from a DNA sequence.

    Args:
        sequence (str): The DNA sequence.
        k (int): Length of each k-mer.

    Returns:
        List[Tuple[str, int]]: List of tuples containing k-mers and their counts.
    """
    kmer_freq = defaultdict(int)  # Initialize a dictionary to count k-mer frequencies
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_freq[kmer] += 1  # Increment the count for this k-mer
    return list(kmer_freq.items())  # Convert to list of tuples

def hash_kmer(kmer: str, seed: int = 0) -> int:
    """
    Hash a k-mer using MD5 and return an integer hash value.

    Args:
        kmer (str): The k-mer to hash.
        seed (int): Seed value to vary the hash function.

    Returns:
        int: Integer hash value of the k-mer.
    """
    # Combine k-mer with seed to create different hash functions
    combined = f"{kmer}_{seed}"
    hash_object = hashlib.md5(combined.encode('utf-8'))
    # Convert hexadecimal digest to integer
    return int(hash_object.hexdigest(), 16)

def weighted_minhash(kmers: List[Tuple[str, int]], num_hashes: int = 100) -> List[int]:
    """
    Generate a Weighted MinHash sketch for a list of k-mers with associated weights.

    Args:
        kmers (List[Tuple[str, int]]): List of tuples containing k-mers and their weights (frequencies).
        num_hashes (int): Number of independent hash functions to use.

    Returns:
        List[int]: MinHash sketch represented as a list of minimum hash values.
    """
    sketch = []  # Initialize the MinHash sketch list
    for hash_idx in range(num_hashes):
        min_hash = float('inf')  # Initialize minimum hash value for this hash function
        for kmer, weight in kmers:
            # Generate hash value for k-mer with current hash function (seed)
            hash_val = hash_kmer(kmer, seed=hash_idx)
            # Incorporate weight by dividing hash value (higher weight reduces hash_val)
            weighted_hash = hash_val / weight
            if weighted_hash < min_hash:
                min_hash = weighted_hash  # Update minimum hash value
        sketch.append(int(min_hash))  # Append the minimum hash value to the sketch
    return sketch

def minhash_similarity(sketch1: List[int], sketch2: List[int]) -> float:
    """
    Estimate Jaccard similarity between two Weighted MinHash sketches.
    Args:
        sketch1 (List[int]): First MinHash sketch.
        sketch2 (List[int]): Second MinHash sketch.
    Returns:
        float: Estimated Jaccard similarity.
    """
    assert len(sketch1) == len(sketch2), "Sketches must be of the same length."
    match_count = sum(1 for a, b in zip(sketch1, sketch2) if a == b)  # Count matching hash values
    return match_count / len(sketch1)  # Calculate similarity

# Example Usage (to be implemented in separate file)
if __name__ == "__main__":
    # Parameters
    k = 3  # Length of k-mers
    num_hashes = 103  # Number of hash functions for MinHash

    # Generate a longer sequence by repeating a pattern
    seq0 = "data/Datasets/test_genome.fasta"
    seq1 = "data/Datasets/test_seqs.fasta"

    # Generate k-mers with their frequencies for both sequences
    kmers1 = generate_kmers(seq0, k)
    kmers2 = generate_kmers(seq1, k)

    #print(f"Number of Genome Strobemers: {len(kmers1)}")
    #print(f"Number of Read Strobemers: {len(kmers2)}")

    # Generate Weighted MinHash sketches for both sequences
    sketch1 = weighted_minhash(kmers1, num_hashes=num_hashes)
    sketch2 = weighted_minhash(kmers2, num_hashes=num_hashes)

    # Compute Jaccard similarity between the two sketches
    similarity = minhash_similarity(sketch1, sketch2)
    print(f"Weighted MinHash Est. Jaccard Similarity (Ʌ/V): {similarity}")
    print(sketch1, sketch2)