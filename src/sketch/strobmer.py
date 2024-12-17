# strobmer.py
import psutil
import time
from threading import Thread

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

import os, sys
import hashlib
from typing import List, Tuple
import tracemalloc
tracemalloc.start()
import time

## Add path 
current_dir = os.path.dirname(os.path.abspath(__file__))    # set current directory as absolute path of location
sys.path.insert(1, os.path.abspath(os.path.join(current_dir, '../utils'))) # Adding 'src/utils' directory to sys path (1=lower path priority than script path)

from utils import RevCom

class StrobmerGenerator:
    """
    A class to generate strobemers from a DNA sequence.
    
    Attributes:
        k (int): Size of each k-mer w/in the strobemer
        l (int): Size of the window within which the second k-mer of the strobemer is selected
        s (int): Spacing between the two k-mers in the strobemer
        include_revcom (bool): Whether to include reverse complement strobemers
    """
    
    def __init__(self, k: int = 5, l: int = 10, s: int = 2, include_revcom: bool = True):
        """
        Initializes the StrobmerGenerator with specified parameters.
        
        Args:
            k (int): Size of each k-mer within the strobemer.
            l (int): Size of the window within which the second k-mer is selected.
            s (int): Spacing between the two k-mers in the strobemer.
            include_revcom (bool): Whether to include reverse complement strobemers.
        """
        self.k = k  # Size of each k-mer
        self.l = l  # Window size for selecting the second k-mer
        self.s = s  # Spacing between the two k-mers
        self.include_revcom = include_revcom  # Whether to include reverse complements
    
    def generate_strobmers(self, sequence: str) -> List[str]:
        """
        Generates strobemers from the given DNA sequence.
        
        A strobemer is defined as a pair of k-mers separated by a spacing 's' within a window 'l'.
        
        Args:
            sequence (str): DNA sequence from which to generate strobemers.
        
        Returns:
            List[str]: A list of strobemer strings.
        """
        strobemers = []  # Initialize an empty list to store strobemers
        seq_length = len(sequence)  # Total length of the sequence
        
        # Iterate over the sequence to extract strobemers
        for i in range(seq_length - 2 * self.k - self.s + 1):
            # Extract the first k-mer
            kmer1 = sequence[i:i + self.k]
            
            # Define the window for the second k-mer
            window_start = i + self.k + self.s
            window_end = window_start + self.l
            
            # Ensure the window does not exceed the sequence length
            if window_end + self.k > seq_length:
                window_end = seq_length - self.k
            
            # Extract the window sequence
            window_seq = sequence[window_start:window_end]
            
            # Find the k-mer within the window that has the minimum hash value
            min_hash = float('inf')
            selected_kmer2 = ''
            for j in range(len(window_seq) - self.k + 1):
                kmer_candidate = window_seq[j:j + self.k]
                hash_val = self._hash_kmer(kmer_candidate)
                if hash_val < min_hash:
                    min_hash = hash_val
                    selected_kmer2 = kmer_candidate
            
            # Combine the two k-mers to form the strobemer
            if selected_kmer2:
                strobemer = f"{kmer1}-{selected_kmer2}"
                strobemers.append(strobemer)
                
                # If reverse complements are included, add the reverse complement strobemer
                if self.include_revcom:
                    revcom_kmer1 = RevCom(kmer1)
                    revcom_kmer2 = RevCom(selected_kmer2)
                    revcom_strobemer = f"{revcom_kmer1}-{revcom_kmer2}"
                    strobemers.append(revcom_strobemer)
        
        return strobemers  # Return the list of generated strobemers
    
    def _hash_kmer(self, kmer: str) -> int:
        """
        Generates a hash value for a given k-mer using MD5.
        
        Args:
            kmer (str): The k-mer string to hash.
        
        Returns:
            int: Integer hash value of the k-mer.
        """
        # Create an MD5 hash object from the k-mer string
        hash_object = hashlib.md5(kmer.encode('utf-8'))
        # Convert the hexadecimal digest to an integer
        return int(hash_object.hexdigest(), 16)

def generate_strobmer_hash(strobemers: List[str], num_hashes: int = 100) -> List[int]:
    """
    Generates a hash-based sketch for a list of strobemers.
    
    This function is similar to MinHash but tailored for strobemers.
    
    Args:
        strobemers (List[str]): List of strobemer strings.
        num_hashes (int): Number of hash functions to use.
    
    Returns:
        List[int]: A list representing the strobemer sketch.
    """
    # Initialize the sketch with maximum hash values
    sketch = [float('inf')] * num_hashes
    
    # Define unique seeds for each hash function
    seeds = [1000003 + i for i in range(num_hashes)]
    
    for strobemer in strobemers:
        for i in range(num_hashes):
            seed = seeds[i]
            # Combine strobemer with seed to simulate different hash functions
            combined = f"{strobemer}_{seed}"
            hash_val = int(hashlib.md5(combined.encode('utf-8')).hexdigest(), 16)
            # Update the sketch with the minimum hash value for each hash function
            if hash_val < sketch[i]:
                sketch[i] = hash_val
    
    return sketch  # Return the generated strobemer sketch

def estimate_jaccard_similarity_sketch(sketch1: List[int], sketch2: List[int]) -> float:
    """
    Estimates the Jaccard similarity between two strobemer sketches.
    
    Args:
        sketch1 (List[int]): The first strobemer sketch.
        sketch2 (List[int]): The second strobemer sketch.
    
    Returns:
        float: The estimated Jaccard similarity between the two sketches.
    """
    # Ensure both sketches have the same length
    if len(sketch1) != len(sketch2):
        raise ValueError("Sketches must be of the same length to estimate similarity.")
    
    # Count the number of matching hash values
    match_count = sum(1 for h1, h2 in zip(sketch1, sketch2) if h1 == h2)
    
    # Calculate the Jaccard similarity
    similarity = match_count / len(sketch1)
    return similarity

# Example Usage
if __name__ == "__main__":
    # Example sequences (to be replaced w/ actual genome and read sequences)
    genome_sequence = "data/Datasets/test_genome.fasta"
    read_sequence =   "data/Datasets/test_seqs.fasta"
    
    # Initialize StrobmerGenerator with desired parameters
    strobmer_gen = StrobmerGenerator(k=5, l=10, s=2, include_revcom=True)
    
    # Generate strobemers for genome and read
    genome_strobemers = strobmer_gen.generate_strobmers(genome_sequence)
    read_strobemers = strobmer_gen.generate_strobmers(read_sequence)
    
    print(f"Number of Genome Strobemers: {len(genome_strobemers)}")
    print(f"Number of Read Strobemers: {len(read_strobemers)}")

    # Generate sketches for genome and read
    genome_sketch = generate_strobmer_hash(genome_strobemers, num_hashes=1000)
    read_sketch = generate_strobmer_hash(read_strobemers, num_hashes=1000)

    # Estimate Jaccard similarity between genome and read sketches
    similarity = estimate_jaccard_similarity_sketch(genome_sketch, read_sketch)
    print(f"Strobmer Estimated Jaccard Similarity: {similarity}")