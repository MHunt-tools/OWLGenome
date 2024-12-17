#!/usr/bin/env/python 

## Import python modules
import argparse
import logging 
import os, sys
from Bio import SeqIO
import numpy as np
import random
from tqdm import tqdm
import itertools
from collections import defaultdict
import pandas as pd
import time

## Add path
modules = [i for i in os.listdir('src') if not i.__contains__('.')]
for module in modules:
    sys.path.insert(1, 'src/' + module)

## Import modules
from utils import *
from index import *
from align import *
from test import *

if __name__ == "__main__":
    # Add arguments
    parser = argparse.ArgumentParser(description='Script to execute OWLGenome read mapper')
    parser.add_argument('-g', '--genome_path', help='Path to genome file', required=True)
    parser.add_argument('-genome_format', '--genome_format', help='Format of genome file (fasta or fastq)', required=False)
    parser.add_argument('-i', '--input_path', help='Path to input (reads) file (fastq)', required=True)
    parser.add_argument('-input_format', '--input_format', help='Format of input (reads) file (fasta or fastq)', required=False)
    parser.add_argument('-o', '--output_path', help='Path to output file (tsv)', required=True)
    parser.add_argument('-k_size', '--k_size', help='Size of kmers for calculation. Default=20', required=False)
    parser.add_argument('-w', '--w', help='Size of window for calculation. Default=3', required=False)
    parser.add_argument('-num_seeds', '--num_seeds', help='Number of Seed Patterns per k_mer Default=3', required=False)
    parser.add_argument('-genome_k_stride', '--genome_k_stride', help='Stride for genome kmer parsing. Default=1', required=False)
    parser.add_argument('-input_k_stride', '--input_k_stride', help='Stride for read kmer parsing. Default=20', required=False)
    parser.add_argument('-ground_truth', '--ground_truth', help='Path to ground truth file for testing', required=False)
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('-l', '--log', action='store_true', help='Enable logging')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite output')
    parser.add_argument('-p', '--precomputed', action='store_true', help='Use if output data is already computed')
    parser.add_argument('-m', '--multiprocessing', action='store_true', help='Use multiprocessing')
    parser.add_argument('-mp', '--max_processes', help='Max number of processes (memory)', required=False)
    parser.add_argument('-minimizer', '--use_minimizers', choices=['none', 'spacedSeeds_window', 'spacedSeeds', 'windowMin'],
                        help='Specify minimizer type: "none", "spacedSeeds_window", "spacedSeeds", or "windowMin" (default)',
                        default='windowMin')
    # Parse arguments
    args = parser.parse_args()

    # Process variables
    genome_path = args.genome_path
    input_path = args.input_path
    output_path = args.output_path
    ground_truth = args.ground_truth
    w = int(args.w) if args.w else 3
    k_size = int(args.k_size) if args.k_size else 23
    num_seeds = int(args.num_seeds) if args.num_seeds else 3
    genome_k_stride = int(args.genome_k_stride) if args.genome_k_stride else 1
    k_stride = int(args.genome_k_stride) if args.genome_k_stride else 1
    input_k_stride = int(args.input_k_stride) if args.input_k_stride else 12
    max_processes = int(args.max_processes) if args.max_processes else 8

    ## Process flags
    verbose = args.verbose if args.verbose else False
    force = args.force if args.force else False
    log = args.log if args.log else False
    precomputed = args.precomputed if args.precomputed else False
    multiprocessing = args.multiprocessing if args.multiprocessing else False
    minimizer = args.use_minimizers 

    ## Check filepaths
    if not os.path.exists(genome_path):
        print('Genome path does not exist. Please provide a valid path to the genome file')
        sys.exit()

    if not os.path.exists(input_path):
        print('Input path does not exist. Please provide a valid path to the input (reads) file')
        sys.exit()

    if os.path.exists(output_path) and not force:
        print('Cowardly refusing to overwrite output path. Use -f to force overwrite')
        sys.exit()

    ## Check file formats
    genome_fmt = args.genome_format if args.genome_format else filepath_to_filetype(genome_path)
    input_fmt = args.input_format if args.input_format else filepath_to_filetype(input_path)

    if not genome_fmt:
        print('Genome file not formatted correctly. Requires fasta or fastq file.')
        sys.exit()

    if not input_fmt:
        print('Input file not formatted correctly. Requires fasta or fastq file.')
        sys.exit()


    if not precomputed:
        ## Get genome sequence (required for match testing)
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        for record in SeqIO.parse(genome_path, genome_fmt):
            genome_seq = str(record.seq.upper())
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for genome sequence extraction: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for genome sequence extraction: {end_cpu_time - start_cpu_time:.2f} seconds")

        ## Make hash table
        if verbose: print('\nMaking genome hash table...')
        start_wall_time = time.time()
        start_cpu_time = time.process_time()

        if minimizer == 'none': hash_table = monitor_function(base_hash_table, log, genome_path, genome_fmt, k_size, genome_k_stride, True)
        elif minimizer == 'spacedSeeds_window': hash_table = monitor_function(spacedSeeds_window, log, genome_path, genome_fmt, k_size, genome_k_stride, num_seeds, w, True)
        elif minimizer == 'spacedSeeds': hash_table = monitor_function(spacedSeeds, log, genome_path, genome_fmt, k_size, genome_k_stride, w, True)
        elif minimizer == 'windowMin': hash_table = monitor_function(windowMin, log, genome_path, genome_fmt, k_size, genome_k_stride, w, True)
        
        end_wall_time = time.time()
        end_cpu_time = time.process_time()

        if log:
            print(f"Wall clock time for hash table construction: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for hash table construction: {end_cpu_time - start_cpu_time:.2f} seconds")

        ## Minimize hash table
        if verbose: print('Minimizing hash table...')
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        hash_table = monitor_function(minimize_hash, log, hash_table)
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for hash table minimization: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for hash table minimization: {end_cpu_time - start_cpu_time:.2f} seconds")
        

        ## Run read mapping with time identifiers
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        if verbose: print('Processing records file...')
        output = monitor_function(process_records_file, log, input_path, input_fmt, 
                                genome_seq, k_size, hash_table, input_k_stride, 
                                use_multiprocessing=multiprocessing, max_processes=max_processes, 
                                )
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        elapsed = end_wall_time - start_wall_time
        if log:
            print(f"Wall clock time for read mapping: {elapsed:.2f} seconds")
            print(f"CPU time for read mapping: {end_cpu_time - start_cpu_time:.2f} seconds\n")


        ## Reshape predicted data
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        calc_df = reshape_process_output_df(output)
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for reshaping predicted data: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for reshaping predicted data: {end_cpu_time - start_cpu_time:.2f} seconds\n")


        ## Write df to file
        if verbose: print('Writing results to tsv...')
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        calc_df.to_csv(output_path, sep='\t')
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for writing output: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for writing output: {end_cpu_time - start_cpu_time:.2f} seconds\n")

        # Call convert_to_sam to write the SAM file
        if log:
            print("Converting results to SAM format...")
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        sam_output_path = output_path.replace('.tsv', '.sam')
        convert_to_sam(
            results=output,
            output_path=sam_output_path,
            reference_name="reference_genome",
            reference_length=len(genome_seq)
        )
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"SAM file successfully created at: {sam_output_path}")
            print(f"Wall clock time for writing output: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for writing output: {end_cpu_time - start_cpu_time:.2f} seconds\n")

        ## Print metrics
        if verbose: print('Done!')
        if log or ground_truth: print(f'\nReads per minute: {(len(calc_df)/elapsed)*60:.2f}')


    if precomputed:
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        calc_df = pd.read_csv(output_path, sep='\t')
        calc_df.pred_start = calc_df.pred_start.astype(float)
        calc_df.pred_end = calc_df.pred_end.astype(float)
        calc_df.set_index('id', inplace=True)
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for loading precomputed data: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for loading precomputed data: {end_cpu_time - start_cpu_time:.2f} seconds")


    if ground_truth:
        ## Load ground truth
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        true_loc_df = process_ground_truth(ground_truth)
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for loading ground truth: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for loading ground truth: {end_cpu_time - start_cpu_time:.2f} seconds")

        ## Calculate precision and recall
        start_wall_time = time.time()
        start_cpu_time = time.process_time()
        true_loc_df.index = true_loc_df.index.astype(str)
        calc_df.index = calc_df.index.astype(str)
        matched_df = pd.concat([true_loc_df, calc_df], axis=1)
        precision, recall = calc_precision_recall(matched_df, distance_cutoff=5)
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        if log:
            print(f"Wall clock time for calculating precision and recall: {end_wall_time - start_wall_time:.2f} seconds")
            print(f"CPU time for calculating precision and recall: {end_cpu_time - start_cpu_time:.2f} seconds")

        ## Print success rate
        print(f'Precision: {precision*100:.2f}%')
        print(f'Recall: {recall*100:.2f}%')
        print('\n')