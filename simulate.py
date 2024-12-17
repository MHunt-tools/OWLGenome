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

## Add path (somebody who works with packages should know a better way?)
modules = [i for i in os.listdir('src') if not i.__contains__('.')]
for module in modules:
    sys.path.insert(1, 'src/' + module)

## Import modules
from utils import *
from index import *
from align import *
from test import *

# Add arguments
parser = argparse.ArgumentParser(description='Script to execute OWLGenome read mapper')
parser.add_argument('-g', '--genome_path', help='Path to genome file (fasta)', required=True)
parser.add_argument('-fasta_out', '--fasta_out', help='Path to output file (fasta)', required=True)
parser.add_argument('-table_out', '--table_out', help='Path to output file (tsv)', required=True)
parser.add_argument('-f', '--force', action='store_true', help='Force overwrite output')

# Parse arguments
args = parser.parse_args()

# Process variables
genome_path = args.genome_path
fasta_out = args.fasta_out
table_out = args.table_out

## Process flags
force = args.force if args.force else False

## Check filepaths
if not os.path.exists(genome_path):
    print('Genome path does not exist. Please provide a valid path to the genome file')
    sys.exit()

if os.path.exists(fasta_out) and not force:
    print('Cowardly refusing to overwrite fasta_out. Use -f to force overwrite')
    sys.exit()

if os.path.exists(table_out) and not force:
    print('Cowardly refusing to overwrite table_out. Use -f to force overwrite')
    sys.exit()

## Get genome sequence (required for match testing)
for record in SeqIO.parse(genome_path, 'fasta'):
    genome_seq = str(record.seq.upper())

## Write simulated seqs to file
num_iter = 10000
write_simulated_seqs(num_iter, genome_seq, fasta_out, table_out,
                        min_read_length = 500,
                        max_read_length = 12000,
                        revcom_prob = .5,
                        corruption_prob = 1,
                        error_rate = 0.15,
                        indel_rate = 0.05, 
                        unmapped_seq_rate = 0.01,
)

print('\nSimulated sequences output to file')