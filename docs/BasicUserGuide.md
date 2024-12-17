# OWLGenome Long Read Mapping Software Basic User Guide

## Table of Contents
- [Introduction](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#introduction)
- [Getting Started](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#getting-started)
- [Cloning the Repository](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#cloning-the-repository)
- [Environment Setup](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#environment-setup)
- [Using the Software](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#using-the-software)
- [Example Usage](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#example-usage)
- [Understanding Inputs](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#understanding-inputs)
- [Understanding Outputs](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#understanding-outputs)
- [Troubleshooting](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/BasicUserGuide.md#troubleshooting)


## Introduction

Welcome to the OWLGenome Long Read Mapping Software User Guide. This guide provides detailed instructions on how to use OWLigner for the efficient and accurate mapping of long-read sequences to reference genomes. OWLGenome Long Read Mapping Software is designed to handle the challenges of long-read sequencing data, such as higher error rates and longer sequence lengths, by implementing advanced algorithms and optimized data structures. 

## Getting Started
Python Version: Ensure you have Python 3.9 or later.

Git (Optional but recommended):
- On Linux (Debian/Ubuntu):
      sudo apt-get update
      sudo apt-get install git
- On MacOS (with Homebrew):
      `brew install git`
- On Windows: Download and run the installer from https://git-scm.com/downloads

## Cloning the Repository
With Git: 
- First, navigate to the directory where you want your code.
- Then, run:
    `git clone https://github.com/MHunt-tools/OWLGenome.git`

Alternative without Git:
- Navigate to the Github page: https://github.com/MHunt-tools/OWLGenome
- Click the green "Code" button.
- Select "Download ZIP".
- Unzip the download file to access the code.

## Environment Setup
Creating a virtual environment:

    python3 -m venv <environment_name>
    
Activating the virtual environment (macOS/Linux):

    source <environment_name>/bin/activate


Activating the virtual environment (macOS/Linux):

    <environment_name>\Scripts\activate
    
Installing dependencies:
    
    pip install -r requirements.txt


## Using the software

Basic command structure:

    python main.py -g <genome_path> -i <input_path> -o <output_path> [options]

Required Arguments
- `-g`, `--genome_path`: Path to the genome file.
- `-i`, `--input_path`: Path to the input reads file (FASTQ or FASTA format).
- `-o`, `--output_path`: Path where the output TSV file will be saved.

Optional Arguments and Flags
- `-genome_format`, `--genome_format`: Format of the genome file (fasta or fastq). If not provided, the format is inferred from the file extension.
- `-input_format`, `--input_format`: Format of the input reads file (fasta or fastq). If not provided, the format is inferred from the file extension.
- `-k_size`, `--k_size`: Size of kmers for calculations. Default is 20.
- `-genome_k_stride`, `--genome_k_stride`: Stride for genome k-mer parsing. Default is 1.
- `-input_k_stride`, `--input_k_stride`: Stride for read k-mer parsing. Default is 5.
- `-ground_truth`, `--ground_truth`: Path to the ground truth file for testing and evaluation.
- `-v`, `--verbose`: Enable verbose output to display detailed processing information.
- `-l`, `--log`: Enable logging of performance metrics and additional details.
- `-f`, `--force`: Force overwrite of the output file if it already exists.
- `-p`, `--precomputed`: Use precomputed data if available. Skips processing steps that have been previously completed.
- `-m`, `--multiprocessing`: Utilize multiprocessing to speed up computation on systems with multiple CPU cores.
- `-w`, `--window`: Window size for minimizer selection. Default is 3.
- `-num_seeds`, `--num_seeds`: Number of seeds used in computations. Default is 3.
- `-mp`, `--max_proceses`: Naximum number of processes to run simultaneously in multiprocessing modes. Default is 8.
- `-minimizer`, `--use_minimizers`: Choice of minimizer usage. Options include: `None`, `spacedSeeds`, `windowMin` and `spacedSeeds_window`. Default is `windowMin`.


## Example Usage: 
Script:
```
python3 main.py -g "data/Datasets/post_midterm/long_reads_ref_genome.fasta" -i "data/Datasets/post_midterm/long_reads_postmidterm.fastq" -o "output.tsv" -v -l -f -ground_truth "data/Datasets/post_midterm/long_reads_postmidterm_ground_truth.txt" -f -m
```

Output:

  ```
  Wall clock time for genome sequence extraction: 0.02 seconds
  CPU time for genome sequence extraction: 0.01 seconds
  
  Making genome hash table...
  Function: windowMin
  Memory before: 127.81 MB, after: 1413.06 MB, used: 1285.25 MB
  Execution time: 6.5988 seconds
  Wall clock time for hash table construction: 6.60 seconds
  CPU time for hash table construction: 6.60 seconds
  
  Minimizing hash table...
  Function: minimize_hash
  Memory before: 1413.06 MB, after: 1418.36 MB, used: 5.30 MB
  Execution time: 0.5122 seconds
  Wall clock time for hash table minimization: 0.51 seconds
  CPU time for hash table minimization: 0.51 seconds
  
  Processing records file...
  Function: process_records_file
  Memory before: 1418.36 MB, after: 2652.58 MB, used: 1234.22 MB
  Execution time: 44.7347 seconds
  Wall clock time for read mapping: 44.73 seconds
  CPU time for read mapping: 31.64 seconds
  
  Wall clock time for reshaping predicted data: 0.60 seconds
  CPU time for reshaping predicted data: 0.60 seconds
  
  Writing results to tsv...
  Wall clock time for writing output: 0.06 seconds
  CPU time for writing output: 0.06 seconds
  
  Converting results to SAM format...
  SAM file successfully created at: output.sam
  Wall clock time for writing output: 64.10 seconds
  CPU time for writing output: 63.86 seconds
  
  Done!
  Reads per minute: 113854.33
  Wall clock time for loading ground truth: 0.02 seconds
  CPU time for loading ground truth: 0.02 seconds
  Wall clock time for calculating precision and recall: 1.14 seconds
  CPU time for calculating precision and recall: 1.14 seconds
  Precision: 93.84%
  Recall: 97.46%
  ```

## Understanding Inputs

Genome File
- Format: FASTA or FASTQ
- Description: The reference genome file against which the reads will be mapped.
- Usage: Provide the path to this file using the `-g` or `--genome_path` argument.

Reads File
- Format: FASTA or FASTQ
- Description: The sequencing reads that need to be mapped to the reference genome.
- Usage: Provide the path to this file using the `-i` or `--input_path` argument.

Ground Truth File (Optional)
- Format: TSV
- Description: Contains the true mapping positions of the reads for evaluation purposes.
- Usage: Provide the path using the `-ground_truth` argument to enable precision and recall calculations.

## Understanding Outputs

Predicted Alignment Results
- Format: SAM
- Description: Contains metadata along with each read's alignment information. To be used with downstream bioinformatics pipelines or tasks.
- Key fields:
  - `QNAME`: The query (read) name.
  - `FLAG:` a bitwise flag indicating alignment properties (e.g. if the read is mapped or unmapped)
  - `RNAME`: The reference sequence name.
  - `POS`: The 1-based leftmost position of the alignment.
  - `MAPQ`: Mapping quality score.
  - `CIGAR`: A compact string describing how the read aligns to the reference.
  - `SEQ` and `QUAL`: The read sequence and corresponding quality scores.

Matched Ground Truth File (If Provided)
- Format: CSV
- Description: Contains both predicted and known mapping positions of each read.
- Key fields:
  - `id`: Read identifier.
  - `start`, `end`: Ground truth start and end positions.
  - `pred_start`, `pred_end`: Predicted start and end positions from the mapper.
 
Performance Metrics File 
- Format: CSV
- Description: Provides a quick quantitative assessement of the mapper's performance.
- Key fields:
  - `Precision`: The proportion of predicted mappings that were accurate.
  - `Recall`: The proportion of true mappable reads that were successfully identified.
  - `Reads per minute (RPM)`: Efficiency metric showing how fast the mapper processed the input data.


## Troubleshooting
- Check that you are using Python 3.9 or later:
    `python --version`
- Verify that all required depencies are installed.
- Review the project's ReadMe or documentation for additional giuidance.
