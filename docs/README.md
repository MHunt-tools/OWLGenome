# OWLGenome: A Genome-Scale Long-Read Mapper
Team members: Dominique Duliepre, Maxwell Hunt, Rachel Ong, RJ Jiao


## Table of Contents
- [Introduction](https://github.com/MHunt-tools/OWLGenome/blob/main/docs/README.md#introduction)
- [Features](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#features)
- [Repository structure](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#repository-structure)
- [Getting started](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#getting-started)
- [Cloning the repository](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#cloning-the-repository)
- [Environmental Setup](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#environment-setup)
- [Usage](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#usage)
- [Understanding Inputs](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#understanding-inputs)
- [Understanding Outputs](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#understanding-outputs)
- [Architectural Design](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#architectural-design)
- [Modules](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#modules)
- [Limitations/Known Issues](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#limitationsknown-issues)
- [Future Enhancements](https://github.com/MHunt-tools/OWLGenome/edit/main/docs/README.md#future-enhancements_=)


## Introduction

In the ever-evolving field of genome sequencing, new and efficient aligners must be developed to match the pace of advancing sequencing technologies. Long-read sequencing offers a unique set of benefits and challenges compared to their short-read counterparts when performing genome assembly. Long reads can span complex genomic regions, which helps resolve large structural variations and repetitive sequences that short reads struggle with. However, long reads require more advanced mappers capable of handling higher error rates, longer sequences, and complex genomic rearrangements. As such, efficient mapping algorithms are essential to handle this computationally intensive task.

Here, we present OWLigner, a genome-scale long-read mapping tool designed for efficient and accurate alignment of long-read sequences. OWLigner achieves this by utilizing kmer minimization for efficient genome hashing, and leveraging a numba JIT enhanced smith-waterman algorithm with the option for multiprocessing.By integrating these advanced methods, OWLigner streamlines the mapping process thus advancing our ability to interpret complex genomes for downstream genomic analyses.


## Features

- Robust Input Configuration: Offers a comprehensive command-line interface for specifying parameters, file paths, and flags.
- Flexible Genome Indexing: Uses advanced k-mer hashing with minimization techniques to create a memory-efficient index of the reference genome, speeding up the mapping process.
- High-performance Read Alignment: Allows fine-tuning of k-mer size, application of window minimizers and spaced seeds, and selection of alignment scoring criteria to meet varied research needs.
- Comprehensive Output and Accuracy Metrics: Converts results into standard formats (e.g. TSV, SAM), provides compatibility with ground truth datasets, and computes precision/recall metrics for performance evaluation.

## Repository Structure

```
├── Midterm Presentation.pptx
├── OWLGenome.png
├── Software Design Document OwlGenome.docx
├── data
│   └── Datasets
│       ├── Long_reads
│       │   ├── long_reads_100_subset.fastq
│       │   ├── long_reads_100_subset_ground_truth.txt
│       │   ├── long_reads_500_subset.fastq
│       │   ├── long_reads_500_subset_ground_truth_1_base.txt
│       │   ├── long_reads_ref_genome.fasta
│       │   └── test_seqs_2_ground_truth.tsv
│       ├── Long_reads.zip
│       ├── test.txt
│       ├── test_genome.fasta
│       ├── test_seqs.fasta
│       ├── test_seqs_2.fasta
│       └── test_seqs_2_ground_truth.tsv
├── docs
│   ├── BasicUserGuide.md
│   └── README.md
├── examples
│   ├── example_script.txt
│   └── output.tsv
├── main.py
├── matched_df.csv
├── output.sam
├── output.tsv
├── requirements.txt
├── results
│   └── output_test.sam
├── scripts
│   └── slurm_script.sh
├── simulate.py
└── src
    ├── __init__.py
    ├── __pycache__
    │   └── __init__.cpython-312.pyc
    ├── alignment
    │   ├── __pycache__
    │   │   ├── align.cpython-312.pyc
    │   │   ├── nw_align.cpython-312.pyc
    │   │   ├── sw_align.cpython-312.pyc
    │   │   ├── sw_align.sw_align_numba-141.py312.1.nbc
    │   │   └── sw_align.sw_align_numba-141.py312.nbi
    │   ├── align.py
    │   ├── nw_align.py
    │   └── sw_align.py
    ├── index
    │   ├── __init__.py
    │   ├── __pycache__
    │   │   └── index.cpython-312.pyc
    │   └── index.py
    ├── sketch
    │   ├── Test.py
    │   ├── scratch.ipynb
    │   ├── strobmer.py
    │   ├── strobmer_Diagrams
    │   │   ├── integrated_strobmer_minhash_flowchart
    │   │   ├── integrated_strobmer_minhash_flowchart.png
    │   │   └── strobmer_FlowChart.ipynb
    │   ├── whtedMinHash.py
    │   └── whtedMinHash_Diagrams
    │       ├── WhtedMinHash_FlowChart
    │       ├── WhtedMinHash_FlowChart.ipynb
    │       └── WhtedMinHash_FlowChart.png
    ├── tests
    │   ├── __init__.py
    │   ├── __pycache__
    │   │   └── test.cpython-312.pyc
    │   └── test.py
    └── utils
        ├── __init__.py
        ├── __pycache__
        │   └── utils.cpython-312.pyc
        └── utils.py

```

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



## Usage

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

### Example Usage

Script:

    python3 main.py -g "data/Datasets/post_midterm/long_reads_ref_genome.fasta" -i "data/Datasets/post_midterm/long_reads_postmidterm.fastq" -o "output.tsv" -v -l -f -ground_truth "data/Datasets/post_midterm/long_reads_postmidterm_ground_truth.txt" -f -m
    
Output:

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



## Architectural Design

The OwlGenome long read mapping software is divided into four main modules. The Input Module (M1) provides a flexible command-line interface for users to specify input files, parameters and options. Next, the Genome Indexing and Hashing Module (M2) builds an efficient data structure to index the reference genome for rapid read kmer lookup. This is followed by the Read Alignment Module (M3), which maps the reads to the genome via the previously-created index and local alignments. Finally, the Output Processing and Ground Truth Comparison Module (M4) produces a user-friendly output file and performance metrics.


<img width="484" alt="Screenshot 2024-12-06 at 9 22 09 PM" src="https://github.com/user-attachments/assets/68d793e3-c319-42c9-8fc5-769f56f6089d">


## Modules

### Input Module
The Input Module is responsible for handling all aspects of user input and initial configuration. It provides a command-line interface for users to specify parameters, file paths, and options that control the behavior of the read alignment software. This module also ensures that all inputs are valid, correctly formatted and accessible, thus setting the stage for successful execution of the following modules.

- Command Line Parsing: Reads and interprets command-line arguments provided by the user, including required inputs and optional parameters.
- Variable Initialization: Assigns user-provided values to internal variables that will be used throughout the software.
- Flag Processing: Detects and sets boolean flags that alter the program’s execution flow, such as verbose output or multiprocessing options.
- Input Validation:
    - File Path Validation: Checks that the provided genome and read input file paths exist and are accessible.
    - File Format Validation: Ensures that input files are in supported formats (e.g. FASTA or FASTQ) and can be correctly parsed.
- Error Handling: Provides meaningful error messages before exiting the program if validation fails.



### Genome Indexing and Hashing Module
The Genome Indexing and Hashing Module focuses on creating an efficient index of the reference genome to enable rapid and accurate mapping of reads. By generating a hash table of k-mers from the genome, the module facilitates the quick lookup of sequences during the alignment process. It also implements minimization techniques to optimize memory usage without sacrificing performance.

- Genome Sequence Loading: Reads the reference genome sequence from the validated file path.
- K-mer Generation: Processes the genome sequence to generate k-mers of a specified length (k_size), including reverse complements.
- Hash Table Creation: Constructs a hash table that maps k-mer sequences to their positions and strands in the genome.
- Minimization Techniques:
    - Window Minimizers: Selects representative k-mers within sliding windows to reduce redundancy.
    - Spaced Seeds: Applies patterns to k-mers to focus on specific positions, reducing the total number of seeds.
    - Combined Minimizers: Utilizes both window minimizers and spaced seeds for enhanced efficiency.
- Hash Statistical Minimization: Further reduces the hash table size by removing k-mers that occur too frequently, based on statistical or set thresholds.


### Long Read Mapping Module
The Long Read Mapping Module is dedicated to mapping sequencing reads to the reference genome using the index created by the Genome Indexing and Hashing Module. It performs the core computational tasks of identifying read kmer matches in the genome index, chaining matched kmers, and producing full local alignments.

- Read Sequence Loading: Reads sequencing reads from the input file validated by the Input Module.
- Read K-mer Generation: Generates k-mers from the reads to facilitate mapping.
- K-mer Mapping: Uses the genome hash table to find matching positions for read k-mers in the genome.
- ptimal K-mer Chain Identification: Constructs an optimal chain of mapped k-mers with a Hough transform to determine the post probable alignment region.
- Alignment start and end calculation: Uses the optimal chain and numba JIT enhanced local alignment to identify where in the genome the read starts and ends
- Full local alignment: Generates a full local alignment of the read to the genome using the previously calculated start and end indices, and numba JIT enhanced local alignment.
- Alignment Scoring: Calculates alignment scores to assess the quality of alignments.
- Multiprocessing Support: Optionally utilizes multiple processors to parallelize read alignment for improved performance.


### Output Processing and Ground Truth Comparison Module
The Output Processing and Ground Truth Comparison Module is designed to convert raw alignment data into standard, user-accessible formats and assess the quality of the resulting alignments. It ensures that the alignment results are properly structured for downstream analyses, optionally compares these results against known truth data, and provides quantitative metrics to evaluate the accuracy of the mapping process.

- Alignment Result Structuring:
    - Data Reshaping: Converts raw alignment tuples into structured, tabular formats (e.g. pandas DataFrames) to facilitate further data manipulation and analysis.
    - Format conversion: Transforms processed alignments into standard bioinformatics formats like TSV and SAM for broad compatibility.
- Data Output Management: Writes the processed results into user-defined file formats at specified output paths, ensuring that the data is reliably accessible for subsequent workflows.
- Ground Truth Integration (Optional):
    - Reference Data Loading: If provided, reads the ground truth file to obtain known accurate mapping positions for comparison.
    - Alignment comparison: Matches predicted read mapping results against ground truth to enable accuracy assessments.
- Performance Evaluation:
    - Metric Calculation: Computes precision and recall metrics to quantify the quality of the alignment process, using configurable criteria to judge correctness.
    - Reporting: Presents performance statistics and for further analysis and interpretation.





## Limitations/Known Issues
The OWLGenome long read mapper demonstrates notable efficiency and precision in mapping genomic reads, but several limitations have been identified that constrain its performance and scalability. Such challenges are outlined below:

- Inefficient Use of Resource Intensive Python Data Structures: Many of the processes in OWLGenome rely on python data structures that are likely hurting the efficiency of the program. For example, our genome hash table is being stored and accessed via a python dictionary, and sequence data is handled as python strings for the majority of functions, besides the alignment, in the package. Converting these data structures to NumPy arrays and all strings into ASCII representations would likely significantly improve performance.
- Underutilization of GPU Acceleration: The NOTS/NOTS-X node cluster includes GPU nodes, which were not leveraged during this project. GPUs offer substantial parallel processing capabilities that could accelerate computationally intensive steps, such as alignment and minimizer selection. However, it is unclear how the sometimes extreme unexplained performance drop observed on Cluster nodes compared to the Login Node  might interfere with GPU-based acceleration. Fully exploring GPU capabilities remains an area for improvement.
- Multiprocessing Efficiency: The current implementation loads all sequences into memory and passes them to each process, which could be optimized further. By passing the SeqIO.parse() generator directly to multiprocessing, memory usage and process communication overhead could be reduced. This optimization may also enable better scalability for larger datasets.
- Adjustments for Greater Complexity: The OWLGenome mapper is designed for typical genomic datasets but does not yet handle more complex scenarios efficiently:
- Large Insertions and Deletions: The current implementation does not expect or handle large structural variants effectively.
Advanced Chaining Methods: Incorporating techniques such as Minimap2’s scoring function could improve chaining accuracy and performance.
Sophisticated Alignment Methods: Existing methods do not include advanced features such as affine gap penalties, which are crucial for accurate alignment in certain genomic regions.
- Hash Table and Minimizer Optimization: The use of hash-based minimizers in the genome index introduces minor inefficiencies. Minimizers are selected based on lexicographically smallest hash values, and although functional, this process could benefit from optimization to reduce redundant computations and memory usage.


## Future Enhancements:
The OWLGenome long read mapper has successfully hit the silver benchmark, but further advancements over the next month can still be made to improve the alignment speed and accuracy. This section highlights key areas for improvement and future development.
- Enhanced Multiprocessing on NOTS / NOTS-X HPC Nodes: Troubleshooting the observed volatilities when multiprocessing on the NOTS / NOTS-X nodes is a critical priority. Optimizing how sequences are loaded and processed across cores, including direct utilization of SeqIO.parse() generators in multiprocessing workflows, could substantially reduce memory overhead and improve runtime performance. Ensuring compatibility with the unique architecture of HPC nodes will also be essential to fully leverage these resources. Additionally, incorporating multiprocessing into other parts of our pipeline, such as the genome-hash index generation, could also be impactful.
- Leveraging GPU Acceleration: Integrating GPU nodes into the workflow represents an exciting avenue for speed improvements. GPUs excel at parallelizing computationally intensive tasks, such as k-mer hashing, minimizer selection, and alignment calculations. However, testing will be required to verify whether GPU acceleration will be affected by the unexplained RPM drop observed on the NOTS / NOTS-X cluster compared to the login node. Frameworks such as CUDA or libraries like Numba and CuPy could be further explored to accelerate critical components.
- Algorithmic Performance Enhancements: Performance optimization is an ongoing goal, with several key areas identified for improvement:
    - Conversion to Efficient Data Structures: Transitioning all k-mers to hashed representations and replacing lists with NumPy arrays with ASCII representation will allow for faster sorting, searching, and indexing operations.
    - Pipeline Profiling Tools: Incorporating tools like cProfile and line_profiler will enable precise identification of bottlenecks in the code, allowing targeted optimizations that yield the greatest runtime reductions.
    - Implementation of Advanced Chaining Techniques: Exploring more efficient chaining algorithms, such as those used in Minimap2, could enhance alignment speed and reduce computational overhead.
- Incorporating Hyperparameter Optimization: Integrating hyperparameter optimization frameworks, such as GridSearch or Bayesian optimization, will allow the selection of the best parameters for k-mer size, stride, window size, and alignment thresholds. Automating this process ensures the pipeline is tailored to the dataset at hand, maximizing both development speed, and recall and specificity.
- Advanced Alignment Methods: Future versions of the mapper should incorporate advanced alignment techniques to handle complex genomic regions:
  - Affine Gap Penalties: Implementing affine gap penalties will improve alignment accuracy in regions with large insertions or deletions, addressing a significant limitation of the current model.
- Translating software package into C: Rewriting the entire software package into code with less overhead, such as C, will be critical to achieving maximal read alignment performance. 
- Scalability for Larger Genomic Data: As genomic datasets grow in complexity and size, the mapper must evolve to maintain efficiency. Scaling the pipeline to handle larger reference genomes and read datasets will require:
    - Optimizing Disk I/O: Minimizing disk read/write bottlenecks, especially on large cluster environments, by using memory-mapped files or more efficient storage formats.
    - Distributed Computing: Exploring distributed computing frameworks, such as Dask or Apache Spark, to parallelize tasks across multiple nodes in an HPC or cloud environment.
- Future-Proofing the Workflow: Given a month of development time with a new fasta/fastq dataset, the mapper’s development should incorporate flexible design principles to ensure adaptability to emerging technologies.
    - Support for New Sequencing Technologies: The pipeline should be updated to accommodate newer sequencing platforms with different read characteristics, such as ultra-long reads or highly accurate short reads.
    - Modular Design: Developing modular components will allow for easier integration of novel algorithms and tools as they are developed in the field.


