#!/bin/bash -l

#SBATCH --job-name=map
#SBATCH --time=00:30:00
#SBATCH --nodes=1 --exclusive
#SBATCH --mem=10G
#SBATCH --partition=commons

init=`date`
echo "Starting at $init"

python main.py -g data/Datasets/Long_reads/long_reads_ref_genome.fasta -i data/Datasets/Long_reads/long_reads_500_subset.fastq -o output.tsv --ground_truth data/Datasets/Long_reads/long_reads_500_subset_ground_truth_1_base.txt -v -l -f

wait
finish=`date`
echo "Job finished at $finish"
