#!/bin/bash
#SBATCH --job-name=hisat2index# Job name
#SBATCH --nodes=1
#SBATCH --ntasks=16 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=4:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=40gb # Memory limit
#SBATCH --output=build_indexes_hisat2.out
#SBATCH --error=build_indexes_hisat2.err


########################################################################
################### LOAD NECESSARY MODULES #############################

module load hisat2/2.1.0

IN=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/ensembl
FA=Homo_sapiens.GRCh38.dna.primary_assembly.fa


mkdir -p ${IN}/hisat2_index

hisat2-build -p 16 ${IN}/${FA} ${IN}/hisat2_index/genome

