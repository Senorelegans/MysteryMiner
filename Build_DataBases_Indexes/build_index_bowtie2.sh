#!/bin/bash
#SBATCH --job-name=bowtie2index# Job name
#SBATCH --nodes=1
#SBATCH --ntasks=16 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=4:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=40gb # Memory limit
#SBATCH --output=build_indexes_bowtie2.out
#SBATCH --error=build_indexes_bowtie2.err


########################################################################
################### LOAD NECESSARY MODULES #############################
ml bowtie/2.2.9

IN=/Users/mame5141/genomes/mouse/ensembl
FA=Mus_musculus.GRCm38.dna.primary_assembly.fa


mkdir -p ${IN}/index/bowtie2_index

bowtie2-build --threads 16 ${IN}/${FA} ${IN}/index/bowtie2_index/genome
