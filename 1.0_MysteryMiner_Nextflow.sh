#!/bin/bash
#SBATCH --job-name=NAME# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=12 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=16:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=20gb # Memory limit
#SBATCH --output=OUT.out
#SBATCH --error=ERR.out


################### SET VARIABLES ######################################

export PATH=~:$PATH
export PATH=~/.local/bin:$PATH
export PATH=~/.local:$PATH


########################################################################
################### LOAD NECESSARY MODULES #############################
module load sra/2.8.0
module load bedtools/2.25.0
module load gcc/7.1.0
module load seqkit/0.9.0
module load fastqc/0.11.8
module load bbmap/38.05
module load igvtools/2.3.75
module load preseq/2.0.3
module load mpich/3.2.1
ml gcc/7.1.0
ml bbmap/38.05
ml magicblast/1.3.0
ml bowtie/2.2.9
ml STAR/2.5.2b
ml ncbi-blast/2.7.1

export PATH=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/biome_tools/seqkitsss/seqkit:$PATH
export PATH=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/biome_tools/SPAdes-3.13.1-Linux/bin:$PATH
export PATH=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/biome_tools/seqtk/:$PATH
export PATH=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/biome_tools/minimap2:$PATH
source /Users/mame5141/.bashrc2


################## PRINT JOB INFO ######################################
DATABASE=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow
MAIN=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow
PROJECT=${MAIN}"/NF_OUT"
BIN=${MAIN}/bin

mkdir -p $PROJECT

printf "Sample ID: $ROOTNAME"
printf "\nDirectory: $PROJECT"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date
printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"


############# Command explanations ####################################
# Be sure to check the nextflow.config for default values!!!
# --kmer_size is kmer_size used in spades
# --or is read orientation
# --genome is star index
# --bowtie2 is bowtie2 index
# -fastqs is path to fastq folder. 
# fastqs should be zipped and renamed R1.fastq.gz for single end, 
# or R1/R2.fastq.gz for paired end
# --singleEnd is set to false if paired end data
#######################################################################
nextflow run ${MAIN}/main.nf -resume \
        -profile slurm \
        -with-report -with-trace -with-timeline -with-dag flowchart.png \
        --bin ${BIN} \
        --workdir ${PROJECT}/temp \
        --outdir ${PROJECT} \
        --email "youremail@gmail.com" \
        --max_cpus 32 \
        --kmer_size 35 \
        --or rf \
        --genome ${DATABASE}/"ensembl" \
        --bowtie2_index ${DATABASE}/"ensembl/bowtie2_index/genome" \
        --ntblastDB ${DATABASE}/"blastDBall/nt" \
        --fastqs ${MAIN}/"fastqs" \
        --sample_table "sample_table.txt" \
        --singleEnd false \
        --unmappedPath ${PROJECT}/"unmapped"

