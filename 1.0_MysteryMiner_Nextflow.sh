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

################### LOAD NECESSARY MODULES #############################
ml gcc/7.1.0
ml fastqc/0.11.8
ml bowtie/2.2.9
ml STAR/2.5.2b
ml ncbi-blast/2.7.1
ml seqkit/0.9.0

# Add programs to path that can't be module loaded
export PATH=biome_tools/SPAdes-3.13.1-Linux/bin:$PATH


################## DATA PATHS ######################################
DATABASE=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow
MAIN=/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow
PROJECT=${MAIN}"/NF_OUT"
BIN=${MAIN}/bin

mkdir -p $PROJECT


##### NEXTFLOW COMMAND #############################
nextflow run ${MAIN}/main.nf -resume \
        -profile slurm \
        -with-report -with-trace -with-timeline -with-dag flowchart.png \
        --email "youremail@gmail.com" \
        --max_cpus 32 \
        --kmer_size 35 \
        --or rf \
        --fastqs ${MAIN}/"fastqs" \
        --sample_table "sample_table.txt" \
        --singleEnd false \
        --STAR_index ${DATABASE}/"ensembl" \
        --bowtie2_index ${DATABASE}/"ensembl/bowtie2_index/genome" \
        --ntblastDB ${DATABASE}/"blastDBall/nt" \        
        --bin ${BIN} \
        --workdir ${PROJECT}/temp \
        --outdir ${PROJECT} \        
        --unmappedPath ${PROJECT}/"unmapped"

############# Command explanations ####################################
# Be sure to check the nextflow.config for default values!!!
# --kmer_size is kmer_size used in spades
# --or is read orientation
# -fastqs is path to fastq folder. 
# fastqs should be zipped and renamed R1.fastq.gz for single end, 
# or R1/R2.fastq.gz for paired end
# --singleEnd is set to false if paired end data
######################################################################










# ################################################################################################
# # Test with chlamydia dataset
# MAIN=/Users/mame5141/2019/biome/GIT
# PROJECT=${MAIN}"/NF_OUT"
# BIN=${MAIN}/bin

# mkdir -p $PROJECT
# nextflow run ${MAIN}/main.nf -resume \
#         -profile slurm \
#         -with-report -with-trace -with-timeline -with-dag flowchart.png \
#         --bin ${BIN} \
#         --workdir ${PROJECT}/temp \
#         --outdir ${PROJECT} \
#         --email "youremail@gmail.com" \
#         --max_cpus 64 \
#         --kmer_size 55 \
#         --or fr \
#         --genome ${DATABASE}/"ensembl" \
#         --mapper_index ${DATABASE}/"ensembl/hisat2_index/genome" \
#         --bowtie2_index ${DATABASE}/"ensembl/bowtie2_index/genome" \
#         --ntblastDB ${DATABASE}/"blastDBall_old/nt" \
#         --fastqs "/Users/mame5141/2019/biome/chlamydia/RNAseq-Biome-master/fastqs" \
#         --sample_table "sample_table.txt" \
#         --unmappedPath ${PROJECT}/"unmapped"
