#!/bin/bash
#SBATCH --job-name=NAME# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=markomelnicksupercomputer@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=12 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=16:05:00 # Time limit hrs:min:sec
#SBATCH --partition short
#SBATCH --mem=20gb # Memory limit
#SBATCH --output=starOUT.out
#SBATCH --error=starERR.out



########################################################################
################### LOAD NECESSARY MODULES #############################
ml STAR/2.5.2b


IN=/Users/mame5141/genomes/mouse/ensembl
FA=${IN}/Mus_musculus.GRCm38.dna.primary_assembly.fa
GTF=${IN}/Mus_musculus.GRCm38.98.gtf



mkdir -p ${IN}/index/STAR_ensembl_99/


STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeFastaFiles ${FA} \
--genomeDir ${IN}/index/STAR_ensembl_99/ \
--sjdbGTFfile ${GTF} \
--sjdbOverhang 99