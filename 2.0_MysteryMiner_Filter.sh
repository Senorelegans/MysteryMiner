#!/bin/bash
#SBATCH --job-name=Filter# Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=youremail@gmail.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=10 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=100:05:00 # Time limit hrs:min:sec
#SBATCH --partition long
#SBATCH --mem=120gb # Memory limit
#SBATCH --output=OUT_filter.out
#SBATCH --error=ERR_filter.out



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

module load preseq/2.0.3

ml gcc/7.1.0
ml magicblast/1.3.0
ml bowtie/2.2.9
ml STAR/2.5.2b
ml ncbi-blast/2.7.1




MAIN=/RNAseq-Biome-master
PROJECT=${MAIN}"/NF_OUT"
BIN=${MAIN}/bin
BLASTXDB = "/scratch/Users/mame5141/2019/RNAseq-Biome-Nextflow/blastDB_all_protein/nr" 


####################################################################### 

python ${BIN}/nb_FilterDarkGenome.py \
--col_data ${MAIN}/sample_table.txt \
--Nextflow_Out ${PROJECT}  \
--cpus 32 \
--nr_db ${BLASTXDB}

# Normalize the pileup data
python ${BIN}/nb_2.0-PileupNormalize.py --Nextflow_Out ${PROJECT}

# Turn the pileup data into a dataframe
python ${BIN}/nb_3.0-PileupDataFrame.py --Nextflow_Out ${PROJECT} --Nextflow_path ${MAIN}

# Count fastq amounts and contigs
python ${BIN}/nb_4.0-CountContigs.py --Nextflow_Out ${PROJECT} --col_data ${MAIN}/sample_table.txt

# Run a generic T-Test for each condition comparison and taxonomy level
python ${BIN}/nb_5FilterTtest.py --Nextflow_Out ${PROJECT} --Nextflow_path ${MAIN}

# Optional notebook for a Custom T-Test in file 6.0-CustomTtest.ipynb
# Open this in a jupyter notebook once the other scripts are finished
# This will allow you to create custom queries and remove taxonomies at different levels
# Or select a particular taxonomy to look at. Instructions are in the jupyter notebook.

