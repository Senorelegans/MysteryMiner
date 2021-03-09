# MysteryMiner
MysteryMiner is an RNA-seq pipeline that searches for and quantifies microorganisms (potentially pathogenic or from contamination) between groups of samples (treatment/disease vs control) by utilizing reads that do not map to the host. Mystery Miner was designed using NextFlow with additional downstream python scripts. Below is a high level overview of the pipeline followed by more specific directions of each step below.

## Overview
![Alt text](Figures/fig1.png)

1. QC and trimming.
 Reads are first checked with FastQC and trimmed using Trimmomatic

2. Alignment to host genome.
Reads are aligned using multiple aligners (STAR folowed by Bowtie). First pass mapping is done with STAR to quickly remove the majority of host reads, followed by mapping with Bowtie on a highly sensitive setting. 

3. Assembly and nucleotide blast.
	Non-host (unmapped) reads are assembled into contigs with RNA SPAdes 
	and are identified with BLASTN (regular biome) or remain unknown (dark biome).
	
4. Unknown contigs. Unidentified contigs are filtered for repetitive sequences
 with Dust, filtered by comparing hiearchies of samples (single, group or all) with LAST, 
 and identified with BLASTX (dark biome), or remain unidentified (double dark biome).
	
5. Joint Genome Institute query and Filter. Dark biome and regular biome contigs are assigned a
taxonomy using the Joint Genome Institute (JGI) server and 
filtered to remove further mammalian/host genome contigs

6. Library Normalization and Quantification.
Non-host reads are then mapped to all contigs and 
normalized coverage is calculated for subsequent statistical analysis.


## Files needed to run MysterMiner
	DISCLAIMER: Helper scripts to create indexes and download 
	databases are provided in Build_DataBases_Indexes folder 
	but questions on creating indexes and downloading databases are best 
	answered directly from the websites of the programs as methods can 
	change frequently. 

	Below are all of the files you will need to run MysteryMiner
	
	FASTQ files:
	    must be zipped and renamed with name_R1.fastq.gz for single end
	    or name_R1.fastq.gz and name_R2.fastq.gz for paired end
	    
    Sample Table File:
        Tab seperated file to identify conditions of each sample.
        A template has been provided in sample_table.txt.
        (Headerless file and can handle 2+ conditions)        
        (NOTE: File must be called sample_table.txt and in same folder as main.nf)
	
	Mapping Indexes:
	    Host organism indexes created by STAR and Bowtie. 
		
	Databases:
        nt (nucleotide BLASTN) and nr (non-redundant BLASTX) databases


## Programs (and tested versions) to run MysteryMiner

Programs are imported at the top of each file but listed here for completeness

gcc/7.1.0
fastqc/0.11.8
bowtie/2.2.9
STAR/2.5.2b
ncbi-blast/2.7.1
seqkit/0.9.0

SPAdes-3.13.1-Linux


### Modules and Programs
After creating the databases and preparing the files, 1.0_MysteryMiner_Nextflow.sh is the bash script
that will run the first part of the pipeline.

### Python specific modules


# Part 1. Running NextFlow

After creating the databases and preparing the files, 1.0_MysteryMiner_Nextflow.sh is the bash script
that will run the first part of the pipeline.

Slurm Disclaiimer: MysteryMiner was only tested with Slurm but 
should function similalry on other High Performance Computing Clusters 
with different systems that are able to install NextFlow

NextFlow Disclaimer: Knowledge of NextFlow is not needed, but changing filepaths and some parameters inside
of the main nextflow script (main.nf) is necessary as they would not work properly when used as main input parameters in 
1.0_MysteryMiner_Nextflow.sh. This will be explained as necessary in the steps below.

Additionally, users are encouraged to check each nextflow step in main.nf to set the max_forks (amount of max threads at each step), 
cpus for each step, or swap out trimmers/aligners or add steps so you can design the pipeline to suit your needs.

### 1.0_MysteryMiner_Nextflow.sh

	Input commands:

	Path to your gtf file
	GTF=/gtf/Homo_sapiens.GRCh38.90.gtf

	Paths to folders containing your .bam files and .bedgraph files
	BAMPATH=BAMS_FOLDER
	BEDPATH=BEDS_FOLDER
    Data path variables in BASH
    DATABASE=/RNAseq-Biome-Nextflow
    MAIN=RNAseq-Biome-Nextflow
    
    These don't need to bechanged....
    PROJECT=${MAIN}"/NF_OUT"
    BIN=${MAIN}/bin
    mkdir -p $PROJECT
    
    About each command...    
    nextflow run ${MAIN}/main.nf -resume
            This runs the main.nf script and resumes previous session if it exists
        -profile slurm 
            Selects the slurm profile
        -with-report -with-trace -with-timeline -with-dag flowchart.png
            Reports a directed acyclic graph of nextflow commands
        --email "youremail@gmail.com"
            Put in your email for notifications
        --max_cpus 32
            Max cpus to use in nextflow
        --kmer_size 35
            kmer size to use for spades
        --or rf
            fr or rf for first and reverse strands
        --fastqs ${MAIN}/"fastqs"
            Path to fastq scripts
        --sample_table "sample_table.txt"
            path to sample table. Keep it in the nextflow directory for simplicity
        --singleEnd false
            Choose paired end or single end with true or false
        --STAR_index ${DATABASE}/"ensembl"
            STAR index path
        --bowtie2_index ${DATABASE}/"ensembl/bowtie2_index/genome"
            bowtie2 index path
        --ntblastDB ${DATABASE}/"blastDBall/nt"
            path to nucleotide blast  
        --bin ${BIN}
            path to the bin where scripts are kept
        --workdir ${PROJECT}/temp
            temporary working directory
        --outdir ${PROJECT}
            output directory      
        --unmappedPath ${PROJECT}/"unmapped"
            output for unmapped reads and the majority of output files
    
## Relevant Processes from main.nf
Processes are the commands that run Bash scripts in nextflow.

Since there are main processes for this pipeline. I will show a concise overview of the relevant commands along with file outputs.
Users are encouraged to look in the main.nf to modify commands.
All file output start from the outdir ($PROJECT by default) and will be omitted below for clarity.
unmappedPath (${PROJECT}/unmapped by default) is the location to the majority of the output since we deal with unmapped reads
Simple commands or default commands are also omitted for clarity.
### process fastQC
    Runs fastQC and outputs quality checks into 
	output:
	    Puts quality checked files into qc/fastqc

### process trim
    This is the trim command set up to use Trimmomatic.
    WARNING!!! Set the Java path by hand as this would not work when called from 1.0_MysteryMiner_Nextflow.sh

	commands:
	        java -Xmx4g -jar Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${pathname}_R1.fastq.gz ${pathname}_R2.fastq.gz ${name}_R1.trim.fastq.gz ${name}_fwd_U.fq ${name}_R2.trim.fastq.gz ${name}_fwd_U.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:36

	output:
	    fastq.gz is zipped trimmed file to fastq_trimmed
	    fastq_count.txt is count of reads in FASTQ sent to /qc/fastq_count
	    fastqtrim_count.txt is count of reads in trimmed FASTQ sent to /qc/fastqtrim_count	    
	    
### process star
    Runs the STAR aligner to the host genome. This is a quick first pass to remove the majority of host reads
	commands:
        STAR --genomeDir ${STAR_index}
            STAR index directory
        --readFilesIn ${name}_R1.trim.fastq.gz ${name}_R2.trim.fastq.gz
            File input (Paired end example is shown)
        --readFilesCommand zcat
            Take in a zipped fastq
        --runThreadN 20
            Run with 20 cpus
        --runMode alignReads
            Align the reads
        --outReadsUnmapped Fastx
            Output unmapped fastq files
        --outSAMattributes All
            Output all same attributes
        --outSAMtype SAM
            Output a SAM file instead of BAM
        --outFileNamePrefix ${name}
            Output name preference
        --quantMode GeneCounts
            Run gene counts (Not necessarily needed but easy to keep in)
        --outFilterMultimapNmax 100
            Filter out multimapped reads with over a 100 hits

	output:
        final.out are map stats and put into /qc/star_mapstats
        unmapped FASQTS to the host are put into filter/STAR


### process bowtie2
    This is the 2nd mapping pass after STAR to remove additional host reads
	commands:
        bowtie2 -q
            run bowtie2 in quiet mode
        -p 20
            run with 20 cores
        -x ${bowtie2_index}
            path to bowtie2 index
        -1 ${name}Unmapped.out.mate1
        -2 ${name}Unmapped.out.mate2
            first and 2nd reads
        -S ${name}_bowtie2.sam
            output SAM name
        --un-conc ${name}_unmapped.fastq
            output unmapped reads to a fastq file

	output:
    initial_count.txt is the initial count of fastqs and is in ${unmappedPath}/final/filter/fastq/initial
    the final unmapped fastqs are put into ${unmappedPath}/final_unmapped/single
    bowtie2_count.txt is the count of bowtie reads and is in ${unmappedPath}/final/filter/fastq/bowtie2
	



### process concatCondition and concatAll
    This concatenates samples by group (condition/disease or control) and also by all samples
    So we have three heirarchies single (single samples), group, and All.
    This was suprisingly difficult to do in Nextflow and had to use the Groovy language (Nextflow default language)
    
	commands:
	    Groovy commands 
	    Essentially makes a list of the read files and puts into groups by the condition list
	    ---------------
        files = filelist.toList()
        R1_list = []
        R2_list = []
        for (pair in files) {
        R1_list.add(pair[0])
        R2_list.add(pair[1])
        }

        r1 = R1_list.join(" ")
        r2 = R2_list.join(" ") 	
	    ------------------------------
	
	    Bash commands to concatenate samples
	    cat ${r1} > group_${condition}_unmapped_R1.fastq &
        cat ${r2} > group_${condition}_unmapped_R2.fastq &	        
	

	output:
    concatenated fastq files are put into "${unmappedPath}/final_unmapped/group"
    count.txt files are put into "${unmappedPath}/final/filter/fastq/concatCondition"
	


### process spades
    This is the command that runs spades.
    WARNING!!!! Users should look at this as it is important if you are getting assembly or memory errors.
    When you see ${condition} from this point, it refers to single or group or all in the output
    
	commands:
	
	    Nextflow commands
	    This will attempt to restart spades on a higher memory setting.
	    This is important because the singles don't need as much memory
	    and All the samples grouped together need a lot
	    --------------------------------------------------
	    queue { task.attempt <= 1 ? 'long' : 'highmem' }
	        Which node to select
        memory { task.attempt <= 1 ? '500GB' : '1000GB' }
            Memory to use
        time { task.attempt <= 1 ? '100h' : '300h' }
            Time to use
    
        maxRetries 2
        errorStrategy 'retry'
        --------------------------------------------------

        spades.py
            Run spades
        -t 48
            Run with 48 cores
        -m ${assembler_memory}
            Use memory from above
        --rna
            Run in rna setting
        -k ${kmer_size}
            Run with kmer size from main input
        -1 ${R1}
        -2 ${R2}
            Read names
        -o ${name}_spades
        mv ${name}_spades/transcripts.fasta ${name}_spades.fasta 
            Output spades directory and renaming the file

	output:
	   Output assembled fasta contigs into "${unmappedPath}/spades/${condition}



### blast
    This is the nucleotide blast step
    
	commands:
    blastn -db ${ntblastDB}
        blast to the nt database
    -query ${fa}_withdummy
        file input name
    -max_target_seqs 1
    -max_hsps 1 
        select one max target    
    -outfmt "6 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle"
        incomplete taxonomy output
    -out ${name}_unmapped.tsv
        output file name
    -num_threads 30	
        amount of cpus

	output:
    tsv with blast results and is put into publishDir "${unmappedPath}/blast/${condition}"



### process jgi
    This queries the joint genome institute for a more complete taxonomy
	commands:
	    Runs the python script (nb_jgiJSON.py) that queries jgi.
	    Please look in the jupyter notebook jgiJSON.ipynb if you want to see how each command is run and demo examples.
	    It essentially runs a curl on the server.
	    It also removes human, mice, chordata (after the previous 2), plants, and artifical key terms

	output:
	Outputs jgi contigs, fasta, and dataframes to ${unmappedPath}/final/jgi
	Outputs removed contigs mentioned above to ${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast


### process concatFasta_Bowtie2Index
    This creates bowtie indexes of the filtered contigs files above


### process bowtie_and_pileup 
    This maps reads to the indexes created above and runs pileup to get counts for the contig
	commands:
	    bowtie2 command is similar to the mapping from earlier

        
        samtools view -@ 20 -bS -o ${fout}_${condition}tmp ${fout}_${condition}.sam
            Creates a bam file
        samtools sort -@ 20 ${fout}_${condition}tmp > ${fout}_${condition}.bam
            sorts the bam file
        samtools index ${fout}_${condition}.bam
            indexes the bam file
        samtools mpileup -f ${fa} ${fout}_${condition}.bam > ${fout}_${condition}_pileup.txt
            runs pileup

	output:

        outputs non-host bams to ${unmappedPath}/final/bam/${condition}
        outputs pileup files to ${unmappedPath}/final/pileup/${condition}
    
### process darkDust
    This runs the dust filter (For repetitive sequences) on the dark genome (reads that are not found with nt BLAST)
    
    command:
        
        Extract nodenames for all contigs
            grep "^>" ${fa} > ${name}_nodes.txt
            awk 'sub(/^>/, "")' ${name}_nodes.txt  | sort > ${name}_clean.txt
            cut -f 1 ${blast_tsv} > ${name}_blasthits.txt.tmp
            sort ${name}_blasthits.txt.tmp | uniq > ${name}_blasthits.txt
            comm -23 ${name}_clean.txt ${name}_blasthits.txt > ${name}_no_blast_hits.txt
            seqtk subseq ${fa} ${name}_no_blast_hits.txt > ${name}_darkbiome_initial.fasta
        
    
        Add in dummy fa sequences (needed because it can't run with nothing and removed later)
            cat ${dummy_seq_fa} ${name}_darkbiome_initial.fasta > ${fa}_withdummy
            sdust ${fa}_withdummy > ${name}_dust_mask.bed
            cut -f 1 ${name}_dust_mask.bed > ${name}_dustname.txt
        
            grep "^>" ${fa}_withdummy > ${name}_names.txt
            awk 'sub(/^>/, "")' ${name}_names.txt > ${name}_nameclean.txt
        
            sort ${name}_nameclean.txt > ${name}_namecleans.txt
            sort ${name}_dustname.txt > ${name}_dustnames.txt
        
            comm -23 ${name}_namecleans.txt ${name}_dustnames.txt > ${name}_afterDust.txt
            #extract fasta sequences that match clean contig list.
            seqtk subseq ${fa}_withdummy ${name}_afterDust.txt > ${name}_darkbiome_afterDust.fasta
            cp ${name}_darkbiome_afterDust.fasta ${name}_darkbiome.fasta  
              
	output:
        outputs contigs to publishDir ${unmappedPath}/final/darkbiome/afterDust/${condition} 

# Part 2. 2.0_MysteryMiner_Filter.sh

Nextflow is now complete and we will filter additional contigs using LAST and do the protein blast (BLASTX)
LAST was too difficult to do in Nextflow with the heirarchies which is why we are finishing with python scripts

    Variables for paths: and similar to first script
    MAIN=/RNAseq-Biome-master
    PROJECT=${MAIN}"/NF_OUT"
    BIN=${MAIN}/bin
    BLASTXDB = "RNAseq-Biome-Nextflow/blastDB_all_protein/nr" 

You don't have to change any parameters after this point
The output for this section will be in 
RNAseq-Biome-master/NF_OUT/unmapped/final/darkbiome

### script nb_FilterDarkGenome.py
    Runs LAST and runs BLASTX
    
    This was difficult to do in python and users are encouraged to look in the python file if they are curious as to how it works.
    Essentially it compares the different heirarchies (single,group,all) and retains contigs that are 60% similar to eachother.
    This reduces noise and contigs assembled by error.
    
    It then creates indexes and performs BLASTX on contigs that pass
    
    commands:
        python nb_FilterDarkGenome.py
        --col_data ${MAIN}/sample_table.txt
        --Nextflow_Out ${PROJECT} 
        --cpus 32
        --nr_db ${BLASTXDB}
        
    output:
        
        
### script nb_2.0-PileupNormalize.py
    This script normalizes all of the regular biome and dark biome contigs by library size (From host mapped reads)
    This will scrapt the total mapped (unique and multi) reads
    Create a mapped norm  = (lowest total mapped of all samples) / total mapped
    This then multiplies each contig count by the mapped norm to normalize 
    
### script nb_3.0-PileupDataFrame.py   
    This script aggregates all of the individual pileup files into a dataframe
    
### script nb_4.0-CountContigs.py     
    This script creates dataframes of counts for all of the fastq and contig files from all of the steps above (Part 1 and 2)
    fastq count df is in RNAseq-Biome-master/NF_OUT/unmapped/final/filter/fastq
    regular biome contig count df is in NF_OUT/unmapped/final/filter/contigs/blast
    
### script nb_5FilterTtest.py
    This script runs the T-test and applies FDR correction to p-values
    Comparisons for each perumation of condition are in NF_OUT/FINAL_OUT 




