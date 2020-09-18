#!/usr/bin/env nextflow
/*
========================================================================================
                         BiomeFlow PIPELINE
========================================================================================
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/RNA-seq-Flow
 #### Authors
Marko Melnick <mame5141@colorado.edu>
Patrick Gonzales
Thanks to Margaret Gruca <magr0763@colorado.edu> for commands from the RNA-seq flow pipeline
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        1a. SRA tools -- fasterq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files

    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)

    3. Mapping w/ HISAT2 -- map to genome reference file

*/


def helpMessage() {
    log.info"""
    =========================================
     SteadyFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.       

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --skipBAM                      Skip saving BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Compresses and saves all fastq reads.

    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipRSeQC                    Skip running RSeQC.
        
    Analysis Options:
        --count                        Run RSeQC FPKM count over RefSeq annotated genes.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.plaintext_email = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.rcc = "$baseDir/bin/rcc.py"
jgijson_py = "$baseDir/bin/nb_jgiJSON.py"
dummy_seq_fa = "$baseDir/bin/dummysequences.fa"



sample_table = "$baseDir/sample_table.txt"

params.extract_fastqc_stats = "$baseDir/bin/extract_fastqc_stats.sh"
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")


// Validate inputs
if ( params.unmappedPath ){
unmappedPath = params.unmappedPath
}

// Validate inputs
if ( params.magicblastDB ){
magicblastDB = params.magicblastDB
}

// Validate inputs
if ( params.ntblastDB ){
ntblastDB = params.ntblastDB
}

// Validate inputs
if ( params.kmer_size ){
kmer_size = params.kmer_size
} else {
    kmer_size = 55
}

// Validate inputs
if ( params.or ){
or = params.or
} else {
    or = fr
}

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}


if ( params.mapper_index ){
    mapper_index = params.mapper_index
}

if ( params.bowtie2_index ){
    bowtie2_index = params.bowtie2_index
}




if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



/*
 * Create a channel for input read files
 */

if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_trim = Channel.fromPath(sample_table)
    } else {
        fastq_reads_trim = Channel.fromPath(sample_table)
    }
}

if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs+"*_{R1,R2}.fastq.gz")
                            .map { file -> tuple(file.simpleName, file) }
    } else {
        fastq_reads_qc = Channel
            .fromFilePairs( params.fastqs+"/*_{R1,R2}.fastq.gz", size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching ${params.fastqs}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    }
}


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'BiomeFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.fastqs) summary['Fastqs']   = params.fastqs
if(params.sras) summary['SRAs']       = params.sras
summary['Genome Ref']       = params.genome
summary['Mapper Index']     = params.mapper_index
summary['Bowtie2 Index']    = params.bowtie2_index
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save BAM']         = params.skipBAM ? 'NO' : 'YES'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run Count']        = params.count ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    samtools --version > v_samtools.txt
    bedtools --version > v_bedtools.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


process fastQC {
    validExitStatus 0,1
    tag "$name"
    memory '8 GB'
    maxForks 10
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)     "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)    "qc/fastqc/$filename"
        else if (filename.indexOf("txt") > 0)     "qc/fastqc_stats/$filename"
        else null            
    }
    

    input:
    set val(prefix), file(reads) from fastq_reads_qc


    output:
    file "*.{zip,html}" into fastqc_results
    file "*.fastqc_stats.txt" into fastqc_stats

    script:
    """
    echo ${prefix}
    fastqc $reads
    
    ${params.extract_fastqc_stats} \
        --srr=${prefix} \
        > ${prefix}.fastqc_stats.txt    
    """
}






process trim {
    validExitStatus 0
    tag "$name"
    cpus 20
    maxForks 6
    time '24h'
    memory '20 GB'
    //publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*stats.txt"
    publishDir "${params.outdir}/qc/fastq_count", mode: 'copy', pattern: "*fastq_count.txt"
    publishDir "${params.outdir}/qc/fastqtrim_count", mode: 'copy', pattern: "*fastqtrim_count.txt"    
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }      

    input:
    set val(name), val(condition) from fastq_reads_trim.splitCsv(header: ["name", "condition"],sep:'\t')

    output:
    set val(condition), val(name), file("*.trim.fastq.gz") into trimmed_reads_fastqc, trimmed_reads_for_mapping


    file "*fastq_count.txt" into trim_stats2
    file "*fastqtrim_count.txt" into trim_stats3

    script:

    pathname = params.fastqs+"/"+name

        if (!params.singleEnd) {
        """
        echo ${name}      
        java -Xmx4g -jar /Users/mame5141/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 ${pathname}_R1.fastq.gz ${pathname}_R2.fastq.gz ${name}_R1.trim.fastq.gz ${name}_fwd_U.fq ${name}_R2.trim.fastq.gz ${name}_fwd_U.fq ILLUMINACLIP:/Users/mame5141/programs/Trimmomatic-0.36/adapters.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:36

        gzip -cd ${pathname}_R1.fastq.gz | wc -l > ${name}_R1_fastq_count.txt
        gzip -cd ${name}_R1.trim.fastq.gz | wc -l > ${name}_R1_fastqtrim_count.txt                     
        """
    } else {
        """
        echo ${name}      
        java -Xmx4g -jar /Users/mame5141/programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 20 -phred33 ${pathname}_R1.fastq.gz ${name}_R1.trim.fastq.gz ILLUMINACLIP:/Users/mame5141/programs/Trimmomatic-0.36/adapters.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:36

        gzip -cd ${pathname}_R1.fastq.gz | wc -l > ${name}_R1_fastq_count.txt
        gzip -cd ${name}_R1.trim.fastq.gz | wc -l > ${name}_R1_fastqtrim_count.txt                     
        """        
    }
}




process fastQC_trim {
    validExitStatus 0,1
    tag "$name"
    memory '4 GB'
    maxForks 10
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if (filename.indexOf("zip") > 0)   "qc/fastqc/zips/$filename"
        else if (filename.indexOf("html") > 0)  "qc/fastqc/$filename"

        else null            
    }
    
    when:
    !params.skipFastQC && !params.skipAllQC && !noTrim

    input:
    set val(condition), val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*.{zip,html}" into trimmed_fastqc_results


    script:
    """
    echo ${name}
    
    fastqc ${trimmed_reads}
    """
}

process star {
    tag "$name"
    validExitStatus 0
    cpus 20
    maxForks 6
    memory '100 GB'
    time '24h'

    publishDir "${params.outdir}/qc/star_mapstats", mode: 'copy', pattern: "*.final.out"
    publishDir "${unmappedPath}/filter/STAR", mode: 'copy', pattern: "*.fastq"

    input:
    set val(condition), val(name), file(trimmed_reads) from trimmed_reads_for_mapping


    output:
    set val(name), file("*.sam") into mapped1
    set val(condition), val(name), file("*Unmapped.out*") into unmapped1
    file("*.final.out") into star_mapstats


    // Should we disallow soft clip? sp 3,1 to --no-softclip
    //pops out as ${name}Unmapped.out.mate1 or 2
    script:
    if (!params.singleEnd) {
        """
        echo ${name}

        STAR --genomeDir ${genome} \
        --readFilesIn ${name}_R1.trim.fastq.gz ${name}_R2.trim.fastq.gz \
        --readFilesCommand zcat \
        --runThreadN 20 \
        --runMode alignReads \
        --outReadsUnmapped Fastx \
        --outSAMattributes All \
        --outSAMtype SAM \
        --outFileNamePrefix ${name} \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 100

        rm ${name}_R1.trim.fastq.gz 
        rm ${name}_R2.trim.fastq.gz
        """
    } else {
        """
        echo ${name}

        STAR --genomeDir ${genome} \
        --readFilesIn ${name}_R1.trim.fastq.gz \
        --runThreadN 32 \
        --readFilesCommand zcat \
        --runMode alignReads \
        --outReadsUnmapped Fastx \
        --outSAMattributes All \
        --outSAMtype SAM \
        --outFileNamePrefix ${name} \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 100

        rm ${name}_R1.trim.fastq.gz 
        """
    }
}

// Map with bowtie2 again
process bowtie2 {
    tag "$name"
    validExitStatus 0
    cpus 20
    maxForks 6
    memory '40 GB'
    time '24h'

    publishDir "${unmappedPath}/bowtie2", mode: 'copy', pattern: "*.fastq"
    publishDir "${unmappedPath}/final/filter/fastq/initial", mode: 'copy', pattern: "*count.txt"


    input:
    set val(condition), val(name), file(fastqs) from unmapped1

    output:
    set val(condition), val(name), file("*.fastq") into unmapped_bowtie
    set val(name), file("*count.txt") into count1


    print("Running bowtie2 remap")
    script:

    if (!params.singleEnd) {
        """
        bowtie2 \
        -q \
        -p 20 \
        -x ${bowtie2_index} \
        -1 ${name}Unmapped.out.mate1 \
        -2 ${name}Unmapped.out.mate2 \
        -S ${name}_bowtie2.sam \
        --un-conc ${name}_unmapped.fastq \

        wc -l ${name}Unmapped.out.mate1 > ${name}Unmapped_fastq_initial_count.txt

        rm ${name}_bowtie2.sam
        rm ${name}Unmapped.out.mate1
        rm ${name}Unmapped.out.mate2       
        """

    } else {
        """
        bowtie2 \
        -q \
        -p 20 \
        -x ${bowtie2_index} \
        -U ${name}Unmapped.out.mate1 \
        -S ${name}_bowtie2.sam \
        --un ${name}_unmapped.fastq \

        wc -l ${name}Unmapped.out.mate1 > ${name}Unmapped_fastq_initial_count.txt
        mv ${name}_unmapped.fastq ${name}_unmapped.1.fastq

        rm ${name}_bowtie2.sam
        rm ${name}Unmapped.out.mate1
        """
    }


}


process magicBlast {
    tag "$name"
    validExitStatus 0
    cpus 20
    maxForks 6
    memory '100 GB'
    time '24h'

    publishDir "${unmappedPath}/final_unmapped/single", mode: 'copy', pattern: "*{R1,R2}.fastq"
    publishDir "${unmappedPath}/filter/magicblast", mode: 'copy', pattern: "*.tsv"
    publishDir "${unmappedPath}/final/filter/fastq/bowtie2", mode: 'copy', pattern: "*bowtie2_count.txt"
    publishDir "${unmappedPath}/final/filter/fastq/magicblast", mode: 'copy', pattern: "*magicblast_count.txt"

    input:
    set val(condition), val(name), file(fastqs) from unmapped_bowtie

    output:
    set val("single"),  val(name), file("*.fastq") into unmapped_single
    file("*.fastq") into unmapped_single_final_bowtie
    set val(condition), val(name), file("*.fastq") into unmapped_single_copy
    file("*.txt") into dummy_text
    script:

    if (!params.singleEnd) {
        """
        # Grab last 4000 lines of unmapped fastq so assembly wont break after magic blast
        head -n +4000 ${name}_unmapped.1.fastq > ${name}_unmapped.1.fastq_rest
        tail -n -4000 ${name}_unmapped.1.fastq > ${name}_unmapped.1.fastq_small
        head -n +4000 ${name}_unmapped.2.fastq > ${name}_unmapped.2.fastq_rest
        tail -n -4000 ${name}_unmapped.2.fastq > ${name}_unmapped.2.fastq_small      


        magicblast \
        -query ${name}_unmapped.1.fastq_rest \
        -query_mate ${name}_unmapped.2.fastq_rest \
        -db ${magicblastDB} \
        -infmt fastq \
        -outfmt tabular \
        -num_threads 20 \
        -no_unaligned \
        -out ${name}_humanblast.tsv
        
        
        
        wc -l ${name}_unmapped.1.fastq > ${name}Unmapped_fastq_bowtie2_count.txt

        tail -n +4 ${name}_humanblast.tsv > ${name}_notail.tsv
        cut -f 1 ${name}_notail.tsv > ${name}_col1.tsv
        sort ${name}_col1.tsv > ${name}_sorted.tsv
        uniq ${name}_sorted.tsv > ${name}_humanfastqID.txt
        seqkit grep ${name}_unmapped.1.fastq -nv -f ${name}_humanfastqID.txt > ${name}_unmapped_R1.fastq.tmp
        seqkit grep ${name}_unmapped.2.fastq -nv -f ${name}_humanfastqID.txt > ${name}_unmapped_R2.fastq.tmp

        wc -l ${name}_unmapped_R1.fastq.tmp > ${name}Unmapped_fastq_magicblast_count.txt
        
        #Add back in some unmapped reads
        cat ${name}_unmapped_R1.fastq.tmp ${name}_unmapped.1.fastq_small  > ${name}_unmapped_R1.fastq
        cat ${name}_unmapped_R2.fastq.tmp ${name}_unmapped.2.fastq_small  > ${name}_unmapped_R2.fastq

        rm ${name}_humanblast.tsv
        rm ${name}_notail.tsv
        rm ${name}_humanfastqID.txt
        rm ${name}_col1.tsv
        rm ${name}_sorted.tsv

        rm ${name}_unmapped.1.fastq
        rm ${name}_unmapped.1.fastq_rest
        rm ${name}_unmapped.1.fastq_small
        rm ${name}_unmapped_R1.fastq.tmp

        rm ${name}_unmapped.2.fastq
        rm ${name}_unmapped.2.fastq_rest
        rm ${name}_unmapped.2.fastq_small
        rm ${name}_unmapped_R2.fastq.tmp


        """

    } else {
        """
        head -n +4000 ${name}_unmapped.1.fastq > ${name}_unmapped.1.fastq_rest
        tail -n -4000 ${name}_unmapped.1.fastq > ${name}_unmapped.1.fastq_small

        magicblast \
        -query ${name}_unmapped.1.fastq_rest \
        -db ${magicblastDB} \
        -infmt fastq \
        -outfmt tabular \
        -num_threads 20 \
        -no_unaligned \
        -out ${name}_humanblast.tsv

        wc -l ${name}_unmapped.1.fastq > ${name}Unmapped_fastq_bowtie2_count.txt

        tail -n +4 ${name}_humanblast.tsv > ${name}_notail.tsv
        cut -f 1 ${name}_notail.tsv > ${name}_col1.tsv
        sort ${name}_col1.tsv > ${name}_sorted.tsv
        uniq ${name}_sorted.tsv > ${name}_humanfastqID.txt
        seqkit grep ${name}_unmapped.1.fastq -nv -f ${name}_humanfastqID.txt > ${name}_unmapped_R1.fastq.tmp
        echo "dummy" >> ${name}_unmapped_R2.fastq

        wc -l ${name}_unmapped_R1.fastq.tmp > ${name}Unmapped_fastq_magicblast_count.txt

        #Add back in some unmapped reads
        cat ${name}_unmapped_R1.fastq.tmp ${name}_unmapped.1.fastq_small > ${name}_unmapped_R1.fastq





        rm ${name}_humanblast.tsv
        rm ${name}_notail.tsv
        rm ${name}_humanfastqID.txt
        rm ${name}_col1.tsv
        rm ${name}_sorted.tsv

        rm ${name}_unmapped.1.fastq
        rm ${name}_unmapped.1.fastq_rest
        rm ${name}_unmapped.1.fastq_small
        rm ${name}_unmapped_R1.fastq.tmp


        """
    }

}



process concatCondition {
    validExitStatus 0,1
    publishDir "${unmappedPath}/final_unmapped/group", mode: 'copy', pattern: "*.fastq"
    publishDir "${unmappedPath}/final/filter/fastq/concatCondition", mode: 'copy', pattern: "*_count.txt"


    input:
    set val(condition), val(namelist), val(filelist) from unmapped_single_copy.groupTuple()

    output:
    set val("group"), val("group_${condition}"), file("*.fastq") into unmapped_group
    file("*_unmapped_R1.fastq") into unmapped_R1_for_ALL
    file("*_unmapped_R2.fastq") into unmapped_R2_for_ALL
    file("*.txt") into textdummy1

    script:

    if (!params.singleEnd) {    
        files = filelist.toList()
        R1_list = []
        R2_list = []
        for (pair in files) {
        R1_list.add(pair[0])
        R2_list.add(pair[1])
        }

        r1 = R1_list.join(" ")
        r2 = R2_list.join(" ") 
        
        """
        cat ${r1} > group_${condition}_unmapped_R1.fastq &
        cat ${r2} > group_${condition}_unmapped_R2.fastq &

        wc -l group_${condition}_unmapped_R1.fastq > group_${condition}_unmapped_fastq_count.txt
        wait
        """
    } else {
        files = filelist.toList()
        R1_list = []
        for (pair in files) {
            R1_list.add(pair[0])
            }
        r1 = R1_list.join(" ")

        """
        cat ${r1} > group_${condition}_unmapped_R1.fastq
        echo "dummy" >> group_${condition}_unmapped_R2.fastq
        wc -l group_${condition}_unmapped_R1.fastq > group_${condition}_unmapped_fastq_count.txt
        """        


    }
}

process concatALL {
    validExitStatus 0,1

    publishDir "${unmappedPath}/final_unmapped/all", mode: 'copy', pattern: "*.fastq"
    publishDir "${unmappedPath}/final/filter/fastq/concatAll", mode: 'copy', pattern: "*_count.txt"

    input:
    file(R1) from unmapped_R1_for_ALL.collect()
    file(R2) from unmapped_R2_for_ALL.collect()


    output:
    set val("all"), val("all"), file("*.fastq") into unmapped_all
    file("*.txt") into textdummy2
    script:


    if (!params.singleEnd) {    
        """
        cat ${R1} > all_unmapped_R1.fastq &
        cat ${R2} > all_unmapped_R2.fastq &

        wc -l all_unmapped_R1.fastq > all_unmapped_fastq_count.txt
        wait
        """
    } else {
       """
       cat ${R1} > all_unmapped_R1.fastq
       echo "dummy" >> all_unmapped_R2.fastq
       wc -l all_unmapped_R1.fastq > all_unmapped_fastq_count.txt 
       """
    }
}


process spades {
    tag "$name"
    validExitStatus 0
    maxForks 6
    cpus 48
    // queue 'long'
    // memory '500GB'
    // time '100h'
    queue { task.attempt <= 1 ? 'long' : 'highmem' }
    memory { task.attempt <= 1 ? '500GB' : '1000GB' }
    time { task.attempt <= 1 ? '100h' : '300h' }

    maxRetries 2
    errorStrategy 'retry'




    publishDir "${unmappedPath}/spades/${condition}", mode: 'copy', pattern: "*"

    input:
    // set val(condition), val(name), val(filelist) from unmapped_single
    set val(condition), val(name), val(filelist) from unmapped_all.mix(unmapped_group,unmapped_single)

    output:
    set val(condition), val(name), file("*spades.fasta") into spades
    stdout result // capture this because spades prints everything

    script:

    int assembler_memory = 500
    if (task.attempt > 1) {
        assembler_memory = 1000
    }
    // Spades cant find file using {name} from nextflow for some reason. 
    // instead made list of files



    if (!params.singleEnd) {

        R1 = filelist.toList()[0]
        R2 = filelist.toList()[1]

        """
        source /Users/mame5141/.bashrc
        spades.py \
        -t 48 \
        -m ${assembler_memory} \
        --rna \
        -k ${kmer_size} \
        -1 ${R1} \
        -2 ${R2} \
        -o ${name}_spades
        mv ${name}_spades/transcripts.fasta ${name}_spades.fasta 
        source /Users/mame5141/.bashrc2
        """



    } else {
        R1 = filelist.toList()[0]
        """
        source /Users/mame5141/.bashrc
        spades.py \
        -t 48 \
        -m ${assembler_memory} \
        --rna \
        -k ${kmer_size} \
        -s ${R1} \
        -o ${name}_spades
        mv ${name}_spades/transcripts.fasta ${name}_spades.fasta
        source /Users/mame5141/.bashrc2 
        """
        }
}




process blast {   

    tag "$name"
    validExitStatus 0,1
    maxForks 6
    cpus 30
    queue "long"
    time { task.attempt <= 1 ? '90h' : '300h' }

    maxRetries 2
    errorStrategy 'retry'

    publishDir "${unmappedPath}/blast/${condition}", mode: 'copy', pattern: "*.tsv"

    input:
    set val(condition), val(name), file(fa) from spades
    
    output:
    set val(condition), val(name), file(fa), file("*.tsv") into jgi, darkbiome


    script:
    """

    # Add in dummy fa sequences
    cat ${dummy_seq_fa} ${fa} > ${fa}_withdummy

    blastn -db ${ntblastDB} \
    -query ${fa}_withdummy \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -outfmt "6 qseqid sseqid pident evalue staxids sscinames scomnames sskingdoms stitle" \
    -out ${name}_unmapped.tsv \
    -num_threads 30
    """
}




process jgi {
    //echo true
    tag "$name"
    validExitStatus 0,1
    cpus 6
    memory '40 GB'
    time '10h'
    maxForks 8


    publishDir "${unmappedPath}/final/jgi_json/${condition}", mode: 'copy', pattern: "*json"
    publishDir "${unmappedPath}/final/jgi_fasta/${condition}", mode: 'copy', pattern: "*_jgi_filtered.fasta"
    publishDir "${unmappedPath}/final/jgi_df/${condition}", mode: 'copy', pattern: "*_jgi_filtered_df.txt"
    publishDir "${unmappedPath}/final/filter/contigs/blast/initial_contigs/${condition}", mode: 'copy', pattern: "*initial_contigs.txt" 
    publishDir "${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast_9606/${condition}", mode: 'copy', pattern: "*9606.txt" 
    publishDir "${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast_10090/${condition}", mode: 'copy', pattern: "*10090.txt"
    publishDir "${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast_Chordata/${condition}", mode: 'copy', pattern: "*Chordata.txt" 
    publishDir "${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast_Viridiplantae/${condition}", mode: 'copy', pattern: "*Viridiplantae.txt"  
    publishDir "${unmappedPath}/final/filter/contigs/blast/Removed_afterBlast_Artificial/${condition}", mode: 'copy', pattern: "*Artificial.txt"  

    input:
    set val(condition), val(name), file(fa), file(blast_tsv) from jgi
    
    output:
    set val(condition), val(name), file(fa), file("*_jgi.json") into filter
    file("*.txt") into txtout
    file("*_jgi_filtered.fasta") into jgi_fasta

    script:
    """
    python ${jgijson_py} \
    --file ${blast_tsv} \
    --fasta ${fa} 
    """
}




process concatFasta_Bowtie2Index {
    //echo true
    tag "$name"
    validExitStatus 0,1
    cpus 20
    memory '40 GB'
    time '10h'
    maxForks 6

    publishDir "${unmappedPath}/final/fasta/${condition}", mode: 'copy', pattern: "${condition}.fasta"
    publishDir "${unmappedPath}/final/bt2_index/${condition}", mode: 'copy', pattern: "*.bt*"


    input:
    set val(condition), val(name), file(fa), file(json) from filter.groupTuple()
    
    output:
    set val(condition), file("*.bt*"), file("${condition}.fasta") into fasta_bt2
    set val(condition), val(name), file("*.fasta"), file("*.bt*") into fasta_dummy_out

    script:

    fas = fa.toList().join(" ")


    script:
    """
    cat ${fas} > ${condition}.fasta
    bowtie2-build --threads 20 ${condition}.fasta ${condition}
    """
}




process bowtie_and_pileup {
    tag "$condition"
    // errorStrategy 'ignore'
    validExitStatus 0
    cpus 20
    memory '40 GB'
    time '24h'
    maxForks 6

    publishDir "${unmappedPath}/final/bam/${condition}", mode: 'copy', pattern: "*.bam"
    publishDir "${unmappedPath}/final/pileup/${condition}", mode: 'copy', pattern: "*_pileup.txt"

    input:
    set val(condition), file(index), file(fa), val(f), val(f2) from fasta_bt2.combine(unmapped_single_final_bowtie)


    output:
    set val(condition), file("*.bam") into final_bowtie
    file("*_pileup.txt") into pileup_dummy_out

    script:

    fin = f.toString().split("_R1")[0]   // get prefix up to ".


    F = file(f)
    parent = F.parent
    base = F.baseName
    fout = base.toString().split("_R1")[0]   // get prefix up to ".
    print("********\n")
    print(condition)
    print(f)
    print(fout)

    if (!params.singleEnd) {
        """
        bowtie2 -q \
        -p 20 \
        -x ${condition} \
        -1 ${fin}_R1.fastq \
        -2 ${fin}_R2.fastq \
        -S ${fout}_${condition}.sam
        samtools view -@ 20 -bS -o ${fout}_${condition}tmp ${fout}_${condition}.sam
        samtools sort -@ 20 ${fout}_${condition}tmp > ${fout}_${condition}.bam
        samtools index ${fout}_${condition}.bam
        samtools mpileup -f ${fa} ${fout}_${condition}.bam > ${fout}_${condition}_pileup.txt


        """
    } else {
        """
        bowtie2 -q \
        -p 20 \
        -x ${condition} \
        -U ${fin}_R1.fastq \
        -S ${fout}_${condition}.sam
        samtools view -@ 20 -bS -o ${fout}_${condition}tmp ${fout}_${condition}.sam
        samtools sort -@ 20 ${fout}_${condition}tmp > ${fout}_${condition}.bam
        samtools index ${fout}_${condition}.bam
        samtools mpileup -f ${fa} ${fout}_${condition}.bam > ${fout}_${condition}_pileup.txt

        """        
    }
}





///////////////////// Start Dark Biome /////////////////////////////////////////
process darkDust {
    tag "$name"
    validExitStatus 0
    cpus 20
    memory '40 GB'
    time '24h'    
    maxForks 8


    publishDir "${unmappedPath}/final/darkbiome/initial_dust_all/${condition}", mode: 'copy', pattern: "*"
    publishDir "${unmappedPath}/final/darkbiome/initial/${condition}", mode: 'copy', pattern: "*darkbiome_initial.fasta" 
    publishDir "${unmappedPath}/final/darkbiome/afterDust/${condition}", mode: 'copy', pattern: "*afterDust.fasta" 


    input:
    set val(condition), val(name), file(fa), file(blast_tsv) from darkbiome

    output:
    set val(condition), val(name), file("*afterDust.fasta"), file(blast_tsv) into dustout
    set val(condition), val(name), file("*"), file(blast_tsv) into darkbiome_out_all
    script:    

    """
    #extract nodenames for all contigs
    grep "^>" ${fa} > ${name}_nodes.txt
    awk 'sub(/^>/, "")' ${name}_nodes.txt  | sort > ${name}_clean.txt
    cut -f 1 ${blast_tsv} > ${name}_blasthits.txt.tmp
    sort ${name}_blasthits.txt.tmp | uniq > ${name}_blasthits.txt
    comm -23 ${name}_clean.txt ${name}_blasthits.txt > ${name}_no_blast_hits.txt
    seqtk subseq ${fa} ${name}_no_blast_hits.txt > ${name}_darkbiome_initial.fasta


    # Add in dummy fa sequences
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
    """            
}



process darkVector {
    tag "$name"
    validExitStatus 0
    cpus 20
    memory '40 GB'
    time '24h'    
    maxForks 8


    publishDir "${unmappedPath}/final/darkbiome/afterVec/${condition}", mode: 'copy', pattern: "*darkbiome.fasta"
    publishDir "${unmappedPath}/final/darkbiome/afterVec_all/${condition}", mode: 'copy', pattern: "*"

    input:
    set val(condition), val(name), file(fa), file(blast_tsv) from dustout

    output:
    set val(condition), val(name), file("*_darkbiome.fasta"), file(blast_tsv) into darkbiome_out2
    set val(condition), val(name), file("*"), file(blast_tsv) into darkbiome_out_all2
    
    script:
    """
    cat ${dummy_seq_fa} ${fa} > ${fa}_withdummy


    blastn -db ${params.univecDB} \
    -query ${fa}_withdummy \
    -max_target_seqs 1 \
    -outfmt 6 \
    -out ${name}_vec.tab \
    -num_threads 20 \
    -reward 1 \
    -penalty -5 \
    -gapopen 3 \
    -gapextend 3 \
    -dust yes \
    -soft_masking true \
    -evalue 700 \
    -searchsp 1750000000000


    cut -f 1 ${name}_vec.tab > ${name}_vec_name.txt
    grep "^>" ${fa}_withdummy > ${name}_name.txt
    awk 'sub(/^>/, "")' ${name}_name.txt > ${name}_name_clean.txt
    sort ${name}_name_clean.txt > ${name}_name_cleansort.txt
    sort ${name}_vec_name.txt > ${name}_vec_name_sort.txt
    comm -23 ${name}_name_cleansort.txt ${name}_vec_name_sort.txt > ${name}_cleancontigs.txt

    seqtk subseq ${fa}_withdummy ${name}_cleancontigs.txt > ${name}_darkbiome.fasta
    """
}


/*
 * STEP 6 - MultiQC
 */
process multiQC {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    //file ('qc/trimstats/*') from trim_stats.collect()
    //file ('qc/star_mapstats/*') from bam_flagstat.collect()
    //file ('qc/rseqc/*') from rseqc_results.collect()
    //file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('qc/star_mapstats/*') from star_mapstats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[SteadyFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[SteadyFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[Biome] Pipeline Complete"

}
