#!/usr/bin/env nextflow
nextflow.preview.dsl=2
// initParamsToDefaults()
if (params.help) exit 0, helpMessage()
// PLATFORM = "ILLUMINA"
_THREADS = 32

// templateDir = "${workflow.projectDir}/lib/bash_templates"
// processParams()
/* Process the parameters and set the environemnt */
params.name = 'Layer Lab DNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>

// Check if genome exists in the config file
if (params.genomes && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

stepList = defineStepList()
step = params.step ? params.step.toLowerCase() : ''

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step ${step}, see --help for more information"

toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'

skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

annoList = defineAnnoList()
annotateTools = params.annotate_tools ? params.annotateTools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(annotateTools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsvPath = params.input
if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

ch_input_sample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    switch (step) {
        case 'mapping': ch_input_sample = extractFastq(tsvFile); break
        case 'recalibrate': ch_input_sample = extractRecal(tsvFile); break
        case 'variantcalling': ch_input_sample = extractBam(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}
// inputSample
// .map{[it[0],it[3],it[4],it[5],it[6]]}
// .set{ch_input_5_col}

(genderMap, statusMap, ch_input_sample) = extractInfos(ch_input_sample)

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null
// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa_index = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa_index ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? params.genomes[params.genome].dbsnp ?: null : null
params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.germline_resource = params.genome && 'mutect2' in tools ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.known_indels = params.genome && 'mapping' in step ? params.genomes[params.genome].known_indels ?: null : null
params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null

// Initialize channels based on params
// ch_acLoci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
// ch_acLociGC = params.ac_lociGC && 'ascat' in tools ? Channel.value(file(params.ac_lociGC)) : "null"
// ch_chrDir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
// ch_chrLength = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
// ch_dbsnp = params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"

ch_acLoci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
ch_acLociGC = params.ac_lociGC && 'ascat' in tools ? Channel.value(file(params.ac_lociGC)) : "null"
ch_chrDir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
ch_chrLength = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
ch_dbsnp = params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"
// ch_dbsnp_index = params.dbsnp_index && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp_index)) : "null"
ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
// ch_dict = params.dict ? Channel.value(file(params.dict)) : "null"
// ch_fasta_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
// ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
// ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_pon = params.pon ? Channel.value(file(params.pon)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"
// knownIndels is currently a list of file for smallGRCh37, so transform it in a channel
li_knownIndels = []
if (params.known_indels && ('mapping' in step)) params.known_indels.each { li_knownIndels.add(file(it)) }
li_knownIndelsIndex = []
if (params.known_indels_index && ('mapping' in step)) params.known_indels_index.each { li_knownIndelsIndex.add(file(it)) }
ch_known_indels = params.known_indels && params.genome == 'smallGRCh37' ? Channel.value(li_knownIndels.collect()) : params.known_indels ? Channel.value(file(params.known_indels)) : "null"
// ch_known_indels_index = params.known_indels_index && params.genome == 'smallGRCh37' ? Channel.value(li_knownIndelsIndex.collect()) : params.known_indels_index ? Channel.value(file(params.known_indels_index)) : "null"
// ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : Channel.empty()
// ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : Channel.empty()

printSummary()

/* Check if the fastq needs to be split into multiple files using the Nextflow splitFastQ operator */
// ch_input_pair_reads = Channel.empty()
// if (params.split_fastq){
//         // newly splitfastq are named based on split, so the name is easier to catch
//     ch_input_pair_reads = ch_input_sample
//         // .splitFastq(by: ÷, compress:true, file:"split", flat: true, compress: true, pe:true)
//         .splitFastq(by: 1_000_000, compress:true, file:"split", pe:true)
//         .map {idPatient, idSample, idRun, reads1, reads2 ->
//             // The split fastq read1 is the 4th element (indexed 3) its name is split_3
//             // The split fastq read2's name is split_4
//             // It's followed by which split it's acutally based on the mother fastq file
//             // Index start at 1
//             // Extracting the index to get a new IdRun
//             splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
//             newIdRun = idRun + "_" + splitIndex
//             // Giving the files a new nice name
//             newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
//             newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
//             [idPatient, idSample, newIdRun, reads1, reads2]}
// }
// if (params.split_fastq){
//     ch_input_pair_reads.
//     subscribe{println (it)}
// }
workflow{
    // First check if various indexes are provided, if not, create them
    BuildFastaFai(ch_fasta)
    BuildBWAindexes(ch_fasta)
    BuildDict(ch_fasta)
    BuildDbsnpIndex(ch_dbsnp)
    BuildGermlineResourceIndex(ch_germline_resource)
    BuildKnownIndelsIndex(ch_known_indels)
    BuildPonIndex(ch_pon)

    ch_fasta_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : BuildFastaFai.out
    // ch_bwa_index = params.bwa_index ? Channel.fromPath(params.bwa_index) : BuildBWAindexes.out
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BuildBWAindexes.out
    ch_dict = params.dict ? Channel.value(file(params.dict)) : BuildDict.out
    
    ch_dbsnp_index = params.dbsnp ? \
        params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : BuildDbsnpIndex.out \
        : "null"
    ch_germline_resource_index = params.germline_resource ? \
        params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) :BuildGermlineResourceIndex.out \
        : "null"
    // ch_known_indels_index = params.known_indels ? \
    //     params.known_indels_index ? Channel.value(file(params.known_indels_index)) : BuildKnownIndelsIndex.out.collect() \
    //     : "null"
    ch_known_indels_index = params.known_indels ? \
        params.known_indels_index ? Channel.value(li_knownIndelsIndex.collect()) : BuildKnownIndelsIndex.out.collect() \
        : "null"
    
    // ch_known_indels_index = params.known_indels_index && params.genome == 'smallGRCh37' ? Channel.value(li_knownIndelsIndex.collect()) : params.known_indels_index ? Channel.value(file(params.known_indels_index)) : "null"
    ch_pon_index = params.pon_index ? Channel.value(file(params.pon_index)) : BuildPonIndex.out

    BuildIntervals(ch_fasta_fai)
    ch_intervals = params.no_intervals ? "null" : \
                    params.intervals && !('annotate' in step) ? \
                    Channel.value(file(params.intervals)) : BuildIntervals.out
    CreateIntervalBeds(ch_intervals)
    // CreateIntervalBeds.out.flatten().subscribe{println(it)}
    ch_bed_intervals = sortBedIntervalsByDescendingDuration(
                            CreateIntervalBeds.out.flatten() )
    
    if (params.no_intervals && step != 'annotate') bedIntervals = Channel.from(file("no_intervals.bed"))
    // ch_bed_intervals.subscribe{println(it)}
    
    // ch_input_sample
    // .subscribe{println(it)}
    // FastQCFQ(ch_input_sample)

    MapReads(ch_input_sample, 
            ch_bwa_index, 
            ch_fasta, 
            ch_fasta_fai)

    // STEP 1.5: MERGING BAM FROM MULTIPLE LANES
    (ch_single_bams, ch_multiple_bams) = 
    MapReads.out.bam_mapped.groupTuple(by:[0, 1])
    .branch{
        single: it[2].size() == 1
        multiple: it[2].size() > 1
    }

    // // ch_multiple_bams.subscribe{ println it}
    ch_single_bams = ch_single_bams.map{
        idPatient, idSample, idRun, bam ->
        [idPatient, idSample, bam]
    }
    MergeBamMapped(ch_multiple_bams)
    ch_merged_bams = MergeBamMapped.out.mix(ch_single_bams)
    IndexBamFile(ch_merged_bams)
    MarkDuplicates(ch_merged_bams)
    
    // TestBaseRecalibrator(MarkDuplicates.out.marked_bams.collect().combine(ch_bed_intervals),
    //     ch_dbsnp,
    //     ch_dbsnp_index,
    //     ch_fasta,
    //     ch_dict,
    //     ch_fasta_fai,
    //     ch_known_indels,
    //     ch_known_indels_index
    // )
    BaseRecalibrator(MarkDuplicates.out.marked_bams.combine(ch_bed_intervals),
        ch_dbsnp,
        ch_dbsnp_index,
        ch_fasta,
        ch_dict,
        ch_fasta_fai,
        ch_known_indels,
        ch_known_indels_index
    )
    table_gather_bqsr_reports = 
        !params.no_intervals ? BaseRecalibrator.out.groupTuple(by:[0, 1]) : BaseRecalibrator.out
    GatherBQSRReports(table_gather_bqsr_reports)
    bam_apply_bqsr = MarkDuplicates.out.marked_bams
                    .join(GatherBQSRReports.out.recal_table, by:[0,1])
    bam_apply_bqsr = bam_apply_bqsr.combine(ch_bed_intervals)
    ApplyBQSR(
        bam_apply_bqsr,
        ch_dict,
        ch_fasta,
        ch_fasta_fai
    )
    bam_merge_bam_recal = ApplyBQSR.out.groupTuple(by:[0, 1])
    // When using intervals, merge (and in the same process index bam files)
    MergeBamRecal(bam_merge_bam_recal)
    // When not using intervals, just index the bam coming from ApplyBQSR
    IndexBamRecal(bam_merge_bam_recal)
    bam_recal = MergeBamRecal.out.bam_recal.mix(IndexBamRecal.out.bam_recal)
    bam_recal_qc = MergeBamRecal.out.bam_recal_qc.mix(IndexBamRecal.out.bam_recal_qc)

    SamtoolsStats(bam_recal_qc)
    bam_BamQC = bam_recal_qc.mix(MapReads.out.bam_mapped_BamQC)
    // bam_BamQC
    // .subscribe{ log.info(it)}
    // BamQC(bam_BamQC,
    //         ch_target_bed)
    bam_HaplotypeCaller = bam_recal.combine(ch_bed_intervals)
    
    HaplotypeCaller(bam_HaplotypeCaller,
        ch_dbsnp,
        ch_dbsnp_index,
        ch_dict,
        ch_fasta,
        ch_fasta_fai
    )
    gvcf_HaplotypeCaller = HaplotypeCaller.out.gvcf_HaplotypeCaller.groupTuple(by:[0, 1, 2])
    if (params.no_gvcf) gvcf_HaplotypeCaller.close()
    // else gvcf_HaplotypeCaller = gvcf_HaplotypeCaller.dump(tag:'GVCF HaplotypeCaller')
    GenotypeGVCFs(HaplotypeCaller.out.gvcf_GenotypeGVCFs,
    ch_dbsnp,
        ch_dbsnp_index,
        ch_dict,
        ch_fasta,
        ch_fasta_fai
    )
    vcf_GenotypeGVCFs = GenotypeGVCFs.out.vcf_GenotypeGVCFs.groupTuple(by:[0, 1, 2])
    // vcf_ConcatenateVCFs = mutect2Output.mix(vcfFreeBayes, vcfGenotypeGVCFs, gvcfHaplotypeCaller)
    vcf_ConcatenateVCFs = vcf_GenotypeGVCFs.mix(gvcf_HaplotypeCaller)
    ConcatVCF(vcf_ConcatenateVCFs,
        ch_fasta_fai,
        ch_target_bed)
    
    // StrelkaSingle(IndexBamRecal.out,
    //     ch_fasta,
    //     ch_fasta_fai,
    //     ch_target_bed
    // )

    // MantaSingle(IndexBamRecal.out,
    //     ch_fasta,
    //     ch_fasta_fai,
    //     ch_target_bed
    // )

   
} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}



// }
/******************************************************************************************/
                                /* Helper functions */
/******************************************************************************************/

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}
// Layer Lab ascii art
   def layerLabAscii() {
    // def ascii_str = 
    return
    '''
    | |                         | |         | |    
    | |     __ _ _   _  ___ _ __| |     __ _| |__ 
    | |    / _` | | | |/ _ \\ '__| |    / _` | '_ \\ 
    | |___| (_| | |_| |  __/ |  | |___| (_| | |_) |
    |______\\__,_|\\__, |\\___|_|  |______\\__,_|_.__/ 
                __/ |                            
               |___/  
    '''
    
  }

def layerLabMessage() {
  // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
     return """ ${c_dim}----------------------------------------------------${c_reset}
     ${c_cyan} LAYER LAB DNA Seq ANALYSIS PIPELINE ${c_reset}
     ${c_dim}----------------------------------------------------${c_reset}
     """
  
}

def helpMessage() {
  // Display help message
    log.info layerLabMessage()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --input sample.tsv -profile fiji

    Mandatory arguments:
        --input                     Path to input TSV file on mapping, recalibrate and variantcalling steps
                                    Multiple TSV files can be specified with quotes
                                    Works also with the path to a directory on mapping step with a single germline sample only
                                    Alternatively, path to VCF input file on annotate step
                                    Multiple VCF files can be specified with quotes
        -profile                    Configuration profile to use
                                    Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test and more

    Options:
        --no_gvcf                    No g.vcf output from HaplotypeCaller
        --no_intervals              Disable usage of intervals
        --nucleotides_per_second      To estimate interval size
                                    Default: 1000.0
        --target_bed                 Target BED file for targeted or whole exome sequencing
        --step                      Specify starting step
                                    Available: Mapping, Recalibrate, VariantCalling, Annotate
                                    Default: Mapping
        --tools                     Specify tools to use for variant calling:
                                    Available: ASCAT, ControlFREEC, FreeBayes, HaplotypeCaller
                                    Manta, mpileup, Mutect2, Strelka, TIDDIT
                                    and/or for annotation:
                                    snpEff, VEP, merge
                                    Default: None
        --skip_qc                   Specify which QC tools to skip when running Sarek
                                    Available: all, bamQC, BCFtools, FastQC, MultiQC, samtools, vcftools, versions
                                    Default: None
        --annotate_tools             Specify from which tools Sarek will look for VCF files to annotate, only for step annotate
                                    Available: HaplotypeCaller, Manta, Mutect2, Strelka, TIDDIT
                                    Default: None
                                    Adds the following tools for --tools: DNAseq, DNAscope and TNscope
        --annotation_cache          Enable the use of cache for annotation, to be used with --snpEff_cache and/or --vep_cache
        --pon                       panel-of-normals VCF (bgzipped, indexed). See: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
        --pon_index                 index of pon panel-of-normals VCF

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
        --bwa_index                  bwa indexes
                                    If none provided, will be generated automatically from the fasta reference
        --dbsnp                     dbsnp file
        --dbsnp_index                dbsnp index
                                    If none provided, will be generated automatically if a dbsnp file is provided
        --dict                      dict from the fasta reference
                                    If none provided, will be generated automatically from the fasta reference
        --fasta                     fasta reference
        --fasta_fai                  reference index
                                    If none provided, will be generated automatically from the fasta reference
        --germline_resource          Germline Resource File
        --germline_esource_index     Germline Resource Index
                                    If none provided, will be generated automatically if a germlineResource file is provided
        --intervals                 intervals
                                    If none provided, will be generated automatically from the fasta reference
                                    Use --no_intervals to disable automatic generation
        --known_indels               knownIndels file
        --known_indels_index          knownIndels index
                                    If none provided, will be generated automatically if a knownIndels file is provided
    Other options:
        --outdir                    The output directory where the results will be saved
        --sequencing_center         Name of sequencing center to be displayed in BAM file
        --multiqc_config            Specify a custom config file for MultiQC
        --monochrome_logs           Logs will be without colors
        --email                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --max_multiqc_email_file_size   Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        -name                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
        --awsqueue                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion                 The AWS Region for your AWS Batch job to run on
    """.stripIndent()
    // println(params.genome)
}

/*
================================================================================
                                PROCESSES
================================================================================
*/

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built

process BuildBWAindexes {
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta}.*")

    when: !(params.bwa_index) && params.fasta && 'mapping' in step

    script:
    """
    bwa index ${fasta}
    """
}



process BuildDict {
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta.baseName}.dict")

    when: !(params.dict) && params.fasta && !('annotate' in step)

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}



process BuildFastaFai {
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta}.fai")

    when: !(params.fasta_fai) && params.fasta && !('annotate' in step)

    script:
    """
    samtools faidx ${fasta}
    """
}



process BuildDbsnpIndex {
    tag {dbsnp}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(dbsnp)

    output:
        file("${dbsnp}.tbi")

    when: !(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${dbsnp}
    """
}


process BuildGermlineResourceIndex {
    tag {germlineResource}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(germlineResource)

    output:
        file("${germlineResource}.tbi")

    when: !(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools

    script:
    """
    tabix -p vcf ${germlineResource}
    """
}

process BuildKnownIndelsIndex {
    tag {knownIndels}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels)

    output:
        file("${knownIndels}.tbi")

    when: !(params.known_indels_index) && params.known_indels && 'mapping' in step

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}


process BuildPonIndex {
    tag {pon}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(pon)

    output:
        file("${pon}.tbi")

    when: !(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${pon}
    """
}


/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)
process BuildIntervals {
  tag {fastaFai}

  publishDir params.outdir

  input:
    file(fastaFai)

  output:
    file("${fastaFai.baseName}.bed")

  when: !(params.intervals) && !('annotate' in step) && !(params.no_intervals)

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

process CreateIntervalBeds {
    tag {intervals.fileName}

    input:
        file(intervals)

    output:
        file '*.bed'

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotides_per_second}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}


process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip}")

    when: !('fastqc' in skipQC)
    
    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}



// process BedToIntervalList {
//     echo true

//     publishDir "${OUT_DIR}/interval_lists/" , mode: 'copy', overwrite: false
    
//     input:
//     file(intervals_bed)

//     output:
//     file(out_file)

//     script:
//     out_file = "${intervals_bed.simpleName}.interval_list"
//    """
//    gatk BedToIntervalList \
//       -I=$intervals_bed \
//       -O=$out_file \
//       -SD=$ref_dict
//    """
// }
// process TestMapReads {
//     label 'cpus_max'
//     echo true
//     tag {idPatient + "-" + idRun}
//     // publishDir "${params.outdir}/Bams/${idSample}/${idSample}_${idRun}", mode: params.publish_dir_mode

//     input:
//         tuple idPatient, idSample, idRun, file(inputFile1), file(inputFile2)
//         file(bwaIndex) 
//         file(fasta) 
//         file(fastaFai)
//     output:
//     script:
//     """
//         echo "idSample: $idSample" 
//     """
//     }
// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM
process MapReads {
    label 'cpus_max'

    tag {idPatient + "-" + idRun}
    publishDir "${params.outdir}/Bams/${idSample}/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file(inputFile1), file(inputFile2)
        file(bwaIndex) 
        file(fasta) 
        file(fastaFai)

    output:
        // tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        // tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam")
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam"), emit : bam_mapped
        tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam"), emit : bam_mapped_BamQC

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    """
        ${convertToFastq}
        bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
        ${input} | \
        samtools sort --threads ${task.cpus} -m 2G - > ${idSample}_${idRun}.bam
    """
}

// STEP 1.5: MERGING BAM FROM MULTIPLE LANES
process MergeBamMapped {
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, idRun, file(bam)

    output:
        tuple idPatient, idSample, file("${idSample}.bam")

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

process IndexBamFile {
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file(bam), file("*.bai")

    // when: !params.knownIndels

    script:
    """
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    """
}

// STEP 2: MARKING DUPLICATES

process MarkDuplicates {
    label 'cpus_16'
    // label 'memory_max'
    // cache false
    tag {idPatient + "-" + idSample}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.bam.metrics" && 'markduplicates' in skipQC) null
            else if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicateMarked/${it}"
        }

    input:
        tuple idPatient, idSample, file("${idSample}.bam")

    output:
        tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), emit: 'marked_bams'
        file ("${idSample}.bam.metrics")

    when: params.known_indels

    script:
    markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    // markdup_java_options =  "\"-Xms" + 16 + "g -Xmx" + 32 + "g\""
    // markdup_java_options =  "\"-Xms16g -Xmx32g\""
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 50000 \
        --INPUT ${idSample}.bam \
        --METRICS_FILE ${idSample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${idSample}.md.bam
    """
}


process BaseRecalibrator {
    // label 'cpus_1'
    label 'cpus_2'
    // cache false
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed)
        file(dbsnp)
        file(dbsnpIndex) 
        file(fasta) 
        file(dict)
        file(fastaFai)
        file(knownIndels)
        file(knownIndelsIndex)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.table")
        // set idPatient, idSample into recalTableTSVnoInt

    when: params.known_indels

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    // intervalsOptions = ""
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${idSample}.recal.table \
        --tmp-dir /tmp \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}
// STEP 3.5: MERGING RECALIBRATION TABLES

process GatherBQSRReports {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'
    echo true
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked", mode: params.publish_dir_mode, overwrite: false

    input:
        tuple idPatient, idSample, file(recal)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.table"), emit: recal_table
        // set idPatient, idSample into recalTableTSV

    when: !(params.no_intervals)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}

// STEP 4: RECALIBRATING

process ApplyBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'
    // label 'cpus_32'
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample }

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed)
        file(dict)
        file(fasta)
        file(fastaFai) 

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.bam")

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${idSample}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    """
}

// STEP 4.5: MERGING THE RECALIBRATED BAM FILES

process MergeBamRecal {
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai"), emit: bam_recal
        tuple idPatient, idSample, file("${idSample}.recal.bam"), emit: bam_recal_qc
        // set idPatient, idSample into bamRecalTSV

    when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    """
}
// STEP 4.5': INDEXING THE RECALIBRATED BAM FILES

process IndexBamRecal {
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file("${idSample}.recal.bam")

    output:
        tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai"), emit: bam_recal
        tuple idPatient, idSample, file("${idSample}.recal.bam"), emit: bam_recal_qc
        
    when: params.no_intervals

    script:
    """
    samtools index ${idSample}.recal.bam
    """
}

// STEP 5: QC

process SamtoolsStats {
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Reports/${idSample}/SamToolsStats", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        file ("${bam}.samtools.stats.out")

    when: !('samtools' in skipQC)

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}

process BamQC {
    // label 'memory_max'
    label 'cpus_16'
    cache false

    tag {idPatient + "-" + idSample}

    // publishDir "${params.outdir}/Reports/${idSample}/bamQC", mode: params.publish_dir_mode
    publishDir "${params.outdir}/Reports/${idSample}/bamQC/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam) 
        file(targetBED)

    output:
        file("${bam.baseName}")

    when: !('bamqc' in skipQC)

    script:
    use_bed = params.target_bed ? "-gff ${targetBED}" : ''
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        $use_bed \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
    """
}


/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {
    label 'memory_singleCPU_task_sq'
    label 'cpus_2'
    // label 'memory_max'
    // label 'cpus_max'

    tag {idSample + "-" + intervalBed.baseName}
    // tag {idSample} 
    publishDir "${params.outdir}/VariantCalling/${idSample}/HaplotypeCaller", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed) 
        // tuple idPatient, idSample, file(bam), file(bai) 
        file(dbsnp)
        file(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
        tuple val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_HaplotypeCaller
        tuple idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_GenotypeGVCFs
        // tuple idPatient, idSample, file("${idSample}.g.vcf")

    when: 'haplotypecaller' in tools

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
         -L ${intervalBed} \
        -D ${dbsnp} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}


// STEP GATK HAPLOTYPECALLER.2

process GenotypeGVCFs {
    tag {idSample + "-" + intervalBed.baseName}

    input:
        tuple idPatient, idSample, file(intervalBed), file(gvcf)
        file(dbsnp)
        file(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
    tuple val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf"), emit: vcf_GenotypeGVCFs

    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile -I ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -V ${gvcf} \
        -O ${intervalBed.baseName}_${idSample}.vcf
    """
}
process ConcatVCF {
    label 'cpus_8'

    tag {variantCaller + "-" + idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idPatient, idSample, file(vcFiles)
        file(fastaFai)
        file(targetBED)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        tuple variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi"), emit: vcf_concatenated

    when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF') 
      outputFile = "HaplotypeCaller_${idSample}.g.vcf"
    else if (variantCaller == "Mutect2") 
      outputFile = "unfiltered_${variantCaller}_${idSample}.vcf"
    else 
      outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options}
    """
}

process PrintCh{
    echo true
    
    input:
    file(in_ch)

    """
    println("PrintCh()")
    cat $in_ch
    """
}
// STEP STRELKA.1 - SINGLE MODE

process StrelkaSingle {
    label 'cpus_max'
    label 'memory_max'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi")

    when: 'strelka' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${idSample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${idSample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${idSample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${idSample}_variants.vcf.gz.tbi
    """
}


// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'cpus_max'
    label 'memory_max'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        tuple val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi")

    when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    status = statusMap[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}


process CramToBam{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/align/bams" , mode: 'copy', overwrite: false
    input:
    file(cram)

    output:
    file(out_file)
    // file("${out_file}.bai")

    script:
    sample = "${cram.simpleName}"
    out_file = "${sample}.bam"
    
    script:
   """
   samtools view -T $ref_fasta -b -o $out_file $cram \
   && samtools index $out_file
   """
}

// process SplitIntervals{
//     echo true
//     tag "$sample"
//     // cache false
//     publishDir "${OUT_DIR}/splitted_intervals/" , mode: 'copy', overwrite: false
    
//     input:
//     file(bed_file)

//     output:
//     file("*.bed")

//     script:
//     """
//     total_lines=`cat $bed_file|wc -l`
//     lines_per_split=`expr \$total_lines / 24`
//     split -l \$lines_per_split $bed_file split_ \
//     --numeric-suffixes=1 \
//     --additional-suffix=.bed \
//     -e
//     """
// }

// process BamToCram{
//     echo true
//     tag "$sample"
//     publishDir "${OUT_DIR}/align/crams" , mode: 'copy', overwrite: false
    
//     input:
//     file(bam)

//     output:
//     file(out_file)
//     file("${out_file}.crai")

//     script:
//     sample = "${bam.simpleName}"
//     out_file = "${sample}.cram"
    
//     script:
//     """
//     samtools view -C $bam -T $ref_fasta  -o $out_file  \
//     && samtools index $out_file
//     """
// }

process CollectInsertSizeMetrics{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/qc/picard_stats" , mode: 'copy', overwrite: false
    
    input:
    file(bam)
    // file(bam_index)

    output:
    file("${sample}_insert_size_metrics.txt")
    file("${sample}_insert_size_metrics.pdf")

    script:
    sample = "${bam.simpleName}"
    """
    gatk --java-options -Xmx32G CollectInsertSizeMetrics --VALIDATION_STRINGENCY=LENIENT \
    -I=$in_bam \
    -O=${sample}_insert_size_metrics.txt \
    -H=${sample}_insert_size_metrics.pdf
    """
}

process CollectAlignmentSummaryMetrics{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/qc/picard_stats" , mode: 'copy', overwrite: false
    
    input:
    file(bam)

    output:
    file("${sample}_alignment_metrics.txt")
    
    script:
    sample = "${bam.simpleName}"
    """
    gatk --java-options -Xmx32G CollectAlignmentSummaryMetrics --VALIDATION_STRINGENCY=LENIENT \
    -I=$in_bam \
    -O=${sample}_alignment_metrics.txt \
    -R=$ref_fasta
    """
}

process CollectHsMetrics{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/qc/picard_stats" , mode: 'copy', overwrite: false
    
    input:
    file(bam)

    output:
    file("${sample}_hs_metrics.txt")

    script:
    sample = "${bam.simpleName}"
    script:
    """
    gatk --java-options -Xmx32G CollectHsMetrics --VALIDATION_STRINGENCY=LENIENT \
    -I=$in_bam \
    -O=${sample}_hs_metrics.txt \
    -TI=$target_intervals \
    -BI=$bait_intervals \
    -R=$ref_fasta
    """
}

// process CreateRecalibrationTable{
//     echo true
//     tag "$sample"
//     publishDir "${OUT_DIR}/recal_table/", mode: 'copy', overwrite: false

//     input:
//     file (sample_cram)
    
//     output:
//     file(sample_cram)
//     file("${sample}.recal.table")
    
//     script:
//     sample = sample_cram.baseName
//     out_file= "${sample}.recal.table"
//     if (params.dataType == 'wgs')
//         """
//         gatk BaseRecalibrator \
//                 -I  $sample_cram\
//                 --known-sites $dbsnp \
//                 --known-sites $known_indels \
//                 -O $out_file \
//                 -R $ref_fasta
//         """
//     else if (params.dataType == 'wes')
//         """
//         gatk BaseRecalibrator \
//                 -I  $sample_cram\
//                 --known-sites $dbsnp \
//                 --known-sites $known_indels \
//                 --intervals $target_intervals \
//                 -O $out_file \
//                 -R $ref_fasta
//         """
//     else
//         error "Invalid data type: ${params.dataType} (see --help)"
// }

// process ApplyBQSR {
//     echo true
//     tag "$sample"
//     publishDir "${OUT_DIR}/BQSR/", mode: 'copy', overwrite: false

//     input:
//     file(sample_cram) 
//     file(recalibrated_file)

    
//     output:
//     file("${sample}.bqsr.bam")
//     file("${sample}.bqsr.bai")

//     script:
//     sample = "${sample_cram.baseName}"
//     out_file = "${sample}.bqsr.bam"

//     if (params.dataType == 'wgs')
//         """
//         gatk ApplyBQSR \
//             -I $sample_cram \
//             -R $ref_fasta \
//             -bqsr $recalibrated_file \
//             -O $out_file 
//         """
//     else if (params.dataType == 'wes')
//         """
//         gatk ApplyBQSR \
//             -I $sample_cram \
//             -R $ref_fasta \
//             -bqsr $recalibrated_file \
//             --intervals $target_intervals \
//             -O $out_file
//         """
//     else
//         error "Invalid data type: ${params.dataType} (see --help)"
// }

// process RunHaplotypeCaller {
//     echo true
//     tag "$sample"
//     publishDir "${OUT_DIR}/haplotype_caller/", mode: 'copy', overwrite: false

//     input:
//     // file(interval_list)
//     file(file_bqsr_bam)
//     // file('*')
    
//     output:
//     file "${sample}.gvcf.gz" 
//     // file "${sample}.gvcf.idx" 
//     // file "${sample}.gvcf.gz.tbi"
    
//     script:
//     sample = file_bqsr_bam.simpleName
//     out_file = "${sample}.gvcf"
//     // if (params.dataType == 'wgs')
//     """
//     gatk --java-options "-Xmx18g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
//         HaplotypeCaller \
//             -I $file_bqsr_bam \
//             --dbsnp $dbsnp \
//             -O $out_file \
//             --emit-ref-confidence GVCF \
//             -R $ref_fasta
//     """
// }


// process GenotypeGVCF{
//     echo true
//     publishDir "${OUT_DIR}/genotyped_gvcf/", mode: 'copy', overwrite: false

//     input:
//     file(gvcf)

//     output:
//     file "$out_file"

//     script:
//     out_file = "${gvcf.baseName}.vcf"
//     """
//     tabix $gvcf
//     gatk --java-options -Xmx64g \
//         GenotypeGVCFs \
//         --dbsnp $dbsnp \
//         -R $ref_fasta \
//         -O $out_file  \
//         -V $gvcf
//     """
// }

// process ConcatVCF{
//     echo true
//     publishDir "${OUT_DIR}/results/", mode: 'copy', overwrite: false
    
//     input:
//     file('*')

//     output:
//     file(out_file)
//     file "${out_file}.gz" 
//     file "${out_file}.gz.tbi"

//     script:
//     out_file='cohort_joint.vcf'
//     // vcf_list=''
//     // chrmList.view{
//     //     vcf_list += "$it "
//     // }
//     vcfs = ''
//     listOfChromosoms.each{
//         vcfs += "${it}.vcf "
//     }
//     """
//     bcftools concat -o ${out_file} -Ov $vcfs
//     bgzip -c ${out_file} > ${out_file}.gz
//     tabix -p vcf ${out_file}.gz
// """
// }

process RunCSQ{
    echo true
    publishDir "${OUT_DIR}/csq/", mode: 'copy', overwrite: false

    input:
    file(vcf)

    output:
    file "$out_file"

    script:
    out_file = "cohort.bcf"    
    """
    bcftools norm -m- -Ov  -f $ref_fasta -w 10000 $vcf |\
    bcftools csq -s - --ncsq 40 -g $ensembl_gene_annotation -l -f $ref_fasta -Ob -o $out_file
    """
}

process VariantEval{
    echo true
    publishDir "${OUT_DIR}/variant_eval/", mode: 'copy', overwrite: false

    input:
    file(cohort_vcf)
    file(vcf_index)

    output:
    file "$out_file"

    script:
    out_file = "cohort.eval.grp"    
    script:
    """ 
    gatk VariantEval --eval $cohort_vcf --comp $dbsnp -R $ref_fasta --output $out_file
    """  
}

process RunMultiQC {
    publishDir "${OUT_DIR}/qc/multiqc", mode: 'copy', overwrite: false

    input:
    file (fastqc:'fastqc/*')
    file ('gatk_base_recalibration/*')
    file ('gatk_variant_eval/*')
    
    output:
    file '*multiqc_report.html'
    file '*_data'
    file '.command.err'
    val prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename  . 2>&1
    """
}

/******************************************************************************************/
                                /* Helper functions */
/******************************************************************************************/


def printSummary(){
    /*
    ================================================================================
                                    PRINTING SUMMARY
    ================================================================================
    */

    // Header log info
    log.info layerLabMessage()
    def summary = [:]
    if (workflow.revision)          summary['Pipeline Release']    = workflow.revision
    summary['Run Name']          = custom_runName ?: workflow.runName
    summary['Max Resources']     = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
    if (workflow.containerEngine)   summary['Container']         = "${workflow.containerEngine} - ${workflow.container}"
    if (params.input)               summary['Input']             = params.input
    if (params.target_bed)           summary['Target BED']        = params.target_bed
    if (step)                       summary['Step']              = step
    if (params.tools)               summary['Tools']             = tools.join(', ')
    if (params.skip_qc)              summary['QC tools skip']     = skipQC.join(', ')

    if (params.no_intervals && step != 'annotate') summary['Intervals']         = 'Do not use'
    if ('haplotypecaller' in tools)                summary['GVCF']              = params.noGVCF ? 'No' : 'Yes'
    // if ('strelka' in tools && 'manta' in tools )   summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
    if (params.sequencing_center)                  summary['Sequenced by']      = params.sequencing_center
    if (params.pon && 'mutect2' in tools)          summary['Panel of normals']  = params.pon

    // summary['Save Genome Index'] = params.saveGenomeIndex ? 'Yes' : 'No'
    // summary['Nucleotides/s']     = params.nucleotidesPerSecond
    summary['Output dir']        = params.outdir
    summary['Launch dir']        = workflow.launchDir
    summary['Working dir']       = workflow.workDir
    summary['Script dir']        = workflow.projectDir
    summary['User']              = workflow.userName
    summary['genome']            = params.genome

    if (params.fasta)                 summary['fasta']                 = params.fasta
    if (params.fasta_fai)              summary['fasta_fai']              = params.fasta_fai
    if (params.dict)                  summary['dict']                  = params.dict
    if (params.bwa_index)              summary['bwa_index']              = params.bwa_index
    if (params.germline_resource)      summary['germline_resource']      = params.germline_resource
    if (params.germline_resource_index) summary['germline_resource_index'] = params.germline_resource_index
    if (params.intervals)             summary['intervals']             = params.intervals
    if (params.ac_loci)                summary['ac_loci']                = params.ac_loci
    if (params.ac_loci_gc)              summary['ac_loci_gc']              = params.ac_loci_gc
    if (params.chr_dir)                summary['chr_dir']                = params.chr_dir
    if (params.chr_length)             summary['chr_length']             = params.chr_length
    if (params.dbsnp)                 summary['dbsnp']                 = params.dbsnp
    if (params.dbsnp_index)            summary['dbsnp_index']            = params.dbsnp_index
    if (params.known_indels)           summary['known_indels']           = params.known_indels
    if (params.known_indels_index)      summary['known_indels_index']      = params.known_indels_index


    if (workflow.profile == 'awsbatch') {
        summary['AWS Region']        = params.awsregion
        summary['AWS Queue']         = params.awsqueue
    }
    summary['Config Profile'] = workflow.profile
    if (params.config_profile_description)  summary['Config Description']  = params.config_profile_description
    if (params.config_profile_contact)      summary['Config Contact']      = params.config_profile_contact
    if (params.config_profile_url)          summary['Config URL']          = params.config_profile_url
    if (params.email) {
        summary['E-mail Address']        = params.email
        summary['MultiQC maxsize']       = params.maxMultiqcEmailFileSize
    }
    log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
    if (params.monochrome_logs) log.info "----------------------------------------------------"
    else log.info "\033[2m----------------------------------------------------\033[0m"

}
def sortBedIntervalsByDescendingDuration(bedIntervals){
    bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
            final fields = line.split('\t')
            if (fields.size() >= 5) duration += fields[4].toFloat()
            else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotides_per_second
            }
        }
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}
}

//  def getChrmList(){
//     def chrs = (1..22).collect()
//     chrs.addAll(['X', 'Y', 'MT'])
//     return chrs
//   }

//   def getChrmListHg38(){
//     def chrs = (1..22).collect()
//     chrs.addAll(['X', 'Y'])
//     return chrs.collect {"chr$it"}
//   }


// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Check if params.item exists and return params.genomes[params.genome].item otherwise
def checkParamReturnFile(item) {
    // Handle deprecation
    if (params.genomeDict && item == "dict") return file(params.genomeDict)
    if (params.genomeFile && item == "fasta") return file(params.genomeFile)
    if (params.genomeIndex && item == "fasta_fai") return file(params.genomeIndex)

    params."${item}" = params.genomes[params.genome]."${item}"
    return file(params."${item}")
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
        'HaplotypeCaller',
        'Manta',
        'Mutect2',
        'Strelka',
        'TIDDIT'
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'bamqc',
        'bcftools',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'sentieon',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'annotate',
        'mapping',
        'recalibrate',
        'variantcalling'
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'ascat',
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

def pr_ch(ch){
    println("pr_ch() printing channel values")
    ch.view{println(it)}
}
// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "subject gender status sample lane fastq1 fastq2"
// or: "subject gender status sample lane bam"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
            if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
        }
        else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
        else "No recognisable extention for input file: ${file1}"

        [idPatient, gender, status, idSample, idRun, file1, file2]
    }
}

// Channeling the TSV file containing Recalibration Tables.
// Format is: "subject gender status sample bam bai recalTables"
def extractRecal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 7)
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

            [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}