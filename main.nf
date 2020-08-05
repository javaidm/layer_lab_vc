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
// model= params.exome ? 'WES' : 'WGS'
model= params.exome ? 'wes' : 'wgs'

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

// Check if bm_dv_against_gatk selected, but deepvariant and haplotypecaller are not in the tools list
if ( ('benchmark_dv_against_hc' in tools) && (!('haplotypecaller' in tools) || !('deepvariant' in tools) ) ) {
    exit 1, "When benchmark_dv_against_hc selected, both HaplotypeCaller and DeepVariant needs to be in the tools list"
}


skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

anno_list = defineAnnoList()
// annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase()} : []
annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim()} : []
if (!checkParameterList(annotate_tools,anno_list)) exit 1, "Unknown tool(s) (${annotate_tools}) to annotate, see --help for more information"

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
        case 'markdups': ch_input_sample = extractMarkDups(tsvFile); break
        case 'recalibrate': ch_input_sample = extractRecal(tsvFile); break
        case 'variantcalling': ch_input_sample = extractBam(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}


(genderMap, statusMap, ch_input_sample) = extractInfos(ch_input_sample)


// ch_input_sample
//     .dump(tag: 'ch_input_sample after extracting info!')
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.fasta_gz = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta_gz ?: null : null
params.fasta_gz_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_gz_fai ?: null : null
params.fasta_gzi = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta_gzi ?: null : null

// Annotation related params (snpEff)
params.snpEff_db = params.genome && ('annotate' in step) ? params.genomes[params.genome].snpEff_db ?: null : null
params.snpEff_cache = params.genome && ('annotate' in step) ? params.genomes[params.genome].snpEff_cache ?: null : null
params.species = params.genome && ('annotate' in step) ? params.genomes[params.genome].species ?: null : null
// VEP
params.vep_cache = params.genome && ('annotate' in step) ? params.genomes[params.genome].vep_cache ?: null : null
params.vep_cache_version = params.genome && ('annotate' in step) ? params.genomes[params.genome].vep_cache_version ?: null : null
// CADD
params.cadd_InDels = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_InDels ?: null : null
params.cadd_InDels_tbi = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_InDels_tbi ?: null : null
params.cadd_WG_SNVs = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_WG_SNVs ?: null : null
params.cadd_WG_SNVs_tbi = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_WG_SNVs_tbi ?: null : null
// CHCO related Pipeline validation using Illumina Hap.py package. GIAB truth sets
params.giab_highconf_vcf = params.genome ? params.genomes[params.genome].giab_highconf_vcf ?: null : null
params.giab_highconf_tbi = params.genome ? params.genomes[params.genome].giab_highconf_tbi ?: null : null
params.giab_highconf_regions = params.genome ? params.genomes[params.genome].giab_highconf_regions ?: null : null
params.chco_highqual_snps = params.genome ? params.genomes[params.genome].chco_highqual_snps ?: null : null


// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa_index = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa_index ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && \
                ('mapping' in step || 'markdups' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) \
                ? params.genomes[params.genome].dbsnp ?: null : null

params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.germline_resource = params.genome ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.known_indels = params.genome && ( 'mapping' in step || 'markdups' in step ) \
                            ? params.genomes[params.genome].known_indels ?: null : null

params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null

// Initialize channels based on params
// ch_acLoci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
// ch_acLociGC = params.ac_lociGC && 'ascat' in tools ? Channel.value(file(params.ac_lociGC)) : "null"
// ch_chrDir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
// ch_chrLength = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
// ch_dbsnp = params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"

ch_acLoci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
ch_acLoci_GC = params.ac_loci_GC && 'ascat' in tools ? Channel.value(file(params.ac_loci_GC)) : "null"
ch_chrDir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
ch_chrLength = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
ch_dbsnp = params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"
// ch_dbsnp_index = params.dbsnp_index && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp_index)) : "null"
ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
// ch_dict = params.dict ? Channel.value(file(params.dict)) : "null"
// ch_fasta_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
// ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
// ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_germline_resource = params.germline_resource &&
                      ('mutect2' in tools || 'mutect2_single' in tools || 'gen_somatic_pon' in tools ) 
                      ? Channel.value(file(params.germline_resource)) : "null"
// if (ch_germline_resource == 'null') error "ch_germline_resource is empty!${params.germline_resource}"
ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"

ch_read_count_pon = params.read_count_pon ? Channel.value(file(params.read_count_pon)) : "null"
ch_somatic_pon = params.somatic_pon ? Channel.value(file(params.somatic_pon)) : "null"
ch_somatic_pon_index = params.somatic_pon_index ? Channel.value(file(params.somatic_pon_index)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"
// padded target, if specified also generate CollectHsMetrics for the recal bams intersected with the padded bed
ch_padded_target_bed = params.padded_target_bed ? Channel.value(file(params.padded_target_bed)) : "null"
ch_bait_bed = params.bait_bed ? Channel.value(file(params.bait_bed)) : "null"
// knownIndels is currently a list of file for smallGRCh37, so transform it in a channel
li_known_indels = []
if (params.known_indels && ('mapping' in step || 'markdups' in step)) params.known_indels.each { li_known_indels.add(file(it)) }

li_known_indels_index = []
if (params.known_indels_index && ('mapping' in step || 'markdups' in step)) params.known_indels_index.each { li_known_indels_index.add(file(it)) }
// ch_known_indels = Channel.empty()
// ch_known_indels_index = Channel.empty()
ch_known_indels = params.known_indels && params.genome == 'smallGRCh37' ? Channel.value(li_known_indels.collect()) 
    : params.known_indels ? Channel.value(file(params.known_indels)) : "null"

ch_known_indels_index = params.known_indels_index && params.genome == 'smallGRCh37' \
    ? Channel.value(li_known_indels_index.collect()) \
    : \
    params.known_indels_index ? Channel.value(file(params.known_indels_index)) : "null"

ch_snpEff_cache = params.snpEff_cache ? Channel.value(file(params.snpEff_cache)) : "null"
// ch_snpeff_cache = params.snpeff_cache ? Channel.fromPath(params.snpeff_cache) : "null"
ch_snpEff_db = params.snpEff_db ? Channel.value(params.snpEff_db) : "null"


ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"
ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"
// Optional files, not defined within the params.genomes[params.genome] scope
ch_cadd_InDels = params.cadd_InDels ? Channel.value(file(params.cadd_InDels)) : "null"
ch_cadd_InDels_tbi = params.cadd_InDels_tbi ? Channel.value(file(params.cadd_InDels_tbi)) : "null"
ch_cadd_WG_SNVs = params.cadd_WG_SNVs ? Channel.value(file(params.cadd_WG_SNVs)) : "null"
ch_cadd_WG_SNVs_tbi = params.cadd_WG_SNVs_tbi ? Channel.value(file(params.cadd_WG_SNVs_tbi)) : "null"

// Optional CHCO files for calculating TP, FP, TN etc against the GIAB
ch_giab_highconf_vcf = params.giab_highconf_vcf ? Channel.value(file(params.giab_highconf_vcf)) : "null"
ch_giab_highconf_tbi = params.giab_highconf_tbi ? Channel.value(file(params.giab_highconf_tbi)) : "null"
ch_giab_highconf_regions = params.giab_highconf_regions ? Channel.value(file(params.giab_highconf_regions)) : "null"
ch_chco_highqual_snps = params.chco_highqual_snps ? Channel.value(file(params.chco_highqual_snps)) : "null"

/* Create channels for various indices. These channels are either filled by the user parameters or 
form inside the build_indices workflow */
ch_fasta_fai = ch_fasta_gz = ch_fasta_gzi = ch_fasta_gz_fai \
= ch_bwa_index = ch_dict = ch_dbsnp_index = ch_germline_resource_index \
=  ch_somatic_pon_index =  Channel.empty()
printSummary()

workflow wf_get_software_versions{
    main:
        GetSoftwareVersions()
    emit:
        GetSoftwareVersions.out
} // end of wf_get_software_versions

workflow wf_build_indexes{
    BuildFastaFai(ch_fasta)
    BuildFastaGz(ch_fasta)
    BuildFastaGzFai(ch_fasta,
                    BuildFastaGz.out)
    BuildFastaGzi(BuildFastaGz.out)

    BuildBWAindexes(ch_fasta)
    BuildDict(ch_fasta)
    BuildDbsnpIndex(ch_dbsnp)
    BuildGermlineResourceIndex(ch_germline_resource)
    BuildKnownIndelsIndex(ch_known_indels)
    BuildSomaticPonIndex(ch_somatic_pon)
    fasta_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : BuildFastaFai.out
    fasta_gz = params.fasta_gz ? Channel.value(file(params.fasta_gz)) : BuildFastaGz.out
    fasta_gz_fai = params.fasta_gz_fai ? Channel.value(file(params.fasta_gz_fai)) : BuildFastaGzFai.out
    fasta_gzi = params.fasta_gzi ? Channel.value(file(params.fasta_gzi)) : BuildFastaGzi.out
    bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BuildBWAindexes.out
    dict = params. dict ? Channel.value(file(params.dict)) :  BuildDict.out

    dbsnp_index = params.dbsnp ? \
                    params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) \
                    : BuildDbsnpIndex.out : "null"        
    
    germline_resource_index = params.germline_resource ? \
                        params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) \
                        : BuildGermlineResourceIndex.out : "null"
    
    
    known_indels_index =   params.known_indels ? \
                            params.known_indels_index ? ch_known_indels_index : \
                            BuildKnownIndelsIndex.out.collect() \
                            : "null"    
    somatic_pon_index = params.somatic_pon_index ? \
                            Channel.value(file(params.somatic_pon_index)) \
                            : BuildSomaticPonIndex.out
    emit:
        fasta_fai = fasta_fai
        fasta_gz = fasta_gz
        fasta_gz_fai = fasta_gz_fai
        fasta_gzi = fasta_gzi
        bwa_index = bwa_index
        dict = dict
        dbsnp_index = dbsnp_index
        germline_resource_index = germline_resource_index
        known_indels_index = known_indels_index
        somatic_pon_index = somatic_pon_index
        
} // end of wf_build_indices

workflow wf_build_intervals{
    take: fasta_fai
    main:
        BuildIntervals(fasta_fai)

        _intervals = params.no_intervals ? "null" : \
                params.intervals && !('annotate' in step) ? \
                Channel.value(file(params.intervals)) : BuildIntervals.out
    
        CreateIntervalBeds(_intervals)
    emit:
        bed_intervals = CreateIntervalBeds.out.flatten()
}

workflow wf_partition_fastq{
    main:
    _pair_reads = Channel.empty()
    // Close the input_pair_read if the starting step is not mapping
    if (step == 'mapping'){           
        if (params.split_fastq){
            PartitionFastQ(ch_input_sample)
            _pair_reads = _pair_reads.mix(PartitionFastQ.out)
            .flatMap{            
                idPatient, idSample, idRun, reads_1, reads_2 ->
                myList= []
                reads_1.each { read_1 -> 
                                split_index = read_1.fileName.toString().minus("r1_split_").minus(".fastq.gz")
                                parent = read_1.parent
                                read_2_fn = read_1.fileName.toString().replace("r1_split_", "r2_split_")
                                read_2 = "${parent}/${read_2_fn}"
                                new_id_run = "${idRun}_${split_index}"
                                myList.add([idPatient, idSample, new_id_run, read_1, file(read_2)])
                            }
                myList     
            }
        }else{
            _pair_reads = ch_input_sample
        }
    }
    emit:
        pair_reads = _pair_reads
} // end of wf_partition_fastq

workflow wf_merge_mapped_reads{
    
    take:
        take: bam_mapped
    main:
        // STEP 1.5: MERGING BAM FROM MULTIPLE LANES
        (single_bams, multiple_bams) = 
        bam_mapped.groupTuple(by:[0, 1])
        .branch{
            _: it[2].size() == 1
            __: it[2].size() > 1
        }

        // // ch_multiple_bams.subscribe{ println it}
        single_bams = single_bams.map{
            idPatient, idSample, idRun, bam ->
            [idPatient, idSample, bam]
        }
        
        // tell the mergeBamMapped to generate the .bams
        multiple_bams = 
        multiple_bams.map{idPatient, idSample, idRun, bams -> 
            [ idPatient, idSample, idRun, '.bam', bams ]
        }

        MergeBamMapped(multiple_bams)
        
        _merged_bams = Channel.empty()
        _merged_bams = MergeBamMapped.out.mix(single_bams)
        IndexBamFile(_merged_bams)
        

        // Creating a TSV file to restart from this step
        _merged_bams
        .map { idPatient, idSample, bamFile ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}.bai"
            bam_file = file(bam)
            bai_file = file(bai)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam_file}\t${bai_file}\n"
        }.collectFile(
            name: 'mapped_bam.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )
        // if params.filter_bams is present, generate the relevant bam as well
        if (params.filter_bams){
            _merged_bams
            .map { idPatient, idSample, bamFile ->
                status = statusMap[idPatient, idSample]
                gender = genderMap[idPatient]
                bam = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}_pq${params.bam_mapping_q}.bam"
                bai = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}_pq${params.bam_mapping_q}.bai"
                "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
            }.collectFile(
                name: "mapped_bam_${params.bam_mapping_q}.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
            )   
        }
        // exit 1, 'leaving early!'
    emit:
        merged_bams = IndexBamFile.out
        merged_bams_qc = _merged_bams
} // end of wf_map_reads

workflow wf_qc_filter_mapped_reads{
    
    take:
       take: _bam_mapped
    main:
        // bam_mapped
        // .dump('bam_mapped_input: ')
        // Run the bam QC filtering
        FilterBamRead1(_bam_mapped)
        FilterBamRead2(_bam_mapped)
        filtered_reads = FilterBamRead1.out.filtered_reads
                        .mix(FilterBamRead2.out.filtered_reads)
                        .groupTuple(by: [0,1,2])
                        .dump(tag: 'filtered_reads_merge: ')
        // Merge Reads1 and Reads2
        MergeFilteredBamReads(filtered_reads)
        // Here we filter out seconday and supplemntary reads
        FilterOutSecondaryAndSupplementaryReads(
            MergeFilteredBamReads.out.filtered_bam)
        filtered_out_bams = FilterOutSecondaryAndSupplementaryReads
                            .out.filtered_bam
                            .groupTuple(by:[0, 1])
                            // .dump(tag: 'qc_bams_for_merging')
        
        // tell the mergeBamMapped to generate the .bams
        filtered_out_bams = 
        filtered_out_bams.map{idPatient, idSample, idRun, bams -> 
            [ idPatient, idSample, idRun, "_pq${params.bam_mapping_q}.bam", bams ]
        }
        .dump(tag: 'qc_bams_for_merging')
        // STEP 1.5: MERGING BAM FROM MULTIPLE LANES
        MergeBamMapped(filtered_out_bams)
        IndexBamFile(MergeBamMapped.out)
    emit:
        merged_bams = IndexBamFile.out
} // end of wf_map_reads

workflow wf_mark_duplicates{
    take: _bams
    main:
        MarkDuplicates(_bams)
    emit:
        dm_bams = MarkDuplicates.out.marked_bams
}


workflow wf_deepvariant{
    take: _dm_bams
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _fasta_gz
    take: _fasta_gz_fai
    take: _fasta_gzi

    main:
        DV_MakeExamples(_dm_bams,
                _fasta,
                _fasta_fai,
                _fasta_gz,
                _fasta_gz_fai,
                _fasta_gzi,
                _target_bed)
        DV_CallVariants(DV_MakeExamples.out.shared_examples)
        DV_PostprocessVariants(
            DV_CallVariants.out.variant_tf_records,
            _fasta_gz,
            _fasta_gz_fai,
            _fasta_gzi
        )

    emit:
        vcf = DV_PostprocessVariants.out.vcf
        vcf_to_annotate = DV_PostprocessVariants.out.vcf_to_annotate
} // end of wf_deepvariant

workflow wf_mpileup{
    take: _int_bam_recal
     take: _fasta
    take: _fasta_fai
    main:
        Mpileup(_int_bam_recal,
                 _fasta,
                _fasta_fai 
        )
        
        MergeMpileup(
            Mpileup.out.groupTuple(by:[0, 1]) 
        )
} // end of wf_mpileup

workflow wf_recal_bams{

    take: _md_bams
    take: _bed_intervals
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    take: _known_indels
    take: _known_indels_index
          
    main:
        // Create a scattered pattern 
        _int_md_bams = _md_bams.combine(_bed_intervals)
        BaseRecalibrator(
            _int_md_bams,
            _fasta,
            _fasta_fai,
            _dict,
            _dbsnp,
            _dbsnp_index,
            _known_indels,
            _known_indels_index
        )

        table_gather_bqsr_reports = 
            !params.no_intervals ? BaseRecalibrator.out.groupTuple(by:[0, 1]) : BaseRecalibrator.out
        GatherBQSRReports(table_gather_bqsr_reports)
        
        bam_apply_bqsr = _md_bams
                        .join(GatherBQSRReports.out.recal_table, by:[0,1])
        // bam_apply_bqsr = Sambamba_MD.out
        //                 .join(GatherBQSRReports.out.recal_table, by:[0,1])

        bam_apply_bqsr = bam_apply_bqsr.combine(_bed_intervals)
        ApplyBQSR(
            bam_apply_bqsr,
             _fasta,
            _fasta_fai,
            _dict
        )

        bam_merge_bam_recal = ApplyBQSR.out.groupTuple(by:[0, 1])
        /* When using intervals, merge (and in the same process index bam files)
            Which one of the MergeBamRecal, or IndexBamRecal runs is controlled by the 'when' clause
            in the processes. We finally mix them together as we know they are mutulaly exclusive.
        */
        MergeBamRecal(bam_merge_bam_recal)
        // When not using intervals, just index the bam coming from ApplyBQSR
        IndexBamRecal(bam_merge_bam_recal)
    
        bam_recal = MergeBamRecal.out.bam_recal.mix(IndexBamRecal.out.bam_recal)
        bam_recal_qc = MergeBamRecal.out.bam_recal_qc.mix(IndexBamRecal.out.bam_recal_qc)
        
        BamRecalOnTarget(bam_recal,
                         ch_padded_target_bed)
        
        // Creating a TSV file to restart from this step
        bam_recal.map { idPatient, idSample, bamFile, baiFile ->
            gender = genderMap[idPatient]
            status = statusMap[idPatient, idSample]
            bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
        }.collectFile(
            name: 'recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

        bam_recal
            .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
                idPatient, idSample, bamFile, baiFile ->
                status = statusMap[idPatient, idSample]
                gender = genderMap[idPatient]
                bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
                bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
                ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
        }
    emit:
        // output structure
        // tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai"), emit: bam_recal
        // tuple idPatient, idSample, file("${idSample}.recal.bam"), emit: bam_recal_qc
        bam_recal = bam_recal
        bam_recal_qc = bam_recal_qc
        bam_recal_on_target_qc = BamRecalOnTarget.out.bam_recal_on_target_qc
} // end of wf_recal_bams

workflow wf_qc_recal_bams{

    take: _bam_raw_qc // tuple idPatient, idSample, file(bam), file(bam // take the raw bams
    take: _bam_recal_qc // tuple idPatient, idSample, file(bam)
    take: _bam_recal_on_target_qc // tuple idPatient, idSample, file(bam)
    take: _target_bed
    take: _bait_bed
    take: _fasta
    take: _fasta_fai
    take: _dict

    main:
        // _bam_raw_qc.dump(tag: 'raw_bam_qc: ')
        SamtoolsStats(_bam_recal_qc)
        CollectAlignmentSummaryMetrics(
            _bam_recal_qc,
            _fasta,
            _fasta_fai,
            _dict
        )
        CollectInsertSizeMetrics(
             _bam_recal_qc
        )
        /* Now for Rebecca's on-target assessment quality metric, 
            we need to run CollectHSMetrics on three type of bams
           
           a) Raw, unmarked, un-recalibrared bams
           b) quality filtered, marked duplicated, reclaibrated bams
           c) bams in above step 'b' intersected with the targets_bed to get the on-target bams

           Skip it if user has specified on_target_assessment in skip_qc param
        */
        _bam_collect_hs_metrics = _bam_recal_qc
        if (!('on_target_assessment' in skipQC)){
            _bam_collect_hs_metrics = _bam_collect_hs_metrics
                                       .mix(_bam_raw_qc) 
            _bam_collect_hs_metrics = _bam_collect_hs_metrics
                                    .mix(_bam_recal_on_target_qc)
        }
                                    
        

        CollectHsMetrics(
            _bam_collect_hs_metrics,
            _fasta,
            _fasta_fai,
            _dict,
            _target_bed,
            _bait_bed,
            // '' // this is the suffix to the output file
        )
        BamQC(
            _bam_recal_qc,
            _target_bed
        )

    emit:
        samtools_stats =  SamtoolsStats.out
        alignment_summary_metrics = CollectAlignmentSummaryMetrics.out
        insert_size_metrics = CollectInsertSizeMetrics.out[0]
        insert_size_metrics_pdf = CollectInsertSizeMetrics.out[1]
        hs_metrics = CollectHsMetrics.out
        bam_qc = BamQC.out
} // end of wf_qc_recal_bams

workflow wf_haplotypecaller{
    take: _int_bam_recal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    // take: _bed_intervls
    main:
        HaplotypeCaller(
            _int_bam_recal,
            _fasta,
            _fasta_fai,
            _dict,
            _dbsnp,
            _dbsnp_index
        )

    emit: 
        gvcf_HC = HaplotypeCaller.out.gvcf_HC
        gvcf_GenotypeGVCFs = HaplotypeCaller.out.gvcf_GenotypeGVCFs
} //  end of wf_haplotypecaller

workflow wf_individually_genotype_gvcf{
    // input structure
    // tuple val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf")
    take: _gvcf_HC
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    take: _target_bed
    main:
        IndividuallyGentoypeGVCF(
            _gvcf_HC,
            _fasta,
            _fasta_fai,
            _dict,
            _dbsnp,
            _dbsnp_index
        )
        
        _vcf_concatVCF = 
            IndividuallyGentoypeGVCF.out.vcf_HaplotypeCaller
            .groupTuple(by: [0,1])
            .map{ idPatient, idSample, vcfs -> 
                ['HaplotypeCaller_Individually_Genotyped', idPatient, idSample, vcfs]
            }
            .dump(tag: "vcf_ConcatVcf")

        ConcatVCF(
                _vcf_concatVCF,
                _fasta_fai,
                _target_bed,
                'HC', // prefix for output files
                'vcf', // extension for the output files
                'HC_individually_genotyped_vcf' // output directory name
                )
        // GvcfToVcf(
        //         ConcatVCF.out.concatenated_vcf_without_index
        //         )
    emit:
        sample_vcf_HC = IndividuallyGentoypeGVCF.out.vcf_HaplotypeCaller
        // sample_vcf_HC = GvcfToVcf.out.vcf_HaplotypeCaller
}



workflow wf_genotype_gvcf{
    // input _scattered_gvcf_HC structure 
    // tuple idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf")
    take: _scattered_gvcf_HC
    take: _target_bed
     take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    // take: _bed_intervls
    main:
    // Transform the haplotypecaller output in the above format to a format where gvcf's are gathered per interval, so carry
    // multiple gvcf (for multiple samples)
    // When grouping gvcf's for per interval basis, we only keep one interval bed as the list carry same intervals for each group
    // This interval then become the key, and a disnguishing element for this joint genotyping
    // patientSampleIdMap will have all the samples names that has a gvcf in a particular group

    mapped_gvcf_GenotypeGVCFs = _scattered_gvcf_HC
                                .map{idPatient, idSample, intervalBed, gvcf ->
                                patientSampleIdMap = [:]
                                patientSampleIdMap['idPatient'] = idPatient
                                patientSampleIdMap['idSample'] = idSample
                                [intervalBed.baseName, intervalBed, patientSampleIdMap , gvcf]}
                                .groupTuple(by:[0])
                                .map{interval_name, liIntervalBed, patientSampleIdMap, liGvcf -> 
                                    [interval_name, liIntervalBed.first(), patientSampleIdMap, liGvcf]
                                    }
                                // .dump(tag: 'Collected HaplotypeCaller output')
    
    
    GenomicsDBImport(mapped_gvcf_GenotypeGVCFs)                              
    GenotypeGVCFs(
        GenomicsDBImport.out,
        _fasta,
        _fasta_fai,
        _dict,
        _dbsnp,
        _dbsnp_index
    )
    // GenotypeGVCFs output
    //tuple val("HaplotypeCaller"),  val(patientSampleIdMap), file(interval_bed), file("vcf"), file("vcf.idx")
    
    // A Cohort vcf
    vcf_cohort_concatenate_vcf = GenotypeGVCFs.out.vcf_GenotypeGVCFs
    .map{caller, li_patient_sample_id_map, interval_bed,  vcf, vcf_idx ->
        [vcf]
    }
    .collect()
    // .dump(tag: 'cohor_vcf: ')

    CohortConcatVCF(
        vcf_cohort_concatenate_vcf,
        _fasta_fai,
        _target_bed,
    )

    // Per Sample vcf (but jointly genotyped)
     ch_select_variants = GenotypeGVCFs.out.vcf_GenotypeGVCFs
    .flatMap{ caller, li_patient_sample_id_map, interval_bed,  vcf, vcf_idx ->
        per_sample_list=[]
        li_patient_sample_id_map.each { entry ->
            per_sample_list.add([caller, entry.idPatient, entry.idSample, interval_bed, vcf, vcf_idx])
        }
        per_sample_list
    }
    // .flatten()
    // .dump(tag: 'ch_select_variants')

    SelectVariants(
        ch_select_variants,
        _fasta,
        _fasta_fai,
        _dict
    )
    // SelectVariants output
    //tuple val("HaplotypeCaller_Jointly_Genotyped"), id_patient, id_sample, file("vcf")
    vcf_ConcatenateVCFs = SelectVariants.out.vcf_SelectVariants.groupTuple(by:[0, 1, 2])
    // Now add the individually called vcfs too
    // vcf_ConcatenateVCFs = vcf_ConcatenateVCFs.mix(vcf_HaplotypeCaller)
    // if (!params.no_gvcf){ // if user specified, noGVCF, skip saving the GVCFs from HaplotypeCaller
    //     vcf_ConcatenateVCFs = vcf_ConcatenateVCFs.mix(gvcf_HaplotypeCaller)
    // }
    // vcf_ConcatenateVCFs.dump('concat_vcf: ')

    ConcatVCF(
        vcf_ConcatenateVCFs,
        _fasta_fai,
        _target_bed,
        'HC_jointly_genotyped', // prefix for output files
        'vcf', // extension for the output files
        'HC_jointly_genotyped_gvcf' // output directory name
        )
    
    // Create a channel to hold GIAB samples for validation using hap.py
    // hc_jointly_genotyped_vcfs = ConcatVCF.out.concatenated_vcf_with_index
    //                       .filter{  "${it[0]}" == 'HaplotypeCaller_Jointly_Genotyped'}
    //                     //   .dump(tag: 'hc_jointly_genotyped_vcfs')

    emit:
    vcfs_with_indexes = ConcatVCF.out.concatenated_vcf_with_index
    vcfs_without_indexes = ConcatVCF.out.concatenated_vcf_without_index
    // cohort_vcf_with_index = CohortConcatVCF.out.cohort_vcf_with_index
    // cohort_vcf_without_index = CohortConcatVCF.out[1]
} // end of wf_haplotypecaller

workflow wf_mutect2_single{
    take: _int_bam_recal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _germline_resource
    take: _germline_resource_index
    take: _target_bed
    main:
        Mutect2Single(
            _int_bam_recal,
            _fasta,
            _fasta_fai,
            _dict,
            _germline_resource,
            _germline_resource_index
        )
        // group scattered vcf's (per interval) by the 
        // idPatient, and idSample

        _concat_vcf = Mutect2Single.out
                      .groupTuple(by: [0, 1])
                      
        _concat_vcf = Channel.from('Mutect2_single')
                      .combine(_concat_vcf)
                    //   .dump(tag: 'vcf_concat_vcf')
        ConcatVCF(
            _concat_vcf,
            _fasta_fai,
            _target_bed,
            'Mu2_single', // prefix for output files
            'vcf', // extension for the output files
            'Mutect2_single_mode' // output directory name
            )
    emit:
        vcf = ConcatVCF.out.concatenated_vcf_with_index
        vcf_without_index = ConcatVCF.out.concatenated_vcf_without_index
} // end of wf_mutect2_single

workflow wf_somatic_pon{
    take:  _vcfs_normal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _germline_resource
    take: _germline_resource_index
    take: _target_bed

    main:
        SomaticPonGenomicsDBImport(
            _vcfs_normal,
            _target_bed
        )

        CreateSomaticPON(
            SomaticPonGenomicsDBImport.out,
            _fasta,
            _fasta_fai,
            _dict,
            _germline_resource,
            _germline_resource_index
            )
    emit:
        
        somatic_pon_gdb = SomaticPonGenomicsDBImport.out
        somatic_pon = CreateSomaticPON.out
        
} // end of wf_somatic_pon

// workflow wf_mutect2_tumor_normal{
//     take: _int_pair_bam
//     take: _fasta
//     take: _fasta_fai
//     take: _dict
//     take: _germline_resource
//     take: _germline_resource_index
//     take: _pon_somatic
//     take: _pon_somatic_index
//     take: _target_bed
//     main:
//         Mutect2(
//             _int_pair_bam,
//             _fasta,
//             _fasta_fai,
//             _dict,
//             _germline_resource,
//             _germline_resource_index,
//             _pon_somatic,
//             _pon_somatic_index
//         )

//         // Gather mutect2 outputs
//         _mutect2_output = Mutect2.out.[0].groupTuple(by: [0,1,2])
//         ConcatVCF(
//             _concat_vcf,
//             _fasta_fai,
//             _target_bed,
//             'Mu2', // prefix for output files
//             'vcf', // extension for the output files
//             'Mutect2_tumor_normal' // output directory name
//             )

//         // Gather mutect2 stats
//         _mutect2_stats = Mutect2.out.[1].groupTuple(by: [0,1,2])
//         MergeMutect2Stats(_mutect2_stats)
//         // Gather pileup summaries for each tumor file to calculate contamination
//         // Summarizes counts of reads that support reference, 
//         // alternate and other alleles for given sites
//         PileupSummariesForMutect2(_int_pair_bam)

//         // Gather pileup summaries
//         _pileupsummaries = PileupSummariesForMutect2.out.groupTuple(by: [0,1])
//         MergePileupSummaries(_pileupsummaries,
//                              _dict)

//         CalculateContamination(MergePileupSummaries.out)

//         /* Now we need to gather the unfiltered vcfs (for TN pairs), and get corresponding
//             merged stats, as well as contamination tables
//             ConcatVCF: 
//             tuple variantCaller, idPatient, idSample, file("${outFile}.gz"), file("${outFile}.gz.tbi"), emit: concatenated_vcf_with_index
//         */
//         _filter_mutect2_calls =         
//         FilterMutect2Calls

//     emit:
//         filtered_vcf = FilterMutect2Calls.out
// } // end of wf_mutect2_tumor_normal

workflow wf_hap_py{
    take: deepvariant_vcfs
    take: haplotypecaller_vcfs
    take: _target_bed
    take: _bait_bed
    take: _fasta
    take: _fasta_fai
    main:
    /* Take all vcf's from deepvariant and haplotypecaller, but filter out all but
        those that has 'GIAB' in the sample name.
        The truth set is going to be the GIAB vcfs obtained from their website 
    
    vcfs_giab structure: ['DeepVariant', 'GIABO1', 'GIABO1_S53', vcf.gz, vcf.gz.tbi]    
        */
    vcfs_giab = Channel.empty()
                    .mix(
                            deepvariant_vcfs,
                            haplotypecaller_vcfs
                        )
                    .filter{
                        // The sampleID is the 3rd component of the channel elements
                        "${it[2]}".contains('GIAB')
                    }
                    .dump(tag: "vcfs_giab")

    /*
    The vcfs_giab structure:
    ['Variant Caller', idPatiend, idSample, vcf.gz, vcf.gz.tbi]
    */

    /*
    hap_py_against_giab structure
    [benchmarked_against, 'source caller', idPatient, idSample, src.vcf.gz, src.vcf.gz.tbi, 
    ['GIAB', 'DeepVariant', 'GIABO1', 'GIABO1_S53', GIABO1_S53.vcf.gz, \
    highconf_PGandRTGphasetransfer.vcf.gz, \
    highconf_PGandRTGphasetransfer.vcf.gz.tbi, \
    highconf_nosomaticdel.bed]
    */
    hap_py_against_giab = Channel.from('GIAB') // this will be what the vcfs are being validated against
                            .combine(vcfs_giab) // above filtered vcfs
                            .combine(ch_giab_highconf_vcf) // truth set vcf
                            .combine(ch_giab_highconf_tbi) // truth set index
                            .combine(ch_giab_highconf_regions) // high quality regions from GIAB
                            .dump(tag: 'hap_py_against_giab')

    vcfs_dv_against_hc = deepvariant_vcfs
                      .join(haplotypecaller_vcfs, by:[1,2])                      
                      .map{
                            idPatient, idSample, val_dv, dv_vcf, dv_tbi, val_hc, hc_vcf, hc_tbi ->
                            [ val_dv, idPatient, idSample, dv_vcf, dv_tbi, hc_vcf, hc_tbi]
                            }
                    //   .dump(tag: 'vcfs_against_hc')
    hap_py_dv_against_hc =  Channel.from('hc') // this will be what the vcfs are being validated against
                            .combine(vcfs_dv_against_hc) // above filtered vcfs
                            .combine(ch_giab_highconf_regions) // high quality regions from GIAB
                            .dump(tag: 'hap_py_dv_against_hc')

    // For this channel to go ahead, both the haplotypecaller and the deepvariant should be in tools
    if (!(('benchmark_dv_against_hc' in tools) && ('benchmark_dv_and_hc_against_giab' in tools))) {
        hap_py_dv_against_hc = Channel.empty()
    }

    hap_py_combined = hap_py_against_giab.mix(hap_py_dv_against_hc)
                       .dump(tag:'hap_py_combined')

    HapPy(hap_py_combined,
          _target_bed,
          _bait_bed,
         _fasta,
         _fasta_fai
        )
    emit:
        HapPy.out
} // end of wf_hap_py

workflow wf_vcf_stats{
    // deepvariant output
    //tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz")
    take: _deepvariant_vcfs
    take: _haplotypecaller_vcfs
    main:
        combined_vcfs = _deepvariant_vcfs
                        .mix(_haplotypecaller_vcfs)
                        // .dump(tag: 'vcfs_stats:')
        BcftoolsStats(combined_vcfs)
        Vcftools(combined_vcfs)
    emit:
        bcfootls_stats = BcftoolsStats.out
        vcfootls_stats = Vcftools.out
}

workflow wf_germline_cnv{
    // deepvariant output
    //tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz")
    take: _bam_recal
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    
    main:
        StrelkaSingle(
            _bam_recal,
            _fasta,
            _fasta_fai,
            _target_bed
    )

        MantaSingle(
            _bam_recal,
            _fasta,
            _fasta_fai,
            _target_bed
        )

        TIDDIT(
            _bam_recal,
             _fasta,
            _fasta_fai,
        )
} // end of wf_germline_cnv

workflow wf_gatk_somatic_cnv{
    take: _md_bam
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _read_count_somatic_pon
    
    main:
        /* GATK Somatic Copy Number related calls */
        /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
        /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/

        PreprocessIntervals(_target_bed,
                            _fasta,
                            _fasta_fai,
                            _dict
                            )

        CollectReadCounts(_md_bam,
                        PreprocessIntervals.out)

        // Create pon if specified
        (normal_read_counts, tumor_read_counts) = 
        CollectReadCounts.out
        .branch{
             _: statusMap[it[0], it[1]] == 0
            __: statusMap[it[0], it[1]] == 1
        }
        
        // normal_read_counts = 
        // CollectReadCounts.out.collect()
        // .filter{
        //      statusMap[it[0], it[1]] == 0
        // }

        normal_read_counts
        .dump(tag: 'normal_read_counts: ')
        
        normal_read_counts_hdf5 = 
                normal_read_counts
                .map{idPatient, idSample, hdf5 -> 
                    [hdf5]
                }
                .dump(tag: 'hdf5: ')
        
        tumor_read_counts
        .dump(tag: 'tumor_read_counts: ')

        CreateReadCountPon(
            normal_read_counts_hdf5.collect()
        )
        
        DenoiseReadCounts(
            CollectReadCounts.out.sample_read_counts,
            _read_count_somatic_pon
        )

        PlotDenoisedCopyRatios(DenoiseReadCounts.out.denoised_cr,
                                _dict)

        denoised_cr_model_segments = DenoiseReadCounts.out.denoised_cr
                                    .map{ idPatien, idSample, std_cr, denoised_cr ->
                                        [idPatien, idSample, denoised_cr]
                                        }

        ModelSegments(denoised_cr_model_segments)
        plot_modeled_segments = ModelSegments.out.modeled_seg
                                .join(DenoiseReadCounts.out.denoised_cr, by: [0,1])
                                .map{idPatient, idSample, cr_seg, model_final_seg, std_cr, denoised_cr -> 
                                    [idPatient, idSample, model_final_seg, denoised_cr]
                                }
                                .dump(tag: 'plot_modeled_segments')

        PlotModeledSegments(plot_modeled_segments,
                            _dict)
        call_cr_seg = ModelSegments.out.modeled_seg                        
                                .map{idPatient, idSample, cr_seg, model_final_seg -> 
                                    [idPatient, idSample, cr_seg]
                                }
                                .dump(tag: 'call_cr_seg')

        CallCopyRatioSegments(call_cr_seg)
} // end of wf_gatk_somatic_cnv

workflow wf_savvy_somatic_cnv{
    take: _md_bam
    main:
     /* SavvyCNV Somatic Copy Number related calls */
    /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
    /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/
    SavvyCNVCoverageSummary(_md_bam)
    SavvyCNV(SavvyCNVCoverageSummary.out.collect())
} // end of wf_germline_cnv

workflow wf_annotate{
    // to do
    take: _md_bam
    take: _fasta
    take: _fasta_fai
    main:
        /* Annotations */
        ch_vcfs_to_annotate = Channel.empty()
        if (step == 'annotate') {
            // ch_vcfs_to_annotate = getVCFsToAnnotate(params.outdir, annotate_tools)
            ch_vcfs_to_annotate = Channel.empty()
            // vcf_no_annotate = Channel.empty()

            if (tsvPath == [] || !tsvPath) {
            // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
            // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
            // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
            // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
            // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
            // This field is used to output final annotated VCFs in the correct directory
                // log.info ("annotate_tools: ${annotate_tools}")
                ch_vcfs_to_annotate = 
                Channel.empty().mix(
                Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
                    .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                // Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
                Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*diploidSV.vcf.gz")
                    .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
                    .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
                    .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
                Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
                    .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
                )
                .filter {
                    annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools)
                }
            } else if (annotate_tools == []) {
            // Annotate user-submitted VCFs
            // If user-submitted, Sarek assume that the idSample should be assumed automatically
            ch_vcfs_to_annotate = Channel.fromPath(tsvPath)
                .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
            } else exit 1, "specify only tools or files to annotate, not both"
        }
        
        // log.info "annotate_tools: ${annotate_tools}"
        // ch_vcfs_to_annotate.dump(tag: 'ch_vcf_to_annotate')
        ch_vcfs_to_annotate = ch_vcfs_to_annotate.mix(
                            ConcatVCF.out.vcf_concatenated_to_annotate,
                            StrelkaSingle.out.map{
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                    [variantcaller, idSample, vcf[1]]
                            },
                            MantaSingle.out.map {
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                [variantcaller, idSample, vcf[2]]
                            },
                            TIDDIT.out.vcfTIDDIT.map {
                                variantcaller, idPatient, idSample, vcf, tbi ->
                                [variantcaller, idSample, vcf]
                                }
                            )

        // ch_vcf_snpEff = ch_vcfs_to_annotate.mix(ConcatVCF.out.vcf_concatenated_to_annotate)
        ch_vcf_snpEff = ch_vcfs_to_annotate
        // ch_vcf_snpEff = Channel.empty()

    ch_vcf_vep = ch_vcf_snpEff.map {
                variantCaller, idSample, vcf ->
                [variantCaller, idSample, vcf, null]
        }
    //     // PrintCh(ch_vcf_snpeff)
        SnpEff( ch_vcf_snpEff,
                ch_snpEff_cache,
                ch_snpEff_db
                )
        CompressVCFsnpEff(SnpEff.out.snpEff_vcf)
        
        VEP(ch_vcf_vep,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_InDels,
            ch_cadd_InDels_tbi,
            ch_cadd_WG_SNVs,
            ch_cadd_WG_SNVs_tbi,
            _fasta,
            _fasta_fai
        )

        VEPmerge(CompressVCFsnpEff.out.compressVCFsnpEffOut,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_InDels,
            ch_cadd_InDels_tbi,
            ch_cadd_WG_SNVs,
            ch_cadd_WG_SNVs_tbi,
            _fasta,
            _fasta_fai
        )
        CompressVCFvep(VEP.out.vep_vcf.mix(
                        VEPmerge.out.vep_vcf_merge)
                        )
} // end of wf_germline_cnv

workflow wf_multiqc{
    // to do
    take: software_versions
    take: bam_qc
    take: fastqc_out 
    take: bcftools_out 
    take: vcftools_out 
    // take: dm_bam_stats
    take: samtools_stats
    take: alignment_summary_metrics
    take: insert_size_metrics
    take: hs_metrics

    main:
        MultiQC(
            Channel.value(""),
            software_versions,
            bam_qc,
            fastqc_out,
            bcftools_out,
            vcftools_out,
            // dm_bam_stats,
            samtools_stats,
            alignment_summary_metrics,
            insert_size_metrics,
            hs_metrics
            // SnpEff.out.snpEff_report,
        )
} // end of wf_germline_cnv

workflow{

    wf_get_software_versions()
    wf_build_indexes()
    ch_fasta_fai =  wf_build_indexes.out.fasta_fai
    // The following three mainly used by the DeepVariant
    ch_fasta_gz =  wf_build_indexes.out.fasta_gz
    ch_fasta_gz_fai =  wf_build_indexes.out.fasta_gz_fai
    ch_fasta_gzi =  wf_build_indexes.out.fasta_gzi
    ch_bwa_index =  wf_build_indexes.out.bwa_index
    ch_dict =  wf_build_indexes.out.dict
    ch_dbsnp_index = wf_build_indexes.out.dbsnp_index
    ch_germline_resource_index = wf_build_indexes.out.germline_resource_index
    ch_known_indels_index = wf_build_indexes.out.known_indels_index.collect()
    ch_somatic_pon_index = wf_build_indexes.out.somatic_pon_index
    
    wf_build_intervals(ch_fasta_fai)
    ch_bed_intervals = wf_build_intervals.out.bed_intervals
    if (params.no_intervals && step != 'annotate') 
        ch_bed_intervals = Channel.from(file("no_intervals.bed"))
    wf_partition_fastq()
    ch_input_pair_reads = wf_partition_fastq.out.pair_reads

    FastQCFQ(ch_input_sample)     
    MapReads(
        ch_input_pair_reads, 
            ch_fasta,
            ch_fasta_fai,
            ch_bwa_index
    )
    wf_merge_mapped_reads(MapReads.out.bam_mapped)
    wf_qc_filter_mapped_reads(MapReads.out.bam_mapped)

    /* 
    Now depending upong the params.filter_bams, either send
    unfiltered or filtered bams for the downstream analysis
    */
    
    ch_merged_bams = Channel.empty()
    if (params.filter_bams)
        ch_merged_bams = wf_qc_filter_mapped_reads.out.merged_bams
    else
        ch_merged_bams = wf_merge_mapped_reads.out.merged_bams
    
    // If the pipeline is being started from the step 'markdups'
    // use bams from the provided tsv
    if (step == 'markdups'){
       ch_merged_bams = ch_input_sample
    }

   // MarkDuplicate Bams using gatk MarkDuplicates
    wf_mark_duplicates(ch_merged_bams)
    // wf_mark_duplicates output
    //tuple idPatient, idSample, file("md.bam"), file("md.bai")
    
    // Filter Bams based upong the bam_mapping_q, and params.filter_bams
    // Rebecca's request
    // wf_filter_bams(wf_mark_duplicates.out.dm_bams)
    
    // Use Duplicated Marked Bams for DeepVariant
    wf_deepvariant(
                wf_mark_duplicates.out.dm_bams,
                ch_target_bed,
                ch_fasta,
                ch_fasta_fai,
                ch_fasta_gz,
                ch_fasta_gz_fai,
                ch_fasta_gzi
                )

    ch_int_dm_bams = wf_mark_duplicates.out.dm_bams.combine(ch_bed_intervals)
    
    wf_mpileup(
        ch_int_dm_bams,
        ch_fasta,
            ch_fasta_fai
    )

    wf_recal_bams(
            wf_mark_duplicates.out.dm_bams, // recalibrated bams
            ch_bed_intervals,
            ch_fasta,
            ch_fasta_fai,         
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_index,
            ch_known_indels,
            ch_known_indels_index
    )
    // ch_int_bam_recal carries a cross product of recalibrated bams and the interval files
    ch_int_bam_recal = wf_recal_bams.out.bam_recal
                    .combine(ch_bed_intervals)
   
    // Run QC metrics generation processes on the recalibrated bams

    wf_qc_recal_bams(
            wf_merge_mapped_reads.out.merged_bams_qc, // raw, unmarked, un-recalibrated bams
            wf_recal_bams.out.bam_recal_qc,
            wf_recal_bams.out.bam_recal_on_target_qc,
            ch_target_bed,
            ch_bait_bed,
            ch_fasta,
            ch_fasta_fai,         
            ch_dict
            )

    // Handle the Mutect2 related workflows

    // separate BAM by status
    // bamNormal = Channel.empty()
    // bamTumor = Channel.empty ()
    (bam_normal, bam_tumor) = 
        wf_recal_bams.out.bam_recal.branch{
            _:  statusMap[it[0], it[1]] == 0
            __: statusMap[it[0], it[1]] == 1
        }
    bam_normal.dump(tag: 'bam_normal: ')
    bam_tumor.dump(tag: 'bam_tumor: ')
    
    // We need to prepare samples going to the mutect2 single sample mode.
    // This depends upon what is included in the params.tools. 
    // Is it just the mutect2_single
    // or gen_somatic_pon or both
    // For generating a somatic pon, we only need normal bams.
    // When running mutect2_single without the intent of generating a somatic_pon,
    // run it on all bams (tumors plus normals)

    ch_int_mutec2_single_bams = Channel.empty()
    if ('mutect2_single' in tools ) // use all bams (normal plus tumor)
        ch_int_mutec2_single_bams = wf_recal_bams.out.bam_recal.combine(ch_bed_intervals)
    else if ('gen_somatic_pon' in tools && !('mutect2_single' in tools)) // use normal only
        ch_int_mutec2_single_bams = bam_normal.combine(ch_bed_intervals)
    // ch_int_mutec2_single_bams.dump(tag: 'ch_int_mutec2_single_bams: ')
    
    wf_mutect2_single(
        ch_int_mutec2_single_bams,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_index,
        ch_target_bed
    )
    // For somatic pon, only filter the mutect2 single vcf's for normal samples
    // wf_mutect_single output tuple is:
    //tuple variantCaller, idPatient, idSample, file(vcf.gz), file(vcf.gz.tbi)
    // so we filter on the second and third element and check against the statusMap
    somatic_pon_vcf = wf_mutect2_single.out.vcf
                    .filter{
                            statusMap[it[1], it[2]] == 0
                        }
                    .collect()
                    .dump(tag: "somatic_pon_vcf")

    wf_somatic_pon(
            somatic_pon_vcf,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_germline_resource,
            ch_germline_resource_index,
            ch_target_bed
    )
    // Mutect2 Tumor/Normal mode
    // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
    // Remapping channel to remove common key idPatient
    pair_bam = bam_normal.cross(bam_tumor).map {
        normal, tumor ->
        [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
    }

    // spread the pair_bam across intervals for parallelization
    ch_int_pair_bam = pair_bam.combine(ch_bed_intervals)
    

    // bam_HaplotypeCaller => [idPatient, idSample, recalibrated bam, bam.bai, single interval to operate on]
    
  
    wf_haplotypecaller(
        ch_int_bam_recal,
        ch_fasta,
        ch_fasta_fai,         
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_index
    )
    // Create individual gvcfs without any genotyping
     gvcf_ConcatVCF = 
            wf_haplotypecaller.out.gvcf_GenotypeGVCFs
            .groupTuple(by: [0,1])
            .map{ idPatient, idSample, interval_beds,  gvcfs -> 
                ['HaplotypeCaller_gvcf', idPatient, idSample, gvcfs]
            }
            .dump(tag: "gvcfs_ConcatVcf")

    ConcatVCF(
                gvcf_ConcatVCF,
                ch_fasta_fai,
                ch_target_bed,
                'HC', // prefix for output files
                'g.vcf', // extension for the output files
                'HC_individually_genotyped_gvcf' // output directory name
                )

    // individually genotype gvcfs
    wf_individually_genotype_gvcf(
        wf_haplotypecaller.out.gvcf_GenotypeGVCFs,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_index,
        ch_target_bed
    )

    wf_genotype_gvcf(
        wf_haplotypecaller.out.gvcf_GenotypeGVCFs,
        ch_target_bed,
         ch_fasta,
        ch_fasta_fai,         
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_index
    )

    wf_hap_py(
        wf_deepvariant.out.vcf,
        wf_genotype_gvcf.out.vcfs_with_indexes,
        ch_target_bed,
        ch_bait_bed,
          ch_fasta,
        ch_fasta_fai
    )
    
    wf_germline_cnv(
        wf_recal_bams.out.bam_recal,
        ch_target_bed,
        ch_fasta,
        ch_fasta_fai
    )

    wf_gatk_somatic_cnv (
        wf_mark_duplicates.out.dm_bams,
        ch_target_bed,
        ch_fasta,
        ch_fasta_fai,         
        ch_dict,
        ch_read_count_pon
    )

    wf_savvy_somatic_cnv(wf_mark_duplicates.out.dm_bams)
    wf_vcf_stats(wf_deepvariant.out.vcf,
        wf_genotype_gvcf.out.vcfs_with_indexes
        )

    wf_multiqc(
        wf_get_software_versions.out,
        FastQCFQ.out.collect().ifEmpty([]),
        wf_qc_recal_bams.out.bam_qc,
        wf_vcf_stats.out.bcfootls_stats,
        wf_vcf_stats.out.vcfootls_stats,
        // wf_mark_duplicates.out.dm_bam_stats,
        wf_qc_recal_bams.out.samtools_stats,
        wf_qc_recal_bams.out.alignment_summary_metrics,
        wf_qc_recal_bams.out.insert_size_metrics, 
        wf_qc_recal_bams.out.hs_metrics.ifEmpty([])
    )
//    wf_annotate()
//    wf_multiqc()
} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}


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
        --target_bed                Target BED file (aka primary/regions/empirical) for targeted  or whole exome sequencing
        --bait_bed                  Bait BED file (aka covered/captured) for targeted or whole exome sequencing (used for GATK CollectHsMetrics)
        --step                      Specify starting step
                                    Available: Mapping, Recalibrate, VariantCalling, Annotate
                                    Default: Mapping
        --tools                     Specify tools to use for variant calling:
                                    Available: ASCAT, ControlFREEC, FreeBayes, HaplotypeCaller, DeepVariant, 
                                    Manta, mpileup, Mutect2, Mutect2_Single, gen_somatic_pon, gen_read_count_pon, 
                                    Strelka, TIDDIT
                                    and/or for annotation:
                                    snpEff, VEP, merge
                                    and for pipline validation (if you have added GIAB samples):
                                    hap_py
                                    Default: None
        --skip_qc                   Specify which QC tools to skip when running the pipeline
                                    Available: all, bamQC, BCFtools, FastQC, MultiQC, samtools, vcftools, versions
                                    Default: None
        --annotate_tools            Specify from which tools the pipeline will look for VCF files to annotate, only for step annotate
                                    Available: HaplotypeCaller, Manta, Mutect2, Strelka, TIDDIT
                                    Default: None
                                    Adds the following tools for --tools: DNAseq, DNAscope and TNscope
        --annotation_cache          Enable the use of cache for annotation, to be used with --snpEff_cache and/or --vep_cache        
        --pon_somatic               panel-of-normals VCF (bgzipped, indexed). See: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
        --pon_somatic_index         index of pon panel-of-normals VCF
        --pon_read_count            panel-of-normals hdf5 file. See https://gatk.broadinstitute.org/hc/en-us/articles/360040510031-CreateReadCountPanelOfNormals
        --filter_bams               Generate additional filter Bams based upon the bam_mapping_q parameter
        --bam_mapping_q             Specify the lower threshold for mapping quality (defaults to 60). All reads below this threshold will be discarded
                                    when generate the additional bams (you must specify the param --filter-bams)
        --remove_supplementary_reads When specified in addition to --filter_bams, will appy the following 
                                    sambamba expression to remove the supplementary reads 'not ([XA] != null or [SA] != null)'
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
        --known_indels              knownIndels file
        --known_indels_index        knownIndels index
                                    If none provided, will be generated automatically if a knownIndels file is provided
        --snpEff_db                 snpeffDb version
        --snpEff_cache              snpeffDb cache path if you have downloaded it locally
        --species                   species for VEP
        --vep_cache                 Path to VEP cache if you have downloaded it locally
        --vep_cache_version         VEP Cache version
    Other options:
        --outdir                    The output directory where the results will be saved
        --sequencing_center         Name of sequencing center to be displayed in BAM file
        --multiqc_config            Specify a custom config file for MultiQC
        --monochrome_logs           Logs will be without colors
        --email                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --max_multiqc_email_file_size   Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        --exome                     Specify when dealing with WES data, used when calling germline variants using Google deepvariant
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
 * Parse software version numbers
 */
process GetSoftwareVersions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    output:
        // file 'software_versions_mqc.yaml', emit: yamlSoftwareVersion
        file 'software_versions_mqc.yaml'

    when: !('versions' in skipQC)

    script:
    """
    bcftools version > v_bcftools.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    configManta.py --version > v_manta.txt 2>&1 || true
    configureStrelkaGermlineWorkflow.py --version > v_strelka.txt 2>&1 || true
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    echo "SNPEFF version"\$(snpEff -h 2>&1) > v_snpeff.txt
    fastqc --version > v_fastqc.txt 2>&1 || true
    freebayes --version > v_freebayes.txt 2>&1 || true
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    R --version &> v_r.txt  || true
    samtools --version &> v_samtools.txt 2>&1 || true
    tiddit &> v_tiddit.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built
process BuildFastaGz {
      tag "${fasta}.gz"
    //   publishDir "$baseDir/sampleDerivatives"

      input:
      file(fasta)

      output:
      file("${fasta}.gz")
      
      when: !(params.fasta_gz)
      script:
      """
      bgzip -c ${fasta} > ${fasta}.gz
      """
}

process BuildFastaGzFai {
    tag "${fasta}.gz.fai"
    // publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta)
    file(fastagz)

    output:
    file("${fasta}.gz.fai")
    when: !(params.fasta_gz_fai)
    script:
    """
    samtools faidx $fastagz
    """
  }
  
process BuildFastaGzi {
    tag "${fasta}.gz.gzi"
    // publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta)

    output:
    file("${fasta}.gz.gzi")
    
    when: !(params.fasta_gzi)

    script:
    """
    bgzip -c -i ${fasta} > ${fasta}.gz
    """
  }

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


process BuildSomaticPonIndex {
    tag {pon}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(pon)

    output:
        file("${pon}.tbi")

    when: !(params.somatic_pon_index) && params.somatic_pon && ('tnscope' in tools || 'mutect2' in tools)

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


process PartitionFastQ {
    // label 'PartitionFastQ'
    label 'cpus_32'

    tag {idPatient + "-" + idRun}

    // publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    // mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        tuple idPatient, idSample, idRun, file("r1_split_*"), file("r2_split_*") 
        // val (tuple_list)

    when: params.split_fastq && step == 'mapping'
    
    script:
    """
    partition.sh \
                in=${idSample}_${idRun}_R1.fastq.gz \
                in2=${idSample}_${idRun}_R2.fastq.gz  \
                out=r1_split_%.fastq.gz \
                out2=r2_split_%.fastq.gz \
                ways=${params.split_fastq}
    """
}

process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip}")

    when: !('fastqc' in skipQC) && (step == 'mapping')
    
    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}



// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM
process MapReads {
    label 'cpus_max'

    tag {idPatient + "-" + idRun}

    input:
        tuple idPatient, idSample, idRun, file(inputFile1), file(inputFile2)
        file(fasta) 
        file(fastaFai)
        file(bwaIndex) 

    output:
        // tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        // tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam")
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam"), emit : bam_mapped
        tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam"), emit : bam_mapped_BamQC
    
    when: !(step in ['recalibrate', 'variantcalling', 'annotate'])
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
// STEP 1.4 FILTERING BAMS ON THE BASIS OF QUALITY SCORES AND UNIQUELY MAPPED READS. ORIGINAL COMMANDS CAME FROM REBECCA BERNARD

process FilterBamRead1 {
    label 'cpus_32'

    tag {idPatient + "-" + idSample + "_" + idRun}

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        
    output:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_filtered_r1.bam"), emit: filtered_reads

    when: params.filter_bams
    script:
    if( params.remove_supplementary_reads)
        """
        sambamba view -t ${task.cpus} -h \
            -F "(first_of_pair and mapping_quality >=${params.bam_mapping_q} \
                and not ([XA] != null or [SA] != null)) \
                or second_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r1.bam"
        """

    else
        """
        sambamba view -t ${task.cpus} -h \
            -F "(first_of_pair and mapping_quality >=${params.bam_mapping_q}) \
                or second_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r1.bam"
        """
}

process FilterBamRead2 {
    label 'cpus_32'

    tag {idPatient + "-" + idSample + "_" + idRun}

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        
    output:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_filtered_r2.bam"), emit: filtered_reads

    when: params.filter_bams
    script:
    if( params.remove_supplementary_reads)
        """
        sambamba view -t ${task.cpus} -h \
            -F "(second_of_pair and mapping_quality >=${params.bam_mapping_q} \
                and not ([XA] != null or [SA] != null)) \
                or first_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r2.bam"
        """
    else
        """
        sambamba view -t ${task.cpus} -h \
            -F "(second_of_pair and mapping_quality >=${params.bam_mapping_q}) \
                or first_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r2.bam"
        """
}

process MergeFilteredBamReads {
    label 'cpus_32'
    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, idRun, file(partial_filtered_bams)

    output:
        tuple idPatient, idSample, idRun,  file(out_bam), file("${out_bam}.bai"), emit: filtered_bam

    when: (params.filter_bams)

    script:
     out_bam = "${idSample}_${idRun}_filtered.bam"
    // bn = file(out_bam).baseName
    // out_bai = "${bn}.bai"
    """
    samtools merge --threads ${task.cpus} -n -c -p merged.bam ${partial_filtered_bams}
    samtools sort --threads ${task.cpus} merged.bam -o ${idSample}_${idRun}_filtered.bam
    #samtools index  --threads ${task.cpus} ${idSample}_${idRun}_filtered.bam
    samtools index   ${idSample}_${idRun}_filtered.bam
    # cleaning
    rm -f merged.bam
    """
}

process FilterOutSecondaryAndSupplementaryReads {
    label 'cpus_32'
    tag {idPatient + "-" + idSample}
    input:
        tuple idPatient, idSample, idRun,  file("${idSample}_${idRun}_filtered.bam"), file(bai)

    output:
        tuple idPatient, idSample, idRun, file(out_bam),  emit: filtered_bam

    when: (params.filter_bams)

    script:
    out_bam = "${idSample}_${idRun}_pq${params.bam_mapping_q}.bam"
    bn = file(out_bam).baseName
    out_bai = "${bn}.bai"
    """
    #!/usr/bin/env python
    import pysam
    import shutil
    import os
    in_bam = "${idSample}_${idRun}_filtered.bam"
    # pysam.index("--threads", ${task.cpus}, in_bam)
    #pysam.index(in_bam)
    
    fd_in_bam = pysam.AlignmentFile(in_bam, "rb")
    fd_out_bam = pysam.AlignmentFile("uniq.bam", "wb", template=fd_in_bam)
    
    # Store reads to the dictionary so only unique reads stay
    reads_dict={}
    for r in fd_in_bam.fetch():
        reads_dict[(r.query_name, r.flag)] = r
    
    # iterate through the dict to write to the output bam file
    for key, read in reads_dict.items():
        fd_out_bam.write(read)

    fd_in_bam.close()
    fd_out_bam.close()
    # Sort the file by coordinates
    pysam.sort("-o", "${out_bam}", "--threads", "${task.cpus}", "uniq.bam")
    pysam.index("${out_bam}")
    shutil.move("${out_bam}.bai", "${out_bai}")
    # Cleanup
    os.unlink("uniq.bam")
    """
}


// STEP 1.5: MERGING BAM FROM MULTIPLE LANES
process MergeBamMapped {
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, idRun, out_suffix, file(bams)
        // tuple idPatient, idSample, idRun, file(bam), val bam_type

    output:
        tuple idPatient, idSample,  file("${idSample}${out_suffix}")

    script:
    // suffix = bams.first().minus("${idSample}_${idRun}")
    // out_file = "${idSample}${suffix}"
    """
    samtools merge --threads ${task.cpus} "${idSample}${out_suffix}" ${bams}
    """
}

process IndexBamFile {
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/Bams/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file(bam), file("${bam.baseName}.bai")

    // when: !params.knownIndels

    script:
    """
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    """
}

// STEP 2: MARKING DUPLICATES

// process Sambamba_MD {
//     label 'cpus_32'
//     // label 'memory_max'
//     // cache false
//     tag {idPatient + "-" + idSample}

//     publishDir "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked/", mode: params.publish_dir_mode

//     input:
//         tuple idPatient, idSample, file("${idSample}.bam")

//     output:
//         tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai")

//     """
//     sambamba markdup -t 20 -p ${idSample}.bam ${idSample}.md.bam
//     samtools index ${idSample}.md.bam
//     mv ${idSample}.md.bam.bai ${idSample}.md.bai
//     """
// }

// process MarkDuplicatesGATK {
//     label 'cpus_max'
//     // label 'memory_max'
//     // cache false
//     tag {idPatient + "-" + idSample}

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "${idSample}.bam.metrics" && 'markduplicates' in skipQC) null
//             else if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
//             else "Preprocessing/${idSample}/DuplicateMarked/${it}"
//         }

//     input:
//         tuple idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bai")
//         // tuple idPatient, idSample, file("${idSample}.bam")

//     output:
//         tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), emit: marked_bams
//         // tuple file ("${idSample}.bam.metrics"), emit: bam_stats
//         file ("${idSample}.bam.metrics")

//     when: params.known_indels

//     script:
//     // markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
//     markdup_java_options =  "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
//     // markdup_java_options =  "\"-Xms" + 16 + "g -Xmx" + 32 + "g\""
//     // markdup_java_options =  "\"-Xms16g -Xmx32g\""
//     """
//     gatk --java-options ${markdup_java_options} \
//         MarkDuplicates \
//         --MAX_RECORDS_IN_RAM 2000000 \
//         --INPUT ${idSample}.bam \
//         --METRICS_FILE ${idSample}.bam.metrics \
//         --TMP_DIR . \
//         --ASSUME_SORT_ORDER coordinate \
//         --CREATE_INDEX true \
//         --OUTPUT ${idSample}.md.bam
//     """
// }



process MarkDuplicates {
    label 'cpus_max'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bai")

    output:
        tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), emit: marked_bams
        // file ("${idSample}.bam.metrics")

    when: !(step in ['recalibrate', 'variantcalling', 'annotate'])

    script:
    """
    samtools sort -n --threads ${task.cpus}  -O SAM  ${idSample}.bam | \
        samblaster -M --ignoreUnmated| \
        samtools sort --threads ${task.cpus}  -O BAM > ${idSample}.md.bam

    samtools index ${idSample}.md.bam && \
        mv ${idSample}.md.bam.bai ${idSample}.md.bai
    """
}

process BaseRecalibrator {
    // label 'cpus_1'
    label 'cpus_8'
    // cache false
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed)
        file(fasta) 
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex) 
        file(knownIndels)
        file(knownIndelsIndex)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.table")
        // set idPatient, idSample into recalTableTSVnoInt

    when: params.known_indels  
        // ('haplotypecaller' in tools || 
        // 'mutect2' in tools ||
        // 'mutect2_single' in tools ||
        // 'gen_somatic_pon' in tools)

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
    label 'cpus_8'
    echo true
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked", mode: params.publish_dir_mode, overwrite: false

    input:
        tuple idPatient, idSample, file(recal)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.table"), emit: recal_table
        // set idPatient, idSample into recalTableTSV

    when: params.known_indels
        // ('haplotypecaller' in tools || 
        //     'mutect2' in tools ||
        //     'mutect2_single' in tools ||
        //     'gen_somatic_pon' in tools)

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
    label 'cpus_8'
    // label 'cpus_32'
    // label 'memory_max'
    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}
    // tag {idPatient + "-" + idSample }

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed)
        file(fasta)
        file(fastaFai) 
        file(dict)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.recal.bam")

    when: params.known_indels
        // ('haplotypecaller' in tools || 
        // 'mutect2' in tools ||
        // 'mutect2_single' in tools ||
        // 'gen_somatic_pon' in tools)

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

process BamRecalOnTarget {
    label 'cpus_32'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/RecalibratedOnTarget", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai")
        file(paddedTargetBed)

    output:
        tuple idPatient, idSample, file("${idSample}.recal.on_target.bam"), file("${idSample}.recal.on_target.bam.bai"), emit: bam_recal_on_target
        tuple idPatient, idSample, file("${idSample}.recal.on_target.bam"), emit: bam_recal_on_target_qc
        
    when: params.padded_target_bed && !('on_target_assessment' in skipQC)

    script:
    """
    bedtools intersect -a ${idSample}.recal.bam -b ${paddedTargetBed} > ${idSample}.recal.on_target.bam
    samtools index ${idSample}.recal.on_target.bam
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
    // cache false

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

process CollectAlignmentSummaryMetrics{
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/alignment_summary/", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam) 
    file(fasta) 
    file(fastaFai)
    file(dict)

    output:
    file("${bam.baseName}_alignment_metrics.txt")
    
    when: ! ('alignment_summary' in skipQC)
    
    script:
    """
    gatk --java-options -Xmx32G CollectAlignmentSummaryMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}_alignment_metrics.txt \
    -R ${fasta}
    """
}

process CollectInsertSizeMetrics{
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/insert_size_metrics/", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam)

    output:
    file("${bam.baseName}_insert_size_metrics.txt")
    file("${bam.baseName}_insert_size_histogram.pdf")
    
    
    when: !('insert_size_metrics' in skipQC)

    script:
    """
    gatk --java-options -Xmx32G CollectInsertSizeMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}_insert_size_metrics.txt \
    -H ${bam.baseName}_insert_size_histogram.pdf 
    """
}

process CollectHsMetrics{
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Reports/${idSample}/hs_metrics/", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam)
    file(fasta) 
    file(fastaFai)
    file(dict)
    file(targetBED)
    file(baitBED)
    // val (output_suffix)

    output:
    file("${bam.baseName}.txt")
    // file("${bam.baseName}_${output_suffix}.txt")
    
    
    when: !('hs_metrics' in skipQC) && params.bait_bed
    script:
    """
    gatk BedToIntervalList -I ${targetBED} -O target.interval_list -SD ${dict}
    gatk BedToIntervalList -I ${baitBED} -O bait.interval_list -SD ${dict}

    gatk --java-options -Xmx32G CollectHsMetrics --VALIDATION_STRINGENCY LENIENT \
    -I ${bam} \
    -O ${bam.baseName}.txt \
    -TI target.interval_list \
    -BI bait.interval_list \
    -R ${fasta}
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
    label 'cpus_8'
    
    tag {idSample + "-" + intervalBed.baseName}
    // tag {idSample} 
    // publishDir "${params.outdir}/VariantCalling/${idSample}/HaplotypeCaller", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed) 
        // tuple idPatient, idSample, file(bam), file(bai) 
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)

    output:
        tuple val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_HC
        // tuple idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_HaplotypeCaller
        tuple idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_GenotypeGVCFs
        // tuple val("${intervalBed.baseName}"), idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_GenotypeGVCFs
        

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

// process GvcfToVcf{
//     label 'memory_singleCPU_task_sq'
//     label 'cpus_2'
//     // label 'memory_max'
//     // label 'cpus_max'

//     tag {idSample + "-" + gvcf.baseName}
//     // tag {idSample} 
//     publishDir "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_vcf", mode: params.publish_dir_mode
//     input:
//         tuple val(variantCaller), idPatient, idSample, file(gvcf)

//     output:
//         // tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file("${gvcf.simpleName}.vcf"), emit: vcf_HaplotypeCaller
//         tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file(out_file), emit: vcf_HaplotypeCaller

//     when: 'haplotypecaller' in tools

//     script:
//     // fn=gvcf.fileName
//     // prefix=fn.minus(".g.vcf")
//     // out_file="${gvcf.fileName}.vcf"
//     prefix="${gvcf.fileName}" - ".g.vcf.gz"
//     // We'll first decompress the gvcf
//     in_file= "${gvcf.fileName}" - ".gz"
//     out_file="${prefix}.vcf"
//     """
//     gzip -d --force ${gvcf}
//     extract_variants < ${in_file} > ${out_file}
//     """
// }

process IndividuallyGentoypeGVCF{
    label 'cpus_8'
    tag {idSample + "-" + gvcf.baseName}
    // tag {idSample} 
    // publishDir "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_vcf", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(intervalsBed), file(gvcf)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)
    output:
        // tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file("${gvcf.simpleName}.vcf"), emit: vcf_HaplotypeCaller
        tuple idPatient, idSample, file(out_file), emit: vcf_HaplotypeCaller

    when: 'haplotypecaller' in tools

    script:
    // fn=gvcf.fileName
    // prefix=fn.minus(".g.vcf")
    // out_file="${gvcf.fileName}.vcf"
    prefix="${gvcf.fileName}" - ".g.vcf"
    // We'll first decompress the gvcf
    // in_file= "${gvcf.fileName}" - ".gz"
    out_file="${prefix}.vcf"
    """
    bgzip  ${gvcf}
    tabix  ${gvcf}.gz
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${intervalsBed} \
        -D ${dbsnp} \
        -V ${gvcf}.gz \
        -O "${out_file}"
    """
}

// STEP GATK HAPLOTYPECALLER.1.5
process GenomicsDBImport {
    label 'cpus_16'
    // echo true
    tag{interval_name}
    // publishDir "${OUT_DIR}/misc/genomicsdb/", mode: 'copy', overwrite: false

    input:
    // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gvcfs)
    tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file(gvcfs)
    
    output:
    tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file ("${interval_name}.gdb")

    when: 'haplotypecaller' in tools

    script:
    sample_map="cohort_samples.map"
    interval_name_with_underscore="${interval_name}_"
    // gDB = chr
    """
    for x in *.g.vcf
    do
        bgzip \$x
        tabix \${x}.gz
    done

    for x in *.g.vcf.gz
    do
        
        base_name=`basename \$x .g.vcf.gz`
        sample=\${base_name#$interval_name_with_underscore}
        echo "\${sample}\t\${x}" >> ${sample_map}
    done
    
    gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${interval_name}.gdb \
    -L $interval_bed \
    --sample-name-map ${sample_map} \
    --reader-threads ${task.cpus}

    """
}

// STEP GATK HAPLOTYPECALLER.2


process GenotypeGVCFs {
    label 'cpus_8'
    tag {interval_bed.baseName}
    input:
        // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gdb)
        tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file(gdb)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)

    output:
    // tuple val("HaplotypeCaller"), list_id_patient, list_id_sample, file("${interval_name}.vcf"), file("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    tuple val("HaplotypeCaller"),  val(patientSampleIdMap), file(interval_bed), file("${interval_name}.vcf"), file("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    // tuple val("HaplotypeCaller_Jointly_Genotyped"), val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file ("${interval_name}.vcf"), file ("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    
    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${interval_bed} \
        -D ${dbsnp} \
        -V gendb://${gdb} \
        --create-output-variant-index \
        -O "${interval_name}.vcf"
    """
}

process SelectVariants {
    label 'cpus_8'
    tag {interval_bed.baseName}
    input:
        // tuple val(caller), val(id_patient), val(id_sample), val(interval_name), file(interval_bed), file (vcf), file (vcf_idx)
        tuple val(caller), val(id_patient), val(id_sample), file(interval_bed), file (vcf), file (vcf_idx)
        file(fasta)
        file(fastaFai)
        file(dict)

    output:
    // tuple val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf"), emit: vcf_GenotypeGVCFs
    // tuple val("HaplotypeCaller"), val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file ("${interval_name}.vcf"), file ("${interval_name}.vcf.idx"), emit: vcf_GenotypeGVCFs
    tuple val("HaplotypeCaller_Jointly_Genotyped"), id_patient, id_sample, file("${interval_bed.baseName}_${id_sample}.vcf"), emit: vcf_SelectVariants
    
    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
            SelectVariants \
            -R ${fasta} \
            -L ${interval_bed} \
            -V ${vcf} \
            -O ${interval_bed.baseName}_${id_sample}.vcf \
            -sn ${id_sample}
    """
}

process ConcatVCF {
    label 'cpus_8'

    tag {variantCaller + "-" + idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/${output_dir}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idPatient, idSample, file(vcFiles)
        file(fastaFai)
        file(targetBED)
        val(output_file_prefix)
        val(output_file_ext)
        val(output_dir)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        // tuple variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi"), emit: concatenated_vcf_with_index
        // tuple variantCaller, idPatient, idSample, file("*_*.vcf.gz"), emit: concatenated_vcf_without_index
        tuple variantCaller, idPatient, idSample, file("${outFile}.gz"), 
            file("${outFile}.gz.tbi"), emit: concatenated_vcf_with_index

        tuple variantCaller, idPatient, idSample, file("${outFile}.gz"), 
            emit: concatenated_vcf_without_index
            
    script:
    outFile =  "${output_file_prefix}_${idSample}.${output_file_ext}"
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outFile} ${options}
    """
}

process CohortConcatVCF {
    label 'cpus_8'

    tag {'CohortConcatVCF'}

    publishDir "${params.outdir}/VariantCalling/HC_cohort_vcf", mode: params.publish_dir_mode

    input:
        file(vcFiles)
        file(fastaFai)
        file(targetBED)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        tuple file("HC_cohort.vcf.gz"), file("HC_cohort.vcf.gz.tbi"), emit: cohort_vcf_with_index
        // tuple file("HC_cohort.vcf.gz"), emit: cohort_vcf_without_index
        file("HC_cohort.vcf.gz")
        // file("HC_cohort.vcf.gz"), file("HC_cohort.vcf.gz.tbi"), emit: cohortvcfwithindex
        // file("HC_cohort.vcf.gz"), emit: cohortvcfwithoutindex

    when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o HC_cohort.vcf ${options}
    """
}


/*
================================================================================
                               Pipeline Benchmarking using Illumina Hap.py
================================================================================
*/
process HapPy {
    label 'cpus_32'

    tag {idSample}

    publishDir "${params.outdir}/Validation/Against_${truth_set_type}/${idSample}/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple truth_set_type, variantCaller, idPatient, idSample, file(vcf), file(tbi), file (truth_set_vcf), file (truth_set_tbi), file(highqual_regions)         
        file(target_bed)
        file(bait_bed)
        file(fasta)
        file(fastaFai)

    output:
        file("Against_${truth_set_type}.${vcf.baseName}.*")

    when: ( 
            ('benchmark_dv_and_hc_against_giab' in tools || 'benchmark_dv_against_hc' in tools )
            && ('haplotypecaller' in tools || 'deepvariant' in tools) 
            && params.giab_highconf_vcf 
            && params.giab_highconf_tbi 
            && params.giab_highconf_regions
        )

    script:
    // bn = "{vcf.baseName}"
    """
    export HGREF=$fasta
    mkdir scratch
    hap.py  \
        ${truth_set_vcf} \
        ${vcf} \
        -f ${highqual_regions} \
        --scratch-prefix scratch \
        --engine vcfeval \
        -T ${target_bed} \
        --threads ${task.cpus} \
        -o Against_${truth_set_type}.${vcf.baseName}
    """
}

/*
================================================================================
                               Google Deepvariant  for Germline variant calling
================================================================================
*/

/********************************************************************
  process make_examples
  Getting bam files and converting them to images ( named examples )
********************************************************************/

process DV_Combined{
    label "cpus_max"
    tag "${bam}"
    publishDir "${params.outdir}/VariantCalling/${idSample}/DeepVariant", mode: params.publish_dir_mode
    
    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file (fasta)
    file (fai )
    file (fastagz)
    file (gzfai)
    file (gzi)
    file (target_bed)
    
    output:
        tuple idPatient, idSample, file("${bam.baseName}.vcf")
        file("*.html")

    when: 'deepvariant' in tools
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${fasta} \
    --reads ${bam} \
     --regions ${target_bed} \
    --output_vcf "${bam.baseName}.vcf" \
    --num_shards ${task.cpus}
    """
}

process DV_MakeExamples{
    label 'cpus_max'
    tag "${bam}"
    publishDir "${params.outdir}/Preprocessing/${idSample}/DV_MakeExamples/", mode: params.publish_dir_mode,
    saveAs: {filename -> "logs/log"}

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file (fasta)
    file (fai )
    file (fastagz)
    file (gzfai)
    file (gzi)
    file (target_bed)

    output:
    tuple idPatient, idSample, file(bam), file('*_shardedExamples'), emit: shared_examples
    when: 'deepvariant' in tools
    script:
    // sharded_folder = "${bam}.tfrecord@${task.cpus}.gz"
    use_regions = params.target_bed ? "--regions ${target_bed}" :  ''
    """
    mkdir logs
    mkdir ${bam.baseName}_shardedExamples
    dv_make_examples.py \
    --cores ${task.cpus} \
    --sample ${bam} \
    --ref ${fastagz} \
    --reads ${bam} \
    ${use_regions} \
    --logdir logs \
    --examples ${bam.baseName}_shardedExamples

    """
}

/********************************************************************
  process call_variants
  Doing the variant calling based on the ML trained model.
********************************************************************/

process DV_CallVariants{
   label 'cpus_max'
  tag "${bam}"

  input:
  tuple idPatient, idSample, file(bam), file(shardedExamples)

  output:
  tuple idPatient, idSample, file(bam), file('*_call_variants_output.tfrecord'), emit: variant_tf_records
  

  when: 'deepvariant' in tools
  script:
  """
    dv_call_variants.py \
    --cores ${task.cpus} \
    --sample ${bam} \
    --outfile ${bam.baseName}_call_variants_output.tfrecord \
    --examples $shardedExamples \
    --model ${model}

  """
}



/********************************************************************
  process postprocess_variants
  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
********************************************************************/

process DV_PostprocessVariants{

  label 'cpus_32'
  tag "${bam}"
  publishDir "${params.outdir}/VariantCalling/${idSample}/DeepVariant", mode: params.publish_dir_mode

  input:
  tuple idPatient, idSample, file(bam), file(call_variants_tfrecord)
//   set file(bam),file('call_variants_output.tfrecord') from called_variants
  file fastagz
  file gzfai
  file gzi

  output:
   tuple val('DeepVariant'), idPatient, idSample, file("${idSample}.vcf.gz"), file("${idSample}.vcf.gz.tbi") , emit: vcf
   tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz") , emit: vcf_to_annotate
//    file("${idSample}.*.html"), emit: html_report
   file('*.html')

  when: 'deepvariant' in tools
  script:
  """
  dv_postprocess_variants.py \
  --ref ${fastagz} \
  --infile ${call_variants_tfrecord} \
  --outfile "${idSample}.vcf" \

  bgzip "${idSample}.vcf"
  tabix -p vcf "${idSample}.vcf.gz"
  """
}

/*
================================================================================
                               Copy Number callers
================================================================================
*/


// STEP MPILEUP.1

process Mpileup {
    label 'memory_singleCPU_2_task'
    tag {idSample + "-" + intervalBed.baseName}
    
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup.gz" ? "VariantCalling/${idSample}/mpileup/${it}" : '' }

    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed)
        file(fasta)
        file(fastaFai)

    output:
        tuple idPatient, idSample, file("${prefix}${idSample}.pileup.gz")

    when: 'controlfreec' in tools || 'mpileup' in tools

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-l ${intervalBed}"
    """
    samtools mpileup \
        -f ${fasta} ${bam} \
        ${intervalsOptions} \
    | bgzip --threads ${task.cpus} -c > ${prefix}${idSample}.pileup.gz
    """
}


// STEP MPILEUP.2 - MERGE

process MergeMpileup {
    tag {idSample}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup.gz" ? "VariantCalling/${idSample}/mpileup/${it}" : '' }

    input:
        tuple idPatient, idSample, file(mpileup)

    output:
        tuple idPatient, idSample, file("${idSample}.pileup.gz")

    when: !(params.no_intervals) && 'controlfreec' in tools || 'mpileup' in tools

    script:
    """
    for i in `ls -1v *.pileup.gz`;
        do zcat \$i >> ${idSample}.pileup
    done

    bgzip --threads ${task.cpus} -c ${idSample}.pileup > ${idSample}.pileup.gz

    rm ${idSample}.pileup
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

// STEP TIDDIT

process TIDDIT {
    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/TIDDIT", mode: params.publish_dir_mode

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta) 
        file(fastaFai)

    output:
        tuple val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi"), emit: vcfTIDDIT
        tuple file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig"), emit: tidditOut

    when: 'tiddit' in tools

    script:
    """
    tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

    tabix TIDDIT_${idSample}.vcf.gz
    """
}


// STEP VCF.QC

process BcftoolsStats {
    label 'cpus_1'

    tag {"${variantCaller} - ${vcf}"}

    publishDir "${params.outdir}/Reports/${idSample}/BCFToolsStats/${variantCaller}", mode: params.publish_dir_mode

    input:
        // tuple variantCaller, idSample, file(vcf)
        tuple variantCaller, idPatient, idSample, file(vcf) , file(vcf_tbi)

    output:
        // file ("*.bcf.tools.stats.out"), emit: bcftools_reports
        file ("*.bcf.tools.stats.out")

    when: !('bcftools' in skipQC)

    script:
    """
    bcftools stats ${vcf} > ${reduceVCF(vcf.fileName)}.bcf.tools.stats.out
    """
}


process Vcftools {
    label 'cpus_1'

    tag {"${variantCaller} - ${vcf}"}

    publishDir "${params.outdir}/Reports/${idSample}/VCFTools/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idPatient, idSample, file(vcf) , file(vcf_tbi)

    output:
        file ("${reduceVCF(vcf.fileName)}.*")

    when: !('vcftools' in skipQC)

    script:
    """
    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${reduceVCF(vcf.fileName)}
    """
}




/*
================================================================================
                                   ANNOTATION
================================================================================
*/

// STEP SNPEFF

process SnpEff {
    tag {"${idSample} - ${variantCaller} - ${vcf}"}
    // cache false
    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${variantCaller}/${it}"
    }

    input:
        tuple variantCaller, idSample, file(vcf) 
        file(dataDir)
        // path(dataDir)
        val snpeffDb

    output:
        tuple file("${reducedVCF}_snpEff.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv"), emit:snpEff_report
        tuple variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf"), emit: snpEff_vcf

    when: 'snpeff' in tools || 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    cache = (params.snpEff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
    // cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir ${dataDir}" : ""
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    mv ${reducedVCF}_snpEff.genes.txt ${reducedVCF}_snpEff.txt
    """
}

// snpeffReport = snpeffReport.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process CompressVCFsnpEff {
    tag {"${idSample} - ${vcf}"}

    publishDir "${params.outdir}/Annotation/${idSample}/snpEff/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idSample, file(vcf)

    output:
        tuple variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: compressVCFsnpEffOut

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}


// STEP VEP.1

process VEP {
    label 'VEP'
    label 'cpus_4'

    tag {"${idSample} - ${variantCaller} - ${vcf}"}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${variantCaller}/${it}"
        else null
    }

    input:
        tuple variantCaller, idSample, file(vcf), file(idx) 
        file(dataDir) 
        val cache_version 
        file(cadd_InDels) 
        file(cadd_InDels_tbi) 
        file(cadd_WG_SNVs) 
        file(cadd_WG_SNVs_tbi) 
        file(fasta)
        file(fasta_fai)
    output:
        tuple variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf"), emit: vep_vcf
        // file("${reducedVCF}_VEP.summary.html") , emit: vep_report
        file("${reducedVCF}_VEP.summary.html") 

    when: ('vep' in tools) && params.cadd_InDels && params.cadd_WG_SNVs

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    // cadd = (params.cadd_cache && params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    cadd = (params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/miniconda/envs/layer_lab_dna_seq/bin/genesplicer,/opt/miniconda/envs/layer_lab_dna_seq/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf \
        --offline
    rm -rf ${reducedVCF}
    """
}

// vepReport = vepReport.dump(tag:'VEP')

// STEP VEP.2 - VEP AFTER SNPEFF

process VEPmerge {
    label 'VEP'
    label 'cpus_4'

    tag {"${idSample} - ${variantCaller} - ${vcf}"}

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${variantCaller}/${it}"
        else null
    }

    input:
        tuple variantCaller, idSample, file(vcf), file(idx) 
        file(dataDir) 
        val cache_version 
        file(cadd_InDels) 
        file(cadd_InDels_tbi) 
        file(cadd_WG_SNVs) 
        file(cadd_WG_SNVs_tbi) 
        file(fasta)
        file(fasta_fai)
    output:
        tuple variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf"), emit: vep_vcf_merge
        // file("${reducedVCF}_VEP.summary.html") , emit: vep_report_merge
        file("${reducedVCF}_VEP.summary.html") 

    when: ('merge' in tools) && params.cadd_InDels && params.cadd_WG_SNVs

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    // cadd = (params.cadd_cache && params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    cadd = (params.cadd_WG_SNVs && params.cadd_InDels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/miniconda/envs/layer_lab_dna_seq/bin/genesplicer,/opt/miniconda/envs/layer_lab_dna_seq/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf \
        --offline

    rm -rf ${reducedVCF}
    """
}


// STEP COMPRESS AND INDEX VCF.2 - VEP

process CompressVCFvep {
    tag {"${idSample} - ${vcf}"}

    publishDir "${params.outdir}/Annotation/${idSample}/VEP/${variantCaller}", 
    mode: params.publish_dir_mode

    input:
        tuple variantCaller, idSample, file(vcf) 

    output:
        tuple variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") 

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

/*
================================================================================
                                     Mutect2Single (Single Sample Mode)
================================================================================
*/

process Mutect2Single{
    tag {idSample + "-" + intervalBed.baseName}
    label 'cpus_16'

    input:
        tuple idPatient, idSample, file(bam), file(bai), 
            file(intervalBed)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(germlineResource)
        file(germlineResourceIndex)

    output:
        tuple idPatient, idSample, file(out_vcf)

    script:
    out_vcf = "${intervalBed.baseName}_${idSample}.vcf"
    """
    # max-mnp-distance is set to 0 to avoid a bug in 
    # next process GenomicsDbImport
    # See https://gatk.broadinstitute.org/hc/en-us/articles/360046224491-CreateSomaticPanelOfNormals-BETA-
    
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta} \
      -I ${bam}  \
      -max-mnp-distance 0 \
      -L ${intervalBed} \
      --germline-resource ${germlineResource} \
      -O ${out_vcf}
    """
}

// STEP GATK GenomicsDBImport
process SomaticPonGenomicsDBImport {
    label 'cpus_32'

    publishDir "${params.outdir}/Preprocessing/Somatic_pon_db", mode: params.publish_dir_mode

    input:
    file("vcfs/*")
    file(targetBED)

    output:
    file("somatic_pon.gdb")

    when: 'gen_somatic_pon' in tools

    script:
    sample_map="cohort_samples.map"
    
    // gDB = chr
    """
    vcfs=' '
    for x in `ls vcfs/*.vcf.gz`
    do
        base_name=`basename \${x}`
        without_ext=\${base_name%.vcf.gz}
        sample_name=\${without_ext##*_}
        echo "\${sample_name}\t\$x" >> $sample_map 
    done

    gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenomicsDBImport  \
    --genomicsdb-workspace-path somatic_pon.gdb \
    -L ${targetBED} \
    --sample-name-map $sample_map \
    --merge-input-intervals \
    --reader-threads ${task.cpus}
    """
}

process CreateSomaticPON{
    label 'cpus_max'
    // label 'memory_max'
     publishDir "${params.outdir}/Preprocessing/Somatic_pon", mode: params.publish_dir_mode

    input:
    file(pon) 
    file(fasta)
    file(fastaFai)
    file(dict)
    file(germlineResource)
    file(germlineResourceIndex)
    
    output:
    tuple file(out_file), file ("${out_file}.tbi")

    when: 'gen_somatic_pon' in tools

    script:
    args_file = "normals_for_pon_vcf.args"
    out_file = "somatic_pon.vcf.gz" 
    pon_db = "gendb://${pon}"
    
    """
     gatk --java-options -Xmx${task.memory.toGiga()}g \
     CreateSomaticPanelOfNormals -R ${fasta} \
     --germline-resource ${germlineResource} \
    -V ${pon_db} \
    -O ${out_file}
    """
}


/*
================================================================================
                                     Mutect2 (Tumor/Normal Mode)
================================================================================
*/

// process Mutect2{
//     tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}
//     label 'cpus_2'

//     input:
//         tuple idPatient, 
//             idSampleNormal, file(bamNormal), file(baiNormal),
//             idSampleTumor, file(bamTumor), file(baiTumor), 
//             file(intervalBed)
//         file(fasta)
//         file(fastaFai)
//         file(dict)
//         file(germlineResource)
//         file(germlineResourceIndex)
//         file(ponSomatic)
//         file(ponSomaticIndex)

//     output:
//         tuple 
//             val("Mutect2"), 
//             idPatient,
//             val("${idSampleTumor}_vs_${idSampleNormal}"),
//             file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf")
        
//         tuple 
//             idPatient,
//             idSampleTumor,
//             idSampleNormal,
//             file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats") optional true

//     when: 'mutect2' in tools

//     script:
//     // please make a panel-of-normals, using at least 40 samples
//     // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
//     PON = params.pon_somatic ? "--panel-of-normals ${ponSomatic}" : ""
//     """
//     # Get raw calls
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//       Mutect2 \
//       -R ${fasta}\
//       -I ${bamTumor}  -tumor ${idSampleTumor} \
//       -I ${bamNormal} -normal ${idSampleNormal} \
//       -L ${intervalBed} \
//       --germline-resource ${germlineResource} \
//       ${PON} \
//       -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// // STEP GATK MUTECT2.2 - MERGING STATS

// process MergeMutect2Stats {
//     label 'cpus_16'
//     tag {idSampleTumor + "_vs_" + idSampleNormal}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Mutect2", mode: params.publishDirMode

//     input:
//         // tuple caller, idPatient, idSampleTumor_vs_idSampleNormal, file(vcfFiles) // corresponding small VCF chunks
//         tuple idPatient, idSampleTumor, idSampleNormal, file(statsFiles)// the actual stats files
//         // file(dict)
//         // file(fasta)
//         // file(fastaFai)
//         // file(germlineResource)
//         // file(germlineResourceIndex)
//         // file(intervals) from ch_intervals

//     output:
//         // tuple idPatient,
//         //     val("${idSampleTumor}_vs_${idSampleNormal}"),
//         file("${idSampleTumor_vs_idSampleNormal}.vcf.gz.stats")

//     when: 'mutect2' in tools

//     script:     
//       stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         MergeMutectStats \
//         ${stats} \
//         -O ${idSampleTumor}_vs_${idSampleNormal}.vcf.gz.stats
//     """
// } // end of MergeMutect2Stats

// // STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

// process PileupSummariesForMutect2 {
//     tag {idSampleTumor + "_vs_" + idSampleNormal + "_" + intervalBed.baseName }
//     label 'cpus_1'

//     input:
//         tuple idPatient, 
//             idSampleNormal, file(bamNormal), file(baiNormal), 
//             idSampleTumor, file(bamTumor), file(baiTumor), 
//             file(intervalBed)
        
//         // set idPatient, idSampleNormal, idSampleTumor, file(statsFile) from intervalStatsFiles
//         file(germlineResource)
//         file(germlineResourceIndex)

//     output:
//         tuple idPatient,
//             idSampleTumor,
//             file("${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table")

//     when: 'mutect2' in tools && params.pon

//     script:
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GetPileupSummaries \
//         -I ${bamTumor} \
//         -V ${germlineResource} \
//         -L ${intervalBed} \
//         -O ${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table
//     """
// } // end of PileupSummariesForMutect2

// // STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

// process MergePileupSummaries {
//     label 'cpus_1'

//     tag {idPatient + "_" + idSampleTumor}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

//     input:
//         tuple idPatient, idSampleTumor, file(pileupSums)
//         file(dict)

//     output:
//         file("${idSampleTumor}_pileupsummaries.table.tsv")

//     when: 'mutect2' in tools
//     script:
//         allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GatherPileupSummaries \
//         --sequence-dictionary ${dict} \
//         ${allPileups} \
//         -O ${idSampleTumor}_pileupsummaries.table.tsv
//     """
// }

// // STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

// process CalculateContamination {
//     label 'cpus_1'

//     tag {idSampleTumor + "_vs_" + idSampleNormal}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

//     input:
//         // tuple idPatient, 
//         //     idSampleNormal, file(bamNormal), file(baiNormal), 
//         //     idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCalculateContamination 
//         file("${idSampleTumor}_pileupsummaries.table")
  
//     output:
//         file("${idSampleTumor}_contamination.table") into contaminationTable

//     when: 'mutect2' in tools && params.pon

//     script:     
//     """
//     # calculate contamination
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         CalculateContamination \
//         -I ${idSampleTumor}_pileupsummaries.table \
//         -O ${idSampleTumor}_contamination.table
//     """
// }

// // STEP GATK MUTECT2.6 - FILTERING CALLS

// process FilterMutect2Calls {
//     label 'cpus_1'

//     tag {idSampleTN}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTN}/${variantCaller}", mode: params.publishDirMode

//     input:
//         // tuple variantCaller, 
//         //     idPatient, idSampleTN, file(unfiltered), file(unfilteredIndex),
//         //     file("${idSampleTN}.vcf.gz.stats"),
//         //     file("${idSampleTN}_contamination.table")
//         tuple variantCaller, 
//             idPatient, idSampleTN, file(unfiltered), file(unfilteredIndex)
//         file("${idSampleTN}.vcf.gz.stats")
//         file("${idSampleTN}_contamination.table")
        
//         file(dict)
//         file(fasta)
//         file(fastaFai)
//         file(germlineResource)
//         file(germlineResourceIndex)
//         // file(intervals) from ch_intervals
        
//     output:
//         tuple val("Mutect2"), idPatient, idSampleTN,
//             file("filtered_${variantCaller}_${idSampleTN}.vcf.gz"),
//             file("filtered_${variantCaller}_${idSampleTN}.vcf.gz.tbi"),
//             file("filtered_${variantCaller}_${idSampleTN}.vcf.gz.filteringStats.tsv")

//     // when: 'mutect2' in tools && params.pon
//     when: 'mutect2' in tools

//     script:
//     """
//     # do the actual filtering
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         FilterMutectCalls \
//         -V ${unfiltered} \
//         --contamination-table ${idSampleTN}_contamination.table \
//         --stats ${idSampleTN}.vcf.gz.stats \
//         -R ${fasta} \
//         -O filtered_${variantCaller}_${idSampleTN}.vcf.gz
//     """
// }


/*
================================================================================
                                     Somatic CNV
================================================================================
*/

process PreprocessIntervals {

    label 'cpus_8'
    
    input:    
        file(intervalBed)
        file(fasta)
        file(fasta_fai)
        file(dict)
    
    output:
        // file("preprocessed_intervals.interval_list"), emit: 'processed_intervals'
        file("preprocessed_intervals.interval_list")

    when: ('gatkcnv' in tools) || ('gen_read_count_pon' in tools)
    
    script:
    intervals_options = params.no_intervals ? "" : "-L ${intervalBed}"
    padding_options =  params.no_intervals ? "--padding 0" : "--padding 250"
    bin_options =  params.no_intervals ? "--bin-length 1000" : "--bin-length 0"

    """
    gatk PreprocessIntervals \
        ${intervals_options} \
        ${padding_options} \
        ${bin_options} \
        -R ${fasta} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O preprocessed_intervals.interval_list
    """
}

process CollectReadCounts {
    label 'cpus_32'
    tag "${idSample}"
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(preprocessed_intervals)

    output:
        tuple idPatient, idSample, file("${idSample}.counts.hdf5"), emit: 'sample_read_counts'

    when: ('gatkcnv' in tools) || ('gen_read_count_pon' in tools)

    script:
    """
    gatk CollectReadCounts \
        -I ${bam} \
        -L ${preprocessed_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${idSample}.counts.hdf5
    """
}

process CreateReadCountPon {
    // echo true
    tag "ReadCountPon"
    
    publishDir "${params.outdir}/Preprocessing/ReadCountPon/", 
    mode: params.publish_dir_mode

    
    input:
    file(read_count_hdf5s)
    
    // file(this_read)

    output:
    file(out_file)

    script:
    when:'gen_read_count_pon' in tools
    // sample = this_read.simpleName
    out_file = "read_count_pon.hdf5"
    params_str = ''
    // Only get the normal samples
    read_count_hdf5s.each{
        params_str = "${params_str} -I ${it}"
    }

    
    """
    gatk CreateReadCountPanelOfNormals \
        $params_str \
        -O $out_file
    """
}

process DenoiseReadCounts {
    label 'cpus_32'
    tag "${idSample}"
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/DenoisedReadCounts/", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file( "${idSample}.counts.hdf5")
        file(read_count_somatic_pon)

    output:
        tuple idPatient, idSample, file(std_copy_ratio), file(denoised_copy_ratio), emit: 'denoised_cr'

    when: 'gatkcnv' in tools

    script:
    std_copy_ratio = "${idSample}.standardizedCR.tsv"
    denoised_copy_ratio = "${idSample}.denoisedCR.tsv"
    pon_option = params.read_count_pon ? "--count-panel-of-normals ${read_count_somatic_pon}" : ""
    """
    
    gatk DenoiseReadCounts \
        -I ${idSample}.counts.hdf5 \
        ${pon_option} \
        --standardized-copy-ratios ${std_copy_ratio} \
        --denoised-copy-ratios ${denoised_copy_ratio}
    """
}

process PlotDenoisedCopyRatios {
    label 'cpus_16'
    tag "${idSample}"
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(std_copy_ratio), file(denoised_copy_ratio)
        file(dict)
    
    output:
        file(out_dir)  

    when: 'gatkcnv' in tools

    script:
    out_dir = "PlotDenoisedReadCounts" 

    """
    mkdir ${out_dir}
    gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios ${std_copy_ratio} \
        --denoised-copy-ratios ${denoised_copy_ratio} \
        --sequence-dictionary ${dict} \
        --output-prefix ${idSample} \
        -O ${out_dir}
    """
}

process ModelSegments {
    label 'cpus_32'
    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}", mode: params.publish_dir_mode
    
    input:
         tuple idPatient, idSample, file(denoised_copy_ratio)

    output:
        tuple idPatient, idSample, file("${out_dir}/${idSample}.cr.seg"), file("${out_dir}/${idSample}.modelFinal.seg"), emit: 'modeled_seg'

    when: 'gatkcnv' in tools

    script:
    out_dir = "ModeledSegments"

    """
    mkdir $out_dir
    gatk ModelSegments \
        --denoised-copy-ratios ${denoised_copy_ratio} \
        --output-prefix ${idSample} \
        -O ${out_dir}
    """
}

process PlotModeledSegments {
    label 'cpus_8'
    tag "${idSample}"
    
    publishDir "${params.outdir}/VariantCalling/${idSample}", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file("${idSample}.modelFinal.seg"), file("${idSample}.denoisedCR.tsv")
        file(dict)
    output:
    file(out_dir)
    
    when: 'gatkcnv' in tools
    script:
    out_dir = "PlotsModeledSegments"
    
    """
    mkdir $out_dir
    gatk PlotModeledSegments \
        --denoised-copy-ratios ${idSample}.denoisedCR.tsv \
        --segments ${idSample}.modelFinal.seg \
        --sequence-dictionary ${dict} \
        --output-prefix ${idSample} \
        -O $out_dir
    """
}

process CallCopyRatioSegments {
   label 'cpus_8'
    tag "${idSample}"
    
    publishDir "${params.outdir}/VariantCalling/${idSample}/CalledCopyRatioSegments", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file("${idSample}.cr.seg")
    
    output:
        file("${idSample}.called.seg")

    when: 'gatkcnv' in tools
    script:
    
    """
    gatk CallCopyRatioSegments \
        -I ${idSample}.cr.seg \
        -O ${idSample}.called.seg
    """
}

// /* SavvyCNV Related processes */

process SavvyCNVCoverageSummary {
   label 'cpus_16'
    tag "${idSample}"
    cache false
    // publishDir "${params.outdir}/VariantCalling/${idSample}/SavvycnvCoverageSummary", mode: params.publish_dir_mode
    
    input:
         tuple idPatient, idSample, file(bam), file(bai)
    
    output:
    file("${idSample}.coverageBinner")

    when: 'savvycnv' in tools

    script:
    
    """
    java -Xmx1g CoverageBinner ${bam} > ${idSample}.coverageBinner
    """
}

process SavvyCNV {
    label 'cpus_16'
    // tag "${idSample}"
    publishDir "${params.outdir}/VariantCalling/", mode: params.publish_dir_mode
    
    input:
        file("*")
    
    output:
        file("SavvycnvResults")
    
    when: 'savvycnv' in tools
    
    script:
    chunk_size = 200000
    
    """
    mkdir -p SavvycnvResults/SavvycnvCoverageSummary
    mkdir  SavvycnvResults/pdfs
    java -Xmx30g SavvyCNV -a -d ${chunk_size} *.coverageBinner > cnv_list.csv 2>log_messages.txt
    cp *.coverageBinner SavvycnvResults/SavvycnvCoverageSummary/
    cp *.cnvs.pdf SavvycnvResults/pdfs
    cp cnv_list.csv log_messages.txt SavvycnvResults/
    """
}

/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

process MultiQC {
    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode
    input:
        file (multiqcConfig) 
        file (versions) 
        file ('bamQC/*') 
        file ('FastQC/*') 
        file ('BCFToolsStats/*') 
        file ('VCFTools/*')
        // file ('MarkDuplicates/*') 
        file ('SamToolsStats/*') 
        file ('CollectAlignmentSummary/*')
        file ('CollectInsertSizeMetrics/*')
        file ('CollectHsMetrics/*')
        // file ('snpEff/*') 

    output:
        set file("*multiqc_report.html"), file("*multiqc_data") 

    when: !('multiqc' in skipQC)

    script:
    """
    multiqc -f -v .
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
    if (params.bait_bed)           summary['BAIT BED']        = params.bait_bed
    if (step)                       summary['Step']              = step
    if (params.tools)               summary['Tools']             = tools.join(', ')
    if (params.skip_qc)              summary['QC tools skip']     = skipQC.join(', ')

    if (params.no_intervals && step != 'annotate') summary['Intervals']         = 'Do not use'
    if ('haplotypecaller' in tools)                summary['GVCF']              = params.no_gvcf ? 'No' : 'Yes'
    // if ('strelka' in tools && 'manta' in tools )   summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
    if (params.sequencing_center)                  summary['Sequenced by']      = params.sequencing_center
    if (params.pon && 'mutect2' in tools)          summary['Panel of normals']  = params.pon
    // if (params.create_read_count_pon)          summary['Panel of normals']  = params.pon

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
    if (params.fasta_gz)              summary['fasta_gz']              = params.fasta_gz
    if (params.fasta_gz_fai)          summary['fasta_gz_fai']        = params.fasta_fai
    if (params.fasta_gzi)              summary['fasta_gzi']              = params.fasta_gzi

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
    if (params.snpEff_db)              summary['snpEff_db']              = params.snpEff_db
    if (params.species)               summary['species']               = params.species
    if (params.vep_cache_version)       summary['vep_cache_version']       = params.vep_cache_version
    if (params.species)               summary['species']               = params.species
    if (params.snpEff_cache)          summary['snpEff_cache']          = params.snpEff_cache
    if (params.vep_cache)             summary['vep_cache']             = params.vep_cache
    if (params.cadd_InDels)            summary['cadd_InDels']          = params.cadd_InDels
    if (params.cadd_WG_SNVs)            summary['cadd_WG_SNVs']          = params.cadd_WG_SNVs
    if (params.giab_highconf_vcf)            summary['GIAB Truth Set']          = params.giab_highconf_vcf
    if (params.chco_highqual_snps)     summary['Children Colorado hig quality SNPs'] = params.chco_highqual_snps
    if (params.cadd_WG_SNVs_tbi)            summary['cadd_WG_SNVs_tbi']          = params.cadd_WG_SNVs_tbi
    

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
    // params.properties.each { log.info "${it.key} -> ${it.value}" }
    // log.info params.dump()
    log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
    if (params.monochrome_logs) log.info "----------------------------------------------------"
    else log.info "\033[2m----------------------------------------------------\033[0m"
    // log.info("ByeBye")
    // println(params)
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
    // println("passed list: $list")
    // println("real list: $realList")
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
        'alignment_summary',
        'markduplicates',
        'multiqc',
        'samtools',
        'hs_metrics',
        'on_target_assessment',
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
        'markdups',
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
        'deepvariant',
        'benchmark_dv_and_hc_against_giab',
        'benchmark_dv_against_hc',
        // 'hap_py',
        'manta',
        'merge',
        'mpileup',
        'gatkcnv',
        'savvycnv',
        'mutect2',
        'mutect2_single',
        'gen_somatic_pon',
        'gen_read_count_pon',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
// def extractBam(tsvFile) {
//     Channel.from(tsvFile)
//         .splitCsv(sep: '\t')
//         .map { row ->
//             checkNumberOfItem(row, 6)
//             def idPatient = row[0]
//             def gender    = row[1]
//             def status    = returnStatus(row[2].toInteger())
//             def idSample  = row[3]
//             def bamFile   = returnFile(row[4])
//             def baiFile   = returnFile(row[5])

//             if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
//             if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

//             return [idPatient, gender, status, idSample, bamFile, baiFile]
//         }
// }

def extractBam(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 6)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        def baiFile   = returnFile(cols[5])
        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        infos.add([idPatient, gender, status, idSample, bamFile, baiFile])
    }
    return Channel.from(infos)
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

def info(message){
    // log.info(message)
    if (params.debug){
        log.info(message)
    }
}

def extractFastq(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 7)
        

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")
        info("cols[6]:${cols[6]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def idRun      = cols[4]
        def file1      = returnFile(cols[5])
        def file2      = "null"
        file2 = returnFile(cols[6])
        if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"

        infos.add([idPatient, gender, status, idSample, idRun, file1, file2])
    }
    return Channel.from(infos)
}

def extractMarkDups(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 6)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        def baiFile   = returnFile(cols[5])
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
    
        infos.add([idPatient, gender, status, idSample, bamFile, baiFile])
    }
    return Channel.from(infos)
}


def extractRecal(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 7)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")
        info("cols[6]:${cols[6]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        def baiFile   = returnFile(cols[5])
        def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"            

        infos.add([idPatient, gender, status, idSample, bamFile, baiFile, recalTable])
    }
    return Channel.from(infos)
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
// getvcfs to annotate
def getVCFsToAnnotate(results_dir, annotate_tools, tsv){
    vcf_to_annotate = Channel.empty()
    // vcf_no_annotate = Channel.empty()

    if (tsvPath == []) {
    // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
    // This field is used to output final annotated VCFs in the correct directory
        println("annotate_tools: ${annotate_tools}")
        vcf_to_ann = 
        Channel.empty().mix(
        Channel.fromPath("${results_dir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
            .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
            .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Mutect2/*.vcf.gz")
            .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
            .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/TIDDIT/*.vcf.gz")
            .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
        )
        .filter {
            annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools)
        }
    } else if (annotate_tools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, Sarek assume that the idSample should be assumed automatically
      vcf_to_annotate = Channel.fromPath(tsvPath)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    //vcfNoAnnotate.close()
    //vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
    return vcf_to_annotate
}