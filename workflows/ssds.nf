#!/usr/bin/env nextflow
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSsds.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check genome for ssds intervals
def ssds_intervals_genome = params.ssds_intervals_genome
if (!ssds_intervals_genome) {
    ssds_intervals_genome = params.genome ? params.genome : null
}

// Check mandatory parameters;
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (params.bwa)   {
    if (file(params.bwa).isFile()){ 
        ch_bwa = file(params.bwa).getParent() 
    } else {
        ch_bwa = file(params.bwa) 
    }
} 

ch_genome_index  = file("${params.fasta}.fai", checkIfExists: true) 
ch_makewin_input = [ [ id:"${params.genome}"],
                     file("${params.fasta}.fai", checkIfExists: true) ]

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    SSDS INTERVAL FILES
========================================================================================
*/
def ch_ssds_intervals

if (params.ssds_interval_beds["${ssds_intervals_genome}"]){
    ch_ssds_intervals = Channel.from(params.ssds_interval_beds["${ssds_intervals_genome}"])
                                   .map {[[id:it.id, 
                                           genome:ssds_intervals_genome], 
                                           file(it.bed, checkIfExists: true)]}
}else{
    System.err.println("WARNING: No SSDS intervals found for SPoT analysis. SPoT analyses will not be performed.")
    ch_ssds_intervals = Channel.empty()
}
                                    
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK }           from '../subworkflows/local/input_check'                  addParams( options: [:] )
include { SORT_AND_INDEX_BAM }    from '../subworkflows/local/sort_and_index_bam'           addParams( options: [:] )
include { MAKE_SSDS_BIGWIGS }     from '../subworkflows/local/make_ssds_bigwigs'            addParams( options: [:] )
include { ALIGN_SS }              from '../subworkflows/local/align_ss'                     addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def fastqc_options                           = modules['fastqc']
def bwa_index_options                        = modules['bwa/index']
def deeptools_plotfingerprint_options        = modules['deeptools/plotfingerprint']
def bedtools_makewindows_options             = modules['bedtools/makewindows']
bedtools_makewindows_options                += [genome: params.genome]

def calculatespot_options                    = modules['calculatespot']

def get_ssds_intervals_options               = modules['get_ssds_intervals']
get_ssds_intervals_options                  += [genome: ssds_intervals_genome]

def multiqc_options                          = modules['multiqc']
multiqc_options.args                        += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

include { FASTQC}                    from '../modules/nf-core/modules/fastqc/main'                        addParams( options: fastqc_options )
include { BWA_INDEX }                from '../modules/nf-core/modules/bwa/index/main'                     addParams( options: bwa_index_options )
include { DEEPTOOLS_PLOTFINGERPRINT} from '../modules/nf-core/modules/deeptools/plotfingerprint/main'     addParams( options: deeptools_plotfingerprint_options )
include { BEDTOOLS_MAKEWINDOWS }     from '../modules/nf-core/modules/bedtools/makewindows/main'          addParams( options: bedtools_makewindows_options )
include { CALCULATESPOT }            from '../modules/local/calculatespot'                                addParams( options: calculatespot_options )
include { MULTIQC}                   from '../modules/local/multiqc/main'                                 addParams( options: multiqc_options )
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SSDS {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))

    //
    // MODULE/SUBWORKFLOW: Index genome if required, then run SSDS alignment sub-workflow    
    //
    if (!params.bwa){
        BWA_INDEX (
            file(params.fasta)
        )
        ch_bwa = BWA_INDEX.out.index
        ch_software_versions = ch_software_versions.mix(BWA_INDEX.out.version.first().ifEmpty(null))
    }
    
    ALIGN_SS(INPUT_CHECK.out.reads, ch_bwa, file(params.fasta))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.bwa_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.samtools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.picard_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.bedtools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(ALIGN_SS.out.parsessds_version.first().ifEmpty(null))
    
    //
    // MODULE/SUBWORKFLOW: Run MAKE_SSDS_BIGWIGS subworkflow to generate SSDS coverage bigwigs; make genomic windows first
    //
    BEDTOOLS_MAKEWINDOWS (
        ch_makewin_input, false 
    )
    ch_software_versions = ch_software_versions.mix(BEDTOOLS_MAKEWINDOWS.out.version.first().ifEmpty(null))
    
    def ch_ssds_beds = ALIGN_SS.out.bed
                        .filter {file(it[1]).size() > 40000}

    MAKE_SSDS_BIGWIGS(ch_ssds_beds, ch_genome_index, BEDTOOLS_MAKEWINDOWS.out.tab)
    ch_software_versions = ch_software_versions.mix(MAKE_SSDS_BIGWIGS.out.version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(MAKE_SSDS_BIGWIGS.out.version2.first().ifEmpty(null))
    
    //
    // MODULE: Run Calculatespot  
    // 
    ch_bams_for_spot = ALIGN_SS.out.bam.map { [["id":it[0].name, "single-end":it[0]."single-end"], it[1], it[2], it[1].name, it[0].type ] }
                                  .groupTuple(by:0)
    
    ch_intervals_for_spot = ch_ssds_intervals.map{[[genome:it[0].genome],it[1]]}.groupTuple(by:0).collect()

    CALCULATESPOT (
        ch_bams_for_spot, ch_genome_index, ch_intervals_for_spot
    )
    ch_software_versions = ch_software_versions.mix(CALCULATESPOT.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Deeptools plotfingerprint
    //
    ch_fingerprint_bams = ALIGN_SS.out.bam.map { [["id":it[0].name, "single-end":it[0]."single-end"], it[1], it[2]] }
                                     .groupTuple(by:0)

    DEEPTOOLS_PLOTFINGERPRINT (
        ch_fingerprint_bams
    )
    ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first().ifEmpty(null))
    
    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSsds.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_SS.out.reports.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CALCULATESPOT.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect{it[1]}.ifEmpty([]))
    
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
