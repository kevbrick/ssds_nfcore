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

if (file(params.fasta).isFile()){ 
    ch_genome_index = file("${params.fasta}.fai", checkIfExists: true) 
}else{
    ch_genome_index = file("${params.fasta}.fai", checkIfExists: true)     
}

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

if (params.ssds_interval_beds[ssds_intervals_genome]){
    ch_ssds_intervals = Channel.from(params.ssds_interval_beds[ssds_intervals_genome])
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
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//

def trimgalore_options                       = modules['trimgalore']
def bwa_mem_options                          = modules['bwa/mem']
def bwa_index_options                        = modules['bwa/index']
def bedtools_sort_options                    = modules['bedtools/sort']
def samtools_index_options                   = modules['samtools/index']
def samtools_idxstats_options                = modules['samtools/idxstats']
def samtools_view_options                    = modules['samtools/view']
def parse_ssds_bam_options                   = modules['parse_ssds_bam']
def picard_sortsam_parsed_options            = modules['picard/sortsam_parsed']
def picard_sortsam_qry_options               = modules['picard/sortsam_qry']
def picard_sortsam_co_options                = modules['picard/sortsam_co']
def picard_markduplicates_options            = modules['picard/markduplicates']
def deeptools_plotfingerprint_options        = modules['deeptools/plotfingerprint']
def ucsc_bedgraphtobigwig_options            = modules['ucsc/bedgraphtobigwig']
def phantompeakqualtools_options             = modules['phantompeakqualtools']
def calculatespot_options                    = modules['calculatespot']

def get_ssds_intervals_options               = modules['get_ssds_intervals']
get_ssds_intervals_options                  += [genome: ssds_intervals_genome]

def bedtools_makewindows_options             = modules['bedtools/makewindows']
bedtools_makewindows_options                += [genome: params.genome]

include { FASTQC}                                   from '../modules/nf-core/modules/fastqc/main'                    addParams( options: modules['fastqc'] )
include { MULTIQC}                                  from '../modules/local/multiqc/main'                             addParams( options: multiqc_options )
include { TRIMGALORE}                               from '../modules/nf-core/modules/trimgalore/main'                addParams( options: trimgalore_options )
include { BWA_INDEX }                               from '../modules/nf-core/modules/bwa/index/main'                 addParams( options: bwa_index_options )
include { BWA_MEM }                                 from '../modules/nf-core/modules/bwa/mem/main'                   addParams( options: bwa_mem_options )
include { BEDTOOLS_SORT as SORT_SSDS_BEDS}          from '../modules/nf-core/modules/bedtools/sort/main'             addParams( options: bedtools_sort_options )
include { SAMTOOLS_INDEX }                          from '../modules/nf-core/modules/samtools/index/main'            addParams( options: samtools_index_options )
include { SAMTOOLS_INDEX as INDEX_SSDS_BAMS}        from '../modules/nf-core/modules/samtools/index/main'            addParams( options: samtools_index_options )
include { SAMTOOLS_IDXSTATS }                       from '../modules/nf-core/modules/samtools/idxstats/main'         addParams( options: samtools_idxstats_options )
include { SAMTOOLS_VIEW as FILTER_SUPP_ALIGNMENTS } from '../modules/nf-core/modules/samtools/view/main'             addParams( options: samtools_view_options )
include { PICARD_SORTSAM as SORT_SSDS_BAMS }        from '../modules/nf-core/modules/picard/sortsam/main'            addParams( options: picard_sortsam_parsed_options )
include { PICARD_SORTSAM as PICARD_SORTSAM_QRY}     from '../modules/nf-core/modules/picard/sortsam/main'            addParams( options: picard_sortsam_qry_options )
include { PICARD_SORTSAM as PICARD_SORTSAM_CO }     from '../modules/nf-core/modules/picard/sortsam/main'            addParams( options: picard_sortsam_co_options )
include { PICARD_MARKDUPLICATES               }     from '../modules/nf-core/modules/picard/markduplicates/main'     addParams( options: picard_markduplicates_options )
include { DEEPTOOLS_PLOTFINGERPRINT           }     from '../modules/nf-core/modules/deeptools/plotfingerprint/main' addParams( options: deeptools_plotfingerprint_options )
include { UCSC_BEDGRAPHTOBIGWIG }                   from '../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'     addParams( options: ucsc_bedgraphtobigwig_options )
include { PHANTOMPEAKQUALTOOLS }                    from '../modules/nf-core/modules/phantompeakqualtools/main'      addParams( options: phantompeakqualtools_options )
include { PARSESSDS }                               from '../modules/local/parse_ssds_bam'                           addParams( options: parse_ssds_bam_options )
include { GENERATE_SSDS_COVERAGE }                  from '../modules/local/generate_SSDS_coverage'                   addParams( options: [publish_files:"false"] )
include { BEDTOOLS_MAKEWINDOWS }                    from '../modules/local/bedtools_makewindows'                     addParams( options: bedtools_makewindows_options )
//include { GET_SSDS_INTERVALS }                      from '../modules/local/get_ssds_intervals'                       addParams( options: get_ssds_intervals_options )
include { CALCULATESPOT }                           from '../modules/local/calculatespot'                            addParams( options: calculatespot_options )

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

    // //
    // // MODULE: Run Get ssds intervals
    // //
    // GET_SSDS_INTERVALS (
    //     ch_ssds_intervals, ch_genome_index
    // )
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Run TrimGalore
    //
    TRIMGALORE (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Bwa_mem    
    //
    if (!params.bwa){
        BWA_INDEX (
            file(params.fasta)
        )
        ch_bwa = BWA_INDEX.out.index
    }
    
    BWA_MEM (
        TRIMGALORE.out.reads, ch_bwa
    )
    ch_software_versions = ch_software_versions.mix(BWA_MEM.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Samtools view    
    //
    FILTER_SUPP_ALIGNMENTS (
        BWA_MEM.out.bam
    )
    ch_software_versions = ch_software_versions.mix(FILTER_SUPP_ALIGNMENTS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Picard_sortsam    
    //
    PICARD_SORTSAM_QRY (
        FILTER_SUPP_ALIGNMENTS.out.bam, 'queryname'
    )
    ch_software_versions = ch_software_versions.mix(PICARD_SORTSAM_QRY.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Picard_sortsam    
    //
    PICARD_SORTSAM_CO (
        FILTER_SUPP_ALIGNMENTS.out.bam, 'coordinate'
    )

    //
    // MODULE: Run Picard_markduplicates    
    //
    PICARD_MARKDUPLICATES(
        PICARD_SORTSAM_CO.out.bam    
    )
    ch_software_versions = ch_software_versions.mix(PICARD_MARKDUPLICATES.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Samtools_index    
    //
    SAMTOOLS_INDEX (
        PICARD_MARKDUPLICATES.out.bam
    )
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_INDEX.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Parse_ssds_bam    
    //
    PARSESSDS(
        PICARD_SORTSAM_QRY.out.bam 
    )
    
    ch_ssds_bam_unsorted = PARSESSDS.out.bam.collect{it[1]}.flatten().map {
                                                 def id   = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$1.$2')
                                                 def name = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$1')
                                                 def type = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$2')
                                                 def meta = [id:"${id}", name:"${name}", type:"${type}", "single-end":"false"];
                                                 return [meta, it]}
                                           .groupTuple()
                                          
                                
    ch_ssds_bed_unsorted = PARSESSDS.out.bed.collect{it[1]}.flatten().map {
                                                 def id   = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$1.$2')
                                                 def name = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$1')
                                                 def type = it.name.replaceFirst('^(.+)_(ssDNA|ssDNA_type2|dsDNA_loconf|dsDNA_hiconf|unclassified)\\..+$','$2')
                                                 def meta = [id:"${id}", name:"${name}", type:"${type}", "single-end":"false"];
                                                 return [meta, it]}
                                           .groupTuple()

    ch_software_versions = ch_software_versions.mix(PARSESSDS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Bedtools_sort
    //
    SORT_SSDS_BEDS(
        ch_ssds_bed_unsorted
    )
    ch_software_versions = ch_software_versions.mix(SORT_SSDS_BEDS.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Generate ssds coverage
    //    
    BEDTOOLS_MAKEWINDOWS (
        [[id:"${params.genome}"], ch_genome_index] 
    )
    ch_software_versions = ch_software_versions.mix(BEDTOOLS_MAKEWINDOWS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Generate ssds coverage
    //
    def ch_ssds_beds = SORT_SSDS_BEDS.out.bed
                        .filter {file(it[1]).size() > 40000}

    GENERATE_SSDS_COVERAGE(
        ch_ssds_beds, ch_genome_index, BEDTOOLS_MAKEWINDOWS.out.bed
    )
    ch_software_versions = ch_software_versions.mix(SORT_SSDS_BEDS.out.version.first().ifEmpty(null))

    ch_coverage_bg = GENERATE_SSDS_COVERAGE.out.bedgraph.flatten().map {
                                                    def id   = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1.$2')
                                                    def name = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1')
                                                    def meta=[id:"${id}", name:"${name}", "single-end":"false"];
                                                    return [meta, it]}
                                           .groupTuple()
    
    //
    // MODULE: Run Ucsc bedgraphtobigwig  
    //    
    UCSC_BEDGRAPHTOBIGWIG (
        ch_coverage_bg, ch_genome_index
    )
    ch_software_versions = ch_software_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Picard_sortsam    
    //
    SORT_SSDS_BAMS (
        ch_ssds_bam_unsorted, 'coordinate'
    )
    ch_software_versions = ch_software_versions.mix(SORT_SSDS_BAMS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Samtools index    
    //
    INDEX_SSDS_BAMS (
        SORT_SSDS_BAMS.out.bam
    )
    ch_software_versions = ch_software_versions.mix(INDEX_SSDS_BAMS.out.version.first().ifEmpty(null))
    
    ch_ssds_bam = SORT_SSDS_BAMS.out.bam.mix(INDEX_SSDS_BAMS.out.bai).groupTuple(by:0)
        .map { [it[0], it[1][0], it[1][1]] }
    
    //
    // MODULE: Run Samtools idxstats    
    //
    SAMTOOLS_IDXSTATS (
        ch_ssds_bam
    )
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_IDXSTATS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Phantompeakqualtools  
    //    
    PHANTOMPEAKQUALTOOLS (
        SORT_SSDS_BAMS.out.bam
    )
    ch_software_versions = ch_software_versions.mix(PHANTOMPEAKQUALTOOLS.out.version.first().ifEmpty(null))
    
    //
    // MODULE: Run Calculatespot  
    // 
    ch_for_spot = ch_ssds_bam.map { [["id":it[0].name, "single-end":it[0]."single-end"], it[1], it[2], it[1].name, it[0].type ] }
                              .groupTuple(by:0)
        
    CALCULATESPOT (
        ch_for_spot, ch_genome_index, ch_ssds_intervals.collect()
    )
    ch_software_versions = ch_software_versions.mix(CALCULATESPOT.out.version.first().ifEmpty(null))
    
    ch_spot_report = CALCULATESPOT.out.report.map {[it[0].name, it[1]]}
                         .collectFile(storeDir:"${params.outdir}/reports") {
                             [ "${it[0]}.SSDS_SPoT_report.txt", it[1] ] }

    //
    // MODULE: Run Deeptools plotfingerprint
    //
    
    ch_ssds_bam.map {[it[1], file(it[1]).size(), file(params.fasta).size()]}
    
    // PLOTFINGERPRINT FAILS IF TOO FEW GENOMIC WINDOWS CONTAIN READS. THUS, WE ONLY RUN THIS 
    // FOR BAMS THAT ARE AT LEAST 75% THE SIZE OF THE GENOME FASTA (A BIT ARBITRARY)
    DEEPTOOLS_PLOTFINGERPRINT (
        ch_ssds_bam.filter {file(it[1]).size() > file(params.fasta).size()*0.75}
    )
    ch_software_versions = ch_software_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.version.first().ifEmpty(null))
    
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
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PARSESSDS.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_spot_report.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.pdf.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.rdata.collect{it[1]}.ifEmpty([]))
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
