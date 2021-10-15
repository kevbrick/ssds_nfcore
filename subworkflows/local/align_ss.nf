//
// SSDS alignment 
//
include { MARKDUPLICATES_AND_INDEX as MD_AND_INDEX}     from '../../subworkflows/local/markduplicates_and_index_bam' addParams( options: [:] )
include { MARKDUPLICATES_AND_INDEX as MD_AND_INDEX_RAW} from '../../subworkflows/local/markduplicates_and_index_bam' addParams( options: [:] )

def modules = params.modules.clone()

def trimgalore_options                       = modules['trimgalore']
def bwa_mem_options                          = modules['bwa/mem']
def samtools_view_remove_supp_options        = modules['samtools/view']
samtools_view_remove_supp_options           += [args: "-F 2048 -hb"]
samtools_view_remove_supp_options           += [suffix: '_noSupp']

def picard_collectmultiplemetrics_options    = modules['picard/collectmultiplemetrics']
def parse_ssds_bam_options                   = modules['parse_ssds_bam']
def picard_sortsam_qry_options               = modules['picard/sortsam_qry']
def bedtools_sort_options                    = modules['bedtools/sort']
def samtools_idxstats_options                = modules['samtools/idxstats']

include { TRIMGALORE }                              from '../../modules/nf-core/modules/trimgalore/main'                    addParams( options: trimgalore_options  )
include { BWA_MEM }                                 from '../../modules/nf-core/modules/bwa/mem/main'                       addParams( options: bwa_mem_options  )
include { SAMTOOLS_VIEW as FILTER_SUPP_ALIGNMENTS}  from '../../modules/nf-core/modules/samtools/view/main'                 addParams( options: samtools_view_remove_supp_options  )
include { PICARD_COLLECTMULTIPLEMETRICS }           from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main' addParams( options: picard_collectmultiplemetrics_options  )
include { PARSESSDS }                               from '../../modules/local/parse_ssds_bam'                               addParams( options: parse_ssds_bam_options )
include { PICARD_SORTSAM as PICARD_SORTSAM_QRY}     from '../../modules/nf-core/modules/picard/sortsam/main'                addParams( options: picard_sortsam_qry_options )
include { BEDTOOLS_SORT}                            from '../../modules/nf-core/modules/bedtools/sort/main'                 addParams( options: bedtools_sort_options )
include { SAMTOOLS_IDXSTATS }                       from '../../modules/nf-core/modules/samtools/idxstats/main'             addParams( options: samtools_idxstats_options )

workflow ALIGN_SS {
    take:
    fastq // channel: [ val(meta), [ bam ] ]
    bwa
    fasta
    
    main:
    
    TRIMGALORE (fastq)
    BWA_MEM (TRIMGALORE.out.reads, bwa)
    FILTER_SUPP_ALIGNMENTS (BWA_MEM.out.bam)
    
    ch_unprocessed = FILTER_SUPP_ALIGNMENTS.out.bam.map {return [[id:"${it[0].id}.unprocessed_SSDS", single_end: it[0].single_end], it[1]]}
    MD_AND_INDEX_RAW(ch_unprocessed, fasta)
    
    PICARD_COLLECTMULTIPLEMETRICS(MD_AND_INDEX_RAW.out.bam, fasta)
    
    PICARD_SORTSAM_QRY (FILTER_SUPP_ALIGNMENTS.out.bam.map {return [[id:it[0].id.replaceAll(".unprocessed_SSDS",""), single_end: it[0].single_end], it[1]]}, 'queryname')
    PARSESSDS(PICARD_SORTSAM_QRY.out.bam)
                                           
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

    BEDTOOLS_SORT(ch_ssds_bed_unsorted)
    
    MD_AND_INDEX (ch_ssds_bam_unsorted, fasta)
    
    ch_ssds_bam = MD_AND_INDEX.out.bam.mix(MD_AND_INDEX.out.bai).groupTuple(by:0)
        .map { [it[0], it[1][0], it[1][1]] }
    
    
    SAMTOOLS_IDXSTATS (ch_ssds_bam)
    
    ch_reports = TRIMGALORE.out.log.collect()
                 .mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics)
                 .mix(MD_AND_INDEX_RAW.out.mdmetrics)
                 .mix(PARSESSDS.out.report)
    
    emit:
    bam     = ch_ssds_bam             // channel: [ val(meta), [ bam ] ]
    bed     = BEDTOOLS_SORT.out.bed   // channel: [ val(meta), [ bai ] ]
    reports = ch_reports
    
    version           = TRIMGALORE.out.version           //    path: *.version.txt
    bwa_version       = BWA_MEM.out.version              //    path: *.version.txt
    samtools_version  = SAMTOOLS_IDXSTATS.out.version    //    path: *.version.txt
    picard_version    = MD_AND_INDEX.out.picard_version  //    path: *.version.txt
    bedtools_version  = BEDTOOLS_SORT.out.version        //    path: *.version.txt
    parsessds_version = PARSESSDS.out.version            //    path: *.version.txt
}
