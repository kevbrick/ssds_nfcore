//
// Sort, mark duplicates, index BAM file and get metrics
//
def modules = params.modules.clone()
def picard_sortsam_co_options                = modules['picard/sortsam_co']
picard_sortsam_co_options                   += [suffix: 'Aligned']

def picard_markduplicates_options            = modules['picard/markduplicates']
def samtools_index_options                   = modules['samtools/index']

include { PICARD_SORTSAM }                from '../../modules/nf-core/modules/picard/sortsam/main'                addParams( options: picard_sortsam_co_options  )
include { PICARD_MARKDUPLICATES }         from '../../modules/nf-core/modules/picard/markduplicates/main'         addParams( options: picard_markduplicates_options  )
include { SAMTOOLS_INDEX }                from '../../modules/nf-core/modules/samtools/index/main'                addParams( options: samtools_index_options  )

workflow MARKDUPLICATES_AND_INDEX {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    fasta
    
    main:
    PICARD_SORTSAM     ( ch_bam, 'coordinate' )

    PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)   

    SAMTOOLS_INDEX (PICARD_MARKDUPLICATES.out.bam)\

    emit:
    bam              = PICARD_MARKDUPLICATES.out.bam              // channel: [ val(meta), [ bam ] ]
    bai              = SAMTOOLS_INDEX.out.bai                     // channel: [ val(meta), [ bai ] ]
    mdmetrics        = PICARD_MARKDUPLICATES.out.metrics          // channel: [ val(meta), [ bam ] ]
    
    version          = PICARD_MARKDUPLICATES.out.version          // path: *.version.txt
    picard_version   = PICARD_MARKDUPLICATES.out.version          // path: *.version.
    samtools_version = SAMTOOLS_INDEX.out.version                 // path: *.version.
}
