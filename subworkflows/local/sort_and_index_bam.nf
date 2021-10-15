//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//
def modules = params.modules.clone()
def picard_sortsam_co_options                = modules['picard/sortsam_co']
def samtools_index_options                   = modules['samtools/index']

include { PICARD_SORTSAM } from '../../modules/nf-core/modules/picard/sortsam/main' addParams( options: picard_sortsam_co_options  )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: samtools_index_options  )

workflow SORT_AND_INDEX_BAM {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
    PICARD_SORTSAM     ( ch_bam, 'coordinate' )
    SAMTOOLS_INDEX     ( PICARD_SORTSAM.out.bam )

    emit:
    bam      = PICARD_SORTSAM.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    version  = PICARD_SORTSAM.out.version       //    path: *.version.txt
}
