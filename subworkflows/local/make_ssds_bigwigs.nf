//
// Make bigwigs from SSDS fragment BED files
//
def modules = params.modules.clone()
def ucsc_bedgraphtobigwig_options                = modules['ucsc/bedgraphtobigwig']

include { UCSC_BEDGRAPHTOBIGWIG }                   from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'         addParams( options: ucsc_bedgraphtobigwig_options )
include { GENERATE_SSDS_COVERAGE }                  from '../../modules/local/generate_SSDS_coverage'                       addParams( options: [publish_files:"false"] )

workflow MAKE_SSDS_BIGWIGS {
    take:
    ssds_beds // channel: [ val(meta), [ bam ] ]
    genome_index
    genome_windows
    
    main:
    GENERATE_SSDS_COVERAGE(ssds_beds, genome_index, genome_windows)

    ch_coverage_bg = GENERATE_SSDS_COVERAGE.out.bedgraph.flatten().map {
                                                    def id   = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1.$2')
                                                    def name = it.name.replaceFirst('^(.+)\\.(FWD|REV|FR|TOT)\\..+$','$1')
                                                    def meta=[id:"${id}", name:"${name}", "single-end":"false"];
                                                    return [meta, it]}
                                           .groupTuple()
    
    UCSC_BEDGRAPHTOBIGWIG (ch_coverage_bg, genome_index)

    emit:
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig    // channel: [ val(meta), [ bai ] ]

    version   = GENERATE_SSDS_COVERAGE.out.version //    path: *.version.txt
    version2  = UCSC_BEDGRAPHTOBIGWIG.out.versions  //    path: *.version.txt
}
