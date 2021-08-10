// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENERATE_SSDS_COVERAGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pybedtools=0.8.2--py27h6a42192_1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1"
    } else {
        container "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1"
    }
    
    input:
    tuple val(meta), path(bed)
    path(genome_index)
    path(windows_bed)

    // when:
    // bed.size() > 40000
        
    output:
    path("*.bedgraph")    , emit: bedgraph
    path  "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    grep -w \\+ ${bed} >${prefix}.FWD.bed
    grep -w \\- ${bed} >${prefix}.REV.bed
    
    nF=`cat ${prefix}.FWD.bed |wc -l`
    nR=`cat ${prefix}.REV.bed |wc -l`
    
    make_SSDS_bedgraphs.py --fwd ${prefix}.FWD.bed \
                           --rev ${prefix}.REV.bed \
                           --g ${genome_index} \
                           --name ${prefix} \
                           --win ${windows_bed}
    
    echo "0.8.2" > ${software}.version.txt
    """
}
