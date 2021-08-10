// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SSDS_INTERVALS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gnu-wget=1.18--h5bf99c6_5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    } else {
        container "quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5"
    }
    
    input:
    tuple val(meta), file(intervals)
    path (index)
        
    output:
    tuple val(meta), path("*.intervals.bed"), emit: bed
    path  "*.version.txt",                    emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo ${prefix}
    
    gunzip -c ${intervals} |cut -f1-3 |sort -k1,1 -k2n,2n >${prefix}.intervals.bed
    
    echo "0.1.0" > ${software}.version.txt
    """
}
