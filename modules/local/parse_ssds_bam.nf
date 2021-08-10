// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PARSESSDS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_files:'false', publish_dir:'aligned', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-480c331443a1d7f4cb82aa41315ac8ea4c9c0b45:3e0fc1ebdf2007459f18c33c65d38d2b031b0052-0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bam')        , emit: bam
    tuple val(meta), path('*.bed')        , emit: bed
    tuple val(meta), path('*_report.txt') , emit: report
    path  "*.version.txt"                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    parse_SSDS_BAM.py --vers >${software}.version.txt
    
    parse_SSDS_BAM.py \\
        --bam $bam \\
        --name $prefix
    """
}
