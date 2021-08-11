// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALCULATESPOT {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
    mode: params.publish_dir_mode,
    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_files:'false', publish_dir:'aligned', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::pybedtools=0.8.2--py27h6a42192_1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1"
    } else {
        container "quay.io/biocontainers/pybedtools:0.8.2--py27h6a42192_1"
    }

    input:
    tuple val(meta), path(reads_bam), path(reads_bai), val(bamids), val(names)
    path(index)
    tuple val(imeta), path(interval_beds)
    
    output:
    tuple val(meta), path('*report.txt') , emit: report
    path  "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    calculate_SPoT.py --vers >${software}.version.txt
    
    calculate_SPoT.py \\
        --reads_bam \"${bamids}\" \\
        --intervals all \\
        --name \"${names}\" \\
        --iname all \\
        --g ${index} \\
        --rand \\
        --o ${meta.id}
    """
}
