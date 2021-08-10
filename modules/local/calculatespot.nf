// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALCULATESPOT {
    tag "${reads_meta.id}:${intervals_meta.id}"
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
    tuple val(reads_meta), path(reads_bam), val(intervals_meta), path(intervals_bed)
    path(index)
    
    output:
    tuple val(reads_meta), path('*report.txt') , emit: report
    path  "*.version.txt"                      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def metaid   = "${reads_meta.type}_${intervals_meta.id}"
    def prefix   = options.suffix ? "${metaid}${options.suffix}" : "${metaid}"
    def iname    = intervals_meta.id.replaceAll("[^A-Za-z0-9]+","_")
    """
    calculate_SPoT.py --vers >${software}.version.txt
    
    calculate_SPoT.py \\
        --reads_bam ${reads_bam} \\
        --intervals ${intervals_bed} \\
        --name ${reads_meta.type} \\
        --iname ${iname} \\
        --g ${index} \\
        --rand \\
        --o ${reads_meta.id}.${iname}
    """
}
