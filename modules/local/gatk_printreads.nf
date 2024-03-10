process GATK_PRINTREADS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/gatk:3.3-0':
        'docker.io/mskcc/gatk:3.3-0' }"

    publishDir "${params.outdir}/${meta.id}/", pattern: "*", mode: params.publish_dir_mode

    input:

    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bqsr)

    output:
    tuple val(meta), path("*.printreads.bam"), path("*.bai")      , emit: bam
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    java \
        -Xms${task.memory.toMega()/4}m \
        -Xmx${task.memory.toGiga()}g \
        -XX:-UseGCOverheadLimit \
        -Djava.io.tmpdir=./tmp \
        -jar \
        /usr/bin/gatk.jar \
        -T \
        PrintReads \
        --input_file \
        ${bam} \
        --reference_sequence \
        ${fasta} \
        --num_cpu_threads_per_data_thread \
        ${task.cpus * 2} \
        --emit_original_quals \
        --BQSR \
        ${bqsr} \
        --baq \
        RECALCULATE \
        --out \
        ${bam.baseName}.printreads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: 3.3-0
        java: 8
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.abra.printreads.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: 3.3-0
        java: 8
    END_VERSIONS
    """
}
