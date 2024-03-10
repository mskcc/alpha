process GATK_BQSR {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/gatk:3.3-0':
        'docker.io/mskcc/gatk:3.3-0' }"

    publishDir "${params.outdir}/${meta.id}/", pattern: "*", mode: params.publish_dir_mode

    input:

    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(known_sites), path(known_sites_index)

    output:
    tuple val(meta), path("*.recal.matrix")               , emit: recal_matrix
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known_sites_list = known_sites.join(" --knownSites ")

    """
    java \
        -Xms${task.memory.toMega()/4}m \
        -Xmx${task.memory.toGiga()}g \
        -XX:-UseGCOverheadLimit \
        -Djava.io.tmpdir=./tmp \
        -jar \
        /usr/bin/gatk.jar \
        -T \
        BaseRecalibrator \
        ${args} \
        --input_file \
        ${tumor_bam} \
        --input_file \
        ${normal_bam} \
        --reference_sequence \
        ${fasta} \
        --num_cpu_threads_per_data_thread \
        ${task.cpus * 2} \
        --knownSites \
        ${known_sites_list} \
        --out \
        ${prefix}.recal.matrix

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
    touch ${prefix}.recal.matrix


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: 3.3-0
        java: 8
    END_VERSIONS
    """
}
