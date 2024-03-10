process PICARD_COLLECT_MULTIPLE_METRICS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/picard:2.9':
        'docker.io/mskcc/picard:2.9' }"

    publishDir "${params.outdir}/${meta.id}/", pattern: "*", mode: params.publish_dir_mode

    input:

    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.quality_by_cycle_metrics")              , emit: qual_metrics
    tuple val(meta), path("*.quality_by_cycle.pdf")                  , emit: qual_pdf
    path "versions.yml"                                              , emit: versions

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
        /usr/bin/picard-tools/picard.jar \
        CollectMultipleMetrics \
        I=${bam} \
        REFERENCE_SEQUENCE=${fasta} \
        PROGRAM=null \
        PROGRAM=MeanQualityByCycle \
        VALIDATION_STRINGENCY=SILENT \
        OUTPUT=${bam.baseName}.qmetrics

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.9
        r: 3.5.1
        java: 8
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qmetrics.quality_by_cycle_metrics
    touch ${prefix}.qmetrics.uality_by_cycle.pdf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.9
        r: 3.5.1
        java: 8
    END_VERSIONS
    """
}
