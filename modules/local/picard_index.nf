process PICARD_INDEX {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/picard:2.9':
        'docker.io/mskcc/picard:2.9' }"

    publishDir "${params.outdir}/${meta.id}/", pattern: "*", mode: params.publish_dir_mode

    input:

    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.abra.bam", includeInputs:true), path("*.bai")      , emit: bam
    path "versions.yml"                                                         , emit: versions

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
        BuildBamIndex \
        I=${bam} \
        O=${bam.baseName}.bai

    cp ${bam.baseName}.bai ${bam.baseName}.bam.bai

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
    touch ${prefix}.bam.bai
    touch ${prefix}.bai


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.9
        r: 3.5.1
        java: 8
    END_VERSIONS
    """
}
