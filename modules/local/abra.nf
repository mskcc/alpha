process ABRA {


    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://mskcc/abra:2.17':
        'docker.io/mskcc/abra:2.17' }"

    publishDir "${params.outdir}/${meta.id}/", pattern: "*", mode: params.publish_dir_mode

    input:
    tuple val(meta),  path(tumor), path(tumor_index), path(normal), path(normal_index), path(targets)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("tumor/*.abra.bam"), path("normal/*.abra.bam")     , emit: bams
    path "versions.yml"                                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus * 2

    """
    mkdir ./tmp
    java \
        -Xms${task.memory.toMega()/4}m \
        -Xmx${task.memory.toGiga()}g \
        -jar \
        /usr/bin/abra.jar \
        --tmpdir \
        ./tmp \
        --threads ${task.cpus * 2} \
        --ref ${fasta} \
        --targets ${targets} \
        --out ${tumor.baseName}.abra.bam,${normal.baseName}.abra.bam \
        --in ${tumor},${normal}

    mkdir tumor
    mkdir normal
    cp ${tumor.baseName}.abra.bam tumor
    cp ${normal.baseName}.abra.bam normal


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra: \$(java -jar /usr/bin/abra.jar 2>&1 | head -1 | sed 's/.*Abra version: //')
        java: \$(java -version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir tumor
    mkdir normal
    touch tumor/${tumor.baseName}.abra.bam
    touch normal/${normal.baseName}.abra.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra: \$(java -jar /usr/bin/abra.jar 2>&1 | head -1 | sed 's/.*Abra version: //')
        java: \$(java -version 2>&1 | head -1)
    END_VERSIONS
    """

}
