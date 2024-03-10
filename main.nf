#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alpha
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mskcc/alpha
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { ABRA } from './modules/local/abra'
include { PICARD_INDEX as PICARD_INDEX_tumor } from './modules/local/picard_index'
include { PICARD_INDEX as PICARD_INDEX_normal } from './modules/local/picard_index'

def getBamIndexPath(bam) {
    bam_index=file(bam.toString().replace(".bam",".bai"))
    if(!bam_index.exists()) {
        bam_index=bam+".bai"
        if(!bam_index) {
            return(null)
        }
    }
    return(bam_index)
}

workflow REALIGN_PAIR {

    take:
        sample_name
        tumor_bam
        normal_bam
        genome
        genome_fasta
        roi_bed

    main:

        ch_versions = Channel.empty()

        genome_index=genome_fasta+".fai"
        tumor_index=getBamIndexPath(tumor_bam)
        normal_index=getBamIndexPath(normal_bam)
        meta=["id":sample_name]

        ABRA(
            tuple(
                meta,
                tumor_bam,tumor_index,
                normal_bam,normal_index,
                roi_bed
                ),
            tuple(
                genome,
                genome_fasta,genome_index
                )
            )
        ch_versions=ch_versions.mix(ABRA.out.versions)

        tumor_bam=ABRA.out.bams.map { new Tuple(it[0],it[1]) }
        normal_bam=ABRA.out.bams.map { new Tuple(it[0],it[2]) }

        PICARD_INDEX_tumor(tumor_bam)
        PICARD_INDEX_normal(normal_bam)

        ch_versions=ch_versions.mix(PICARD_INDEX_tumor.out.versions)
        ch_versions=ch_versions.mix(PICARD_INDEX_normal.out.versions)

    emit:
        tumor_bam = PICARD_INDEX_tumor.out.bam
        normal_bam = PICARD_INDEX_normal.out.bam
        versions = ch_versions

}

workflow BRAVO {

    take:
        sample_name
        tumor_bam
        normal_bam
        genome
        genome_fasta
        roi_bed

    main:
        (tumor_ra_bam, normal_ra_bam, version) = REALIGN_PAIR(sample_name,tumor_bam,normal_bam,genome,genome_fasta,roi_bed)

}

workflow {

    BRAVO(
        params.sample_name,
        file(params.tumor_bam),
        file(params.normal_bam),
        params.genome,
        file(params.fasta),
        file(params.roi_bed)
        )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

