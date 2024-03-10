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
        REALIGN_PAIR(sample_name,tumor_bam,normal_bam,genome,genome_fasta,roi_bed)

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

