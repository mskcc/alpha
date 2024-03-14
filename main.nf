#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alpha
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mskcc/alpha
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { REALIGN_PAIR } from './subworkflows/local/realign_pair.nf'
include { GATK_BQSR } from './modules/local/gatk_bqsr.nf'

workflow RECAL_PAIR {

// tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index)
// tuple val(meta2), path(fasta), path(fai)
// tuple val(meta3), path(known_sites), path(known_sites_index)


    take:
        sample_name
        tumor_ra_bam
        normal_ra_bam
        genome
        genome_fasta
        genome_fasta_index
        ch_known_sites

    main:

        ch_indexed_bams = normal_ra_bam.join(tumor_ra_bam)

        println "====================="
        println "tumor_bam"
        ch_indexed_bams.view()

        GATK_BQSR(
            ch_indexed_bams,
            tuple(
                genome,
                genome_fasta,genome_fasta_index
                ),
            ch_known_sites
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
        ch_known_sites


    main:

        ch_version = Channel.empty()

        tumor_bam_index=getIndexFromPath(tumor_bam)
        normal_bam_index=getIndexFromPath(normal_bam)
        genome_fasta_index=getIndexFromPath(genome_fasta)

        results = REALIGN_PAIR(
                        sample_name,
                        tumor_bam,tumor_bam_index,
                        normal_bam,normal_bam_index,
                        genome,
                        genome_fasta,genome_fasta_index,
                        roi_bed
                        )

        RECAL_PAIR(
            sample_name,
            results.tumor_bam,
            results.normal_bam,
            genome,
            genome_fasta,genome_fasta_index,
            ch_known_sites
        )

}

workflow {

    known_sites=[]
    for(single_site in params.known_sites){
        kn=file(single_site)
        if(kn.exists()) {
            known_sites.add(kn)
        }
    }

    known_sites_index=[]
    for(single_site in known_sites){
        known_sites_index.add(getIndexFromPath(single_site))
    }

    ch_known_sites = Channel.value(["known_sites", known_sites, known_sites_index])

    ch_known_sites.view()

    BRAVO(
        params.sample_name,
        file(params.tumor_bam),
        file(params.normal_bam),
        params.genome,
        file(params.fasta),
        file(params.roi_bed),
        ch_known_sites
        )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Util functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getBamIndexPath(bam) {
    bam_index=file(bam.toString().replace(".bam",".bai"))
    if(!bam_index.exists()) {
        bam_index=bam+".bai"
        if(!bam_index) {
            return(null)
        }
    }
    return(file(bam_index))
}

def getFastaIndexPath(fasta) {
    file(fasta+".fai")
}

def getVcfIndexPath(vcf) {
    file(vcf+".idx")
}

def getIndexFromPath(path) {
    filename=path.toString()
    if(filename.endsWith(".vcf")) {
        return(getVcfIndexPath(path))
    } else if(filename.endsWith(".bam")) {
        return(getBamIndexPath(path))
    } else if(filename.endsWith(".fasta")) {
        return(getFastaIndexPath(path))
    } else {
        return("NA")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
