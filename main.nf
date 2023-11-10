#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowPG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowPG
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { extractFastqPairFromDir } from './NextflowModules/Utils/fastq.nf'

include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'
include { SAMBAMBA_MARKDUP } from './modules/nf-core/sambamba/markdup/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Reference file channels
    ch_genome = Channel.fromPath(params.genome).map {genome -> [genome.getSimpleName(), genome] }

    // Input channel
    ch_fastq = extractFastqPairFromDir(params.input)

    // Mapping
    BWAMEM2_MEM(ch_fastq, ch_genome, true)
    SAMBAMBA_MARKDUP(BWAMEM2_MEM.out.bam.map{ meta, bam -> [ meta - meta.subMap('read_group'), bam ] }.groupTuple())
    SAMTOOLS_INDEX(SAMBAMBA_MARKDUP.out.bam)
    ch_bam_bai = SAMBAMBA_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.bai)

    // Variant calling
    // GATK4_HAPLOTYPECALLER(SAMBAMBA_MARKDUP.out.bam)
    // GATK_GenotypeGVCFs

    // GLIMS output
    // VCF2GLIMS


    // QC
    FASTQC(ch_fastq)
    // Mosdepth
    // VerifyBamID2

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )

}
