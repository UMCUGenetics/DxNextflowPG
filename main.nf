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
include { extractFastqPairFromDir } from './modules/local/utils/fastq.nf'

include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS } from './modules/nf-core/gatk4/genotypegvcfs/main'
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
    ch_bwa_index = Channel.fromPath("${params.genome}*").map {genome -> [genome.getSimpleName(), genome] }.groupTuple().collect()
    ch_genome_fasta = Channel.fromPath("${params.genome}").collect()
    ch_genome_fasta_index = Channel.fromPath("${params.genome}.fai").collect()
    ch_genome_dict = Channel.fromPath("${params.genome_dict}").collect()
    ch_dbsnp = Channel.fromPath("${params.dbsnp}").collect()
    ch_dbsnp_index = Channel.fromPath("${params.dbsnp}.idx").collect()

    ch_intervals = Channel.fromPath("${params.intervals}").collect()

    // Input channel
    ch_fastq = extractFastqPairFromDir(params.input)

    // // Mapping
    BWAMEM2_MEM(ch_fastq, ch_bwa_index, true)
    SAMBAMBA_MARKDUP(BWAMEM2_MEM.out.bam.map{ meta, bam -> [ meta - meta.subMap('rg_id', 'flowcell'), bam ] }.groupTuple())
    SAMTOOLS_INDEX(SAMBAMBA_MARKDUP.out.bam)

    ch_bam_bai = SAMBAMBA_MARKDUP.out.bam.join(SAMTOOLS_INDEX.out.bai)

    // Variant calling
    GATK4_HAPLOTYPECALLER(
        ch_bam_bai.combine(ch_intervals).map{ meta, bam, bai, intervals -> [meta, bam, bai, intervals, [] ] },
        ch_genome_fasta, ch_genome_fasta_index, ch_genome_dict, ch_dbsnp, ch_dbsnp_index
    )
    GATK4_HAPLOTYPECALLER.out.vcf.view()
    GATK4_GENOTYPEGVCFS(
        GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.index).combine(ch_intervals).map{
            meta, vcf, tbi , intervals -> [meta, vcf, tbi, intervals, [] ]
        },
        ch_genome_fasta, ch_genome_fasta_index, ch_genome_dict, ch_dbsnp, ch_dbsnp_index
    )

    // GLIMS output
    // VCF2GLIMS


    // QC
    // FASTQC(ch_fastq)
    // // Mosdepth
    // // VerifyBamID2

    // // MultiQC
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // MULTIQC(
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     Channel.empty().toList(),
    //     Channel.empty().toList()
    // )

}
