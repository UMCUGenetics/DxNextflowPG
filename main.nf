#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowPG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowPG
----------------------------------------------------------------------------------------
*/

//nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate and log parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
log.info paramsSummaryLog(workflow)
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER as FILTER_PYPGX_VCF } from './modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as FILTER_GATK_VCF } from './modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_VIEW } from './modules/nf-core/bcftools/view/main'
include { CALL_STARALLELES } from './modules/local/star_alleles/main'
include { COMBINERESULTS } from './modules/local/combine_outputs/main'
include { COV_QA } from './modules/local/QA/COV_QA'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS } from './modules/nf-core/gatk4/genotypegvcfs/main'
include { MOSDEPTH } from './modules/nf-core/mosdepth/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { PARSE_MOSDEPTH } from './modules/local/utils/parse_mosdepth'
include { PYPGX_CREATEREGIONS } from './modules/local/pypgx/create_regions/main'
include { PYPGX_CREATEINPUTVCF } from './modules/local/pypgx/create_input_vcf/main'
include { PYPGX_PREPAREDEPTHOFCOVERAGE } from './modules/local/pypgx/prepare_depth_of_coverage/main'
include { PYPGX_COMPUTECONTROLSTATISTICS } from './modules/local/pypgx/compute_control_statistics/main'
include { PYPGX_RUNNGSPIPELINE } from './modules/local/pypgx/run_ngs_pipeline/main'
include { SV_QA } from './modules/local/QA/SV_QA'
include { VCF2GLIMS } from './modules/local/vcf2glims/main'
include { VERIFYBAMID_VERIFYBAMID2 } from './modules/nf-core/verifybamid/verifybamid2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create reference file channels, add meta values
    ch_genome_fasta = Channel.fromPath("${params.genome_fasta}")
        .map{ file -> [file.getSimpleName(), file] }
        .collect()

    ch_genome_fasta_index = Channel.fromPath("${params.genome_fasta}.fai")
        .map{ file -> [file.getSimpleName(), file] }
        .collect()

    ch_genome_dict = Channel.fromPath("${params.genome_dict}")
        .map{ file -> [file.getSimpleName(), file] }
        .collect()

    ch_bams_meta = Channel.fromFilePairs(
        "${params.bam_path}/*.{bam,bai}",
        checkIfExists: true) {
            file -> file.name.replaceAll(/.bam|.bai$/,'') }
        .map{ meta, bam_index -> [['id': meta], bam_index[0], bam_index[1]] }

    ch_dbsnp = Channel.fromPath("${params.dbsnp}")
        .map{ file -> [file.getSimpleName(), file] }.collect()

    ch_dbsnp_index = Channel.fromPath("${params.dbsnp}.tbi")
        .map{ file -> [file.getSimpleName(), file] }.collect()

    ch_PGx_genes = Channel.fromList(params.pgx_genes)

    PYPGX_CREATEREGIONS()

    BCFTOOLS_VIEW(
        ch_dbsnp
            .join(ch_dbsnp_index)
        .map{ meta, vcf, tbi -> [[id: meta], vcf,tbi] },
        PYPGX_CREATEREGIONS.out.bed,
        [],
        []
    )

    GATK4_HAPLOTYPECALLER (
        ch_bams_meta
            .combine(PYPGX_CREATEREGIONS.out.bed)
            .map{ meta, bam, bai, intervals ->  [meta,bam,bai,intervals, []]},
        ch_genome_fasta,
        ch_genome_fasta_index,
        ch_genome_dict,
        ch_dbsnp,
        ch_dbsnp_index
    )

    GATK4_GENOTYPEGVCFS(
        GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi)
            .combine(PYPGX_CREATEREGIONS.out.bed)
            .map{ meta, vcf, tbi, intervals -> [meta, vcf, tbi, intervals, []] },
        ch_genome_fasta.map{ meta, file -> [file] },
        ch_genome_fasta_index.map{ meta, file -> [file] },
        ch_genome_dict.map{ meta, file -> [file] },
        ch_dbsnp.map{ meta, file -> [file] },
        ch_dbsnp_index.map{ meta, file -> [file] }
    )

    PYPGX_CREATEINPUTVCF(
        ch_bams_meta,
        ch_genome_fasta
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Variant QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Run mosdepth on PGx gene regions to calculate average coverage
    MOSDEPTH(
        ch_bams_meta
            .combine(PYPGX_CREATEREGIONS.out.bed),
        ch_genome_fasta
    )

    //todo terminate workflow when sample coverage is too low

    PARSE_MOSDEPTH(
        MOSDEPTH.out.summary_txt
    )

    FILTER_PYPGX_VCF(
        PYPGX_CREATEINPUTVCF.out.vcf
    )

    FILTER_GATK_VCF(
        GATK4_GENOTYPEGVCFS.out.vcf
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PyPGx pipeline
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    PYPGX_PREPAREDEPTHOFCOVERAGE(
        ch_bams_meta
    )

    PYPGX_COMPUTECONTROLSTATISTICS(
        ch_bams_meta
    )

    PYPGX_RUNNGSPIPELINE(
        FILTER_PYPGX_VCF.out.vcf
            .join(FILTER_PYPGX_VCF.out.tbi)
            .join(PYPGX_PREPAREDEPTHOFCOVERAGE.out.coverage)
            .join(PYPGX_COMPUTECONTROLSTATISTICS.out.control_stats)
            .combine(ch_PGx_genes),
        params.pypgx_resource_bundle
    )



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Star allele calling translation table
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_excelsheet = Channel.fromPath(params.phenotypes_excel)
        .map{ file -> [[file.getSimpleName()], file] }
        .collect()
    
    CALL_STARALLELES(
        FILTER_GATK_VCF.out.vcf
            .join(FILTER_GATK_VCF.out.tbi)
            .combine(ch_PGx_genes)
            .map{meta, vcf, tbi, gene -> [meta, gene, vcf,tbi]}
            .join(PYPGX_RUNNGSPIPELINE.out.outdir, by: [0,1]),
        ch_excelsheet,
        BCFTOOLS_VIEW.out.vcf //dbSNP subset
            .collect(),
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Finalize / QA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    COMBINERESULTS(
        PYPGX_RUNNGSPIPELINE.out.outdir
            .map { meta, gene, dir -> [gene, dir]}
            .groupTuple()
            .join(CALL_STARALLELES.out.csv
                  .groupTuple())
    )

    SV_QA(COMBINERESULTS.out.csv)

    COV_QA(
        COMBINERESULTS.out.csv,
        PARSE_MOSDEPTH.out.average_pg_coverage
            .map {meta, val -> val}
            .collect()
    )

    // Softare versions
    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
    ch_versions = ch_versions.mix(PYPGX_CREATEINPUTVCF.out.versions)
    ch_versions = ch_versions.mix(PYPGX_PREPAREDEPTHOFCOVERAGE.out.versions)
    ch_versions = ch_versions.mix(PYPGX_COMPUTECONTROLSTATISTICS.out.versions)
    ch_versions = ch_versions.mix(PYPGX_RUNNGSPIPELINE.out.versions)
    ch_versions = ch_versions.mix(CALL_STARALLELES.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
    ch_versions = ch_versions.mix(FILTER_PYPGX_VCF.out.versions)
    ch_versions = ch_versions.mix(FILTER_GATK_VCF.out.versions)
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    // MultiQC
    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )
 }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    def analysis_id = params.outdir.split('/')[-1]
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: analysis_id,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "PG Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/multiqc_report.html")
    } else {
        def subject = "PG Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
