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

include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { MOSDEPTH } from './modules/nf-core/mosdepth/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { VCF2GLIMS } from './modules/local/vcf2glims/main'
include { VERIFYBAMID_VERIFYBAMID2 } from './modules/nf-core/verifybamid/verifybamid2/main'
include { pypgx_prepare; pypgx_run_ngs; combine_results } from './modules/local/pypgx/main'

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

    // ch_genome_fasta_index = Channel.fromPath("${params.genome_fasta}.fai")
    //     .map{ file -> [file.getSimpleName(), file] }
    //     .collect()

    // ch_genome_dict = Channel.fromPath("${params.genome_dict}")
    //     .map{ file -> [file.getSimpleName(), file] }
    //     .collect()

    ch_bams_meta = Channel.fromFilePairs("${params.bam_path}/*.{bam,bai}", checkIfExists: true) { file -> file.name.replaceAll(/.bam|.bai$/,'') }
        .map{ meta, bam_index -> [['id': meta], bam_index[0], bam_index[1]] }


    // pypgx
    ch_PGx_genes = Channel.fromList(params.pgx_genes)

    pypgx_prepare(
        ch_bams_meta,
        ch_genome_fasta,
        params.pypgx_control_gene,
        params.assembly_version
    )

    pypgx_run_ngs(
        pypgx_prepare.out.vcf
            .join(pypgx_prepare.out.coverage)
            .join(pypgx_prepare.out.control_stats)
            .combine(ch_PGx_genes),
        params.pypgx_resource_bundle,
        params.assembly_version
    )

    combine_results(
        pypgx_run_ngs.out.outdir.groupTuple()
    )


    MOSDEPTH(
        ch_bam_idx_meta.map{ meta, bam, bai -> [meta, bam, bai, []] },
        ch_genome_fasta
    )

    // Softare versions
    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
//     ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)
//     ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
//     ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
//     ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)
    ch_versions = ch_versions.mix(pypgx_prepare.out.versions)
    ch_versions = ch_versions.mix(pypgx_run_ngs.out.versions)
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    // MultiQC
    ch_multiqc_files = Channel.empty()

    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.collect{it[1]})
    // ch_multiqc_files = ch_multiqc_files.mix(create_pypgx_output_table.out.csv.collect())
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
