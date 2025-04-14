#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UMCUGenetics/DxNextflowPG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UMCUGenetics/DxNextflowPG
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Validate and log parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
log.info paramsSummaryLog(workflow)
validateParameters()

nextflow.enable.moduleBinaries = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER as FILTER_PYPGX_VCF } from './modules/nf-core/bcftools/filter/main'
include { COMBINERESULTS                      } from './modules/local/combine_outputs/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { MOSDEPTH                            } from './modules/nf-core/mosdepth/main'
include { MULTIQC                             } from './modules/nf-core/multiqc/main'
include { PYPGX_CREATEINPUTVCF                } from './modules/nf-core/pypgx/createinputvcf/main'
include { PYPGX_PREPAREDEPTHOFCOVERAGE        } from './modules/nf-core/pypgx/preparedepthofcoverage/main'
include { PYPGX_COMPUTECONTROLSTATISTICS      } from './modules/nf-core/pypgx/computecontrolstatistics/main'
include { PYPGX_RUNNGSPIPELINE                } from './modules/nf-core/pypgx/runngspipeline/main'
include { SAMTOOLS_INDEX                      } from './modules/nf-core/samtools/index/main'
include { SV_QA                               } from './modules/local/SV_QA/main'
include { VERIFYBAMID_VERIFYBAMID2            } from './modules/nf-core/verifybamid/verifybamid2/main'


include { MAPPING } from './subworkflows/local/mapping/main'

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

    ch_bwa_index = Channel.fromPath("${params.bwa_index}*")
        .map{ file -> [file.getSimpleName(), file] }
        .groupTuple()
        .collect()

    ch_bams_meta = Channel.fromFilePairs(
        "${params.bam_path}/*.{bam,bai}",
        checkIfExists: true) {
            file -> file.name.replaceAll(/.bam|.bai$/,'') }
        .map{ meta, bam_index -> [['id': meta], bam_index[0], bam_index[1]] }

    ch_PGx_resource_bundle = Channel.fromPath(params.pypgx_resource_bundle)
        .map{ file -> [[id: file.getSimpleName()], file] }
        .collect()

    ch_PGx_genes = Channel.fromList(params.pgx_genes)
    ch_assembly_version = Channel.value(params.assembly_version)

    ch_svd = Channel.fromPath(["${params.svd_ud}", "${params.svd_mu}", "${params.svd_bed}"]).collect()




    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Optional read mapping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    MAPPING(
        ch_bwa_index,
        params.fastq_path
    )

    // Merge bam files from (optionally) mapped samples into the channel with bam
    // files that were already mapped
    // ch_bams_meta
        // .concat(MAPPING.out.bam
                    // .join(MAPPING.out.bai)
        // )

    ch_bam_bai = ch_bams_meta


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    PyPGx pipeline
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // pypgx variant calling
    PYPGX_CREATEINPUTVCF(
        ch_bams_meta,
        ch_genome_fasta,
        ch_PGx_genes.collect(),
        ch_assembly_version
    )


    FILTER_PYPGX_VCF(
        PYPGX_CREATEINPUTVCF.out.vcf
    )

    // Coverage depth for each pharmacogene, relevant for SV prediction
    PYPGX_PREPAREDEPTHOFCOVERAGE(
        ch_bams_meta,
        ch_PGx_genes.collect(),
        ch_assembly_version
    )

    // Control statistics to compare pharmacogenes with household gene Vitamin D Receptor
    PYPGX_COMPUTECONTROLSTATISTICS(
        ch_bams_meta,
        params.pgx_control_gene,
        ch_assembly_version
    )


    PYPGX_RUNNGSPIPELINE(
        FILTER_PYPGX_VCF.out.vcf
            .join(FILTER_PYPGX_VCF.out.tbi)
            .join(PYPGX_PREPAREDEPTHOFCOVERAGE.out.coverage)
            .join(PYPGX_COMPUTECONTROLSTATISTICS.out.control_stats)
            .combine(ch_PGx_genes),
        ch_PGx_resource_bundle,
        ch_assembly_version
    )


    // Group samples and produce a summary table for each pharmacogene
    COMBINERESULTS(
        PYPGX_RUNNGSPIPELINE.out.outdir
            .map { meta, gene, dir -> [gene, dir]}
            .groupTuple()
    )






    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    SV_QA(COMBINERESULTS.out.csv)

    if(params.verify_bamid) {
        VERIFYBAMID_VERIFYBAMID2(
            ch_bams_meta,
            ch_svd,
            Channel
                .empty()
                .toList(),
            ch_genome_fasta
                .map{ meta, file -> [file] }
        )
    }



    MOSDEPTH(
        ch_bams_meta
            .map{ meta, bam, bai -> [meta, bam, bai, []] },
        ch_genome_fasta
    )


    // Softare versions
    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_versions = ch_versions.mix(PYPGX_CREATEINPUTVCF.out.versions)
    ch_versions = ch_versions.mix(PYPGX_PREPAREDEPTHOFCOVERAGE.out.versions)
    ch_versions = ch_versions.mix(PYPGX_COMPUTECONTROLSTATISTICS.out.versions)
    ch_versions = ch_versions.mix(PYPGX_RUNNGSPIPELINE.out.versions)
    if(params.verify_bamid) {
        ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)
    }
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(COMBINERESULTS.out.csv.collect())
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})

    if(params.verify_bamid) {
        ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.collect{it[1]})
    }

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
        def subject = "PGx Workflow Successful: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/multiqc_report.html")
    } else {
        def subject = "PGx Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
