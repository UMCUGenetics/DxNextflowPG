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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER as FILTER_PYPGX_VCF } from './modules/nf-core/bcftools/filter/main'
include { COMBINERESULTS                      } from './modules/local/combine_outputs/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS         } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                             } from './modules/nf-core/multiqc/main'
include { PYPGX_CREATEREGIONS                 } from './modules/local/pypgx/create_regions/main'
include { PYPGX_CREATEINPUTVCF                } from './modules/nf-core/pypgx/createinputvcf/main'
include { PYPGX_PREPAREDEPTHOFCOVERAGE        } from './modules/nf-core/pypgx/preparedepthofcoverage/main'
include { PYPGX_COMPUTECONTROLSTATISTICS      } from './modules/nf-core/pypgx/computecontrolstatistics/main'
include { PYPGX_RUNNGSPIPELINE                } from './modules/nf-core/pypgx/runngspipeline/main'

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

    // Coverage depth for each pharmacogene, relevant for SV prediction
    PYPGX_PREPAREDEPTHOFCOVERAGE(
        ch_bams_meta,
        ch_PGx_genes.collect(),
        ch_assembly_version
    )

    // Control statistics to compare pharmacogenes with household gene Vitamin D Receptor
    PYPGX_COMPUTECONTROLSTATISTICS(
        ch_bams_meta,
        "VDR",
        ch_assembly_version
    )


    PYPGX_RUNNGSPIPELINE(
        PYPGX_CREATEINPUTVCF.out.vcf
            .join(PYPGX_CREATEINPUTVCF.out.tbi)
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




    // Softare versions
    ch_versions = channel.empty()
    ch_versions = ch_versions.mix(PYPGX_CREATEINPUTVCF.out.versions)
    ch_versions = ch_versions.mix(PYPGX_PREPAREDEPTHOFCOVERAGE.out.versions)
    ch_versions = ch_versions.mix(PYPGX_COMPUTECONTROLSTATISTICS.out.versions)
    ch_versions = ch_versions.mix(PYPGX_RUNNGSPIPELINE.out.versions)
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(COMBINERESULTS.out.csv.collect())
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )
 }

