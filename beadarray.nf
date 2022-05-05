#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Custom modules
include { VariantGenotypeToPhenotype } from './CustomModules/VariantGenotypeToPhenotype/VariantGenotypeToPhenotype.nf'
include { VersionLog } from './CustomModules/VersionLog/VersionLog.nf'

// Utils modules
include { extractIdatPairFromDir } from './NextflowModules/Utils/idat.nf'
include { ExportParams as Workflow_ExportParams } from './NextflowModules/Utils/workflow.nf'

// Genotyping modules
include { GenCall } from './CustomModules/iaap_cli/1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7/gencall.nf' params(
    bead_pool_manifest_file: "${params.bead_pool_manifest_file}",
    cluster_file: "${params.cluster_file}", 
    gender_estimate_file: "${params.gender_estimate_file}", 
    estimate_gender: false,
    iaap_path: ${params.iaap_path}",
    optional: ""
)
include { GtcToVcf as PICARD_GtcToVcf } from './NextflowModules/Picard/2.26.4--hdfd78af_0/GtcToVcf.nf' params(
    bead_pool_manifest_file: "${params.bead_pool_manifest_file}",
    cluster_file: "${params.cluster_file}", 
    extended_chip_manifest_file: "${params.extended_chip_manifest_file}",
    genome: "${params.genome}",
    optional: ""
)
include { GtcToVcf as Illumina_GtcToVcf } from './CustomModules/IlluminaGtcToVcf/1.2.1/GtcToVcf.nf'

// Contamination modules
include { VcfToAdpc as PICARD_VcfToAdpc } from './NextflowModules/Picard/2.26.4--hdfd78af_0/VcfToAdpc.nf' params(optional: "")
include { VerifyIDIntensity } from './NextflowModules/VerifyIDIntensity/0.0.1--hc90279e_1/VerifyIDIntensity.nf'
include { CreateVerifyIDIntensityContaminationMetricsFile as PICARD_VerifyIDToMetrics } from './NextflowModules/Picard/2.26.4--hdfd78af_0/CreateVerifyIDIntensityContaminationMetricsFile.nf'
include { BafRegress } from './tools/BafRegress/1.0.0/BafRegress.nf' 

// Quality metrics
include { CollectArraysVariantCallingMetrics as PICARD_VariantCallingMetrics } from './NextflowModules/Picard/2.26.4--hdfd78af_0/CollectArraysVariantCallingMetrics.nf' params(
    dbsnp: "$params.dbsnp", 
    call_rate_threshold: "$params.call_rate_threshold",
    output_prefix: "_VC_metrics"
)
include { CollectArraysVariantCallingMetrics as PICARD_VariantCallingMetrics_Intervals } from './NextflowModules/Picard/2.26.4--hdfd78af_0/CollectArraysVariantCallingMetrics.nf' params(
    dbsnp: "$params.dbsnp", 
    call_rate_threshold: "$params.call_rate_threshold",
    output_prefix: "_VC_metrics_subset"
)

// VCF manipulation modules
include { VariantFiltration as GATK_VariantFiltration } from './NextflowModules/GATK/4.2.1.0/VariantFiltration.nf' params(
    genome: "$params.genome", 
    compress: true,
    filter: "$params.gatk_filter",
    optional: ""
)

include { SelectVariants as GATK_SelectVariants } from './NextflowModules/GATK/4.2.1.0/SelectVariants.nf' params(
    genome:"$params.genome",
    compress: true,
    optional: "$params.gatk_select"
)

include { SelectVariants as GATK_SelectVariants_Intervals } from './NextflowModules/GATK/4.2.1.0/SelectVariants.nf' params(
    genome:"$params.genome", 
    compress: true,
    output_prefix: "_select_soi",
    optional: "--intervals $params.intervals_of_interest "
)

include { SelectVariants as GATK_SelectVariants_Autosomes } from './NextflowModules/GATK/4.2.1.0/SelectVariants.nf' params(
    genome:"$params.genome", 
    compress: true,
    output_prefix: "_autosomes",
    optional: "--exclude-filtered \
        --intervals chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    "
)

// Retrieve input data files
def idat_files = extractIdatPairFromDir(params.idat_path) // [sample_id, array_id, grn_path, red_path]

def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Raw idat to Genotypes (VCF format)
    GenCall(idat_files) 
    Illumina_GtcToVcf(GenCall.out.map{sample_id, array_id, gtc_file -> [sample_id, gtc_file]})
    // PICARD_GtcToVcf(GenCall.out.map{sample_id, array_id, gtc_file -> [sample_id, gtc_file]})
    
    // Contamination
    GATK_SelectVariants_Autosomes(Illumina_GtcToVcf.out)
    BafRegress(GATK_SelectVariants_Autosomes.out)
    PICARD_VcfToAdpc(Illumina_GtcToVcf.out)
    VerifyIDIntensity(PICARD_VcfToAdpc.out) 
    PICARD_VerifyIDToMetrics(VerifyIDIntensity.out)
    
    // VariantCallingMetrics
    PICARD_VariantCallingMetrics(Illumina_GtcToVcf.out)

    // Filter and select loci
    GATK_VariantFiltration(Illumina_GtcToVcf.out) 
    GATK_SelectVariants(GATK_VariantFiltration.out) 
    GATK_SelectVariants_Intervals(GATK_SelectVariants.out) 

    // VariantCallingMetrics on subset.
    PICARD_VariantCallingMetrics_Intervals(GATK_SelectVariants_Intervals.out)
    
    // Translate to phenotype
    VariantGenotypeToPhenotype(GATK_SelectVariants_Intervals.out)

    // Repository versions
    VersionLog()
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
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
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    } else {
        def subject = "PG Workflow Failed: ${analysis_id}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}
