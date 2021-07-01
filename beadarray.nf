#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// Retrieve input data files modules
include extractIdatPairFromDir from './NextflowModules/Utils/idat.nf'

// Genotyping modules
include GtcToVcf as PICARD_GtcToVcf from './NextflowModules/Picard/2.25.5/GtcToVcf.nf' params(
    bead_pool_manifest_file: "${params.bead_pool_manifest_file}",
    cluster_file: "${params.cluster_file}", 
    extended_chip_manifest_file: "${params.extended_chip_manifest_file}",
    genome: "${params.genome}",
    optional: ""
    )

// Contamination modules
include VcfToAdpc as PICARD_VcfToAdpc from './NextflowModules/Picard/2.25.5/VcfToAdpc.nf' params(optional: "")
include VerifyIDIntensity from './NextflowModules/VerifyIDIntensity/0.0.1--hc90279e_1/VerifyIDIntensity.nf'
include CreateVerifyIDIntensityContaminationMetricsFile as PICARD_VerifyIDToMetrics from './NextflowModules/Picard/2.25.5/CreateVerifyIDIntensityContaminationMetricsFile.nf'
include BafRegress from './NextflowModules/BafRegress/1.0.0/BafRegress.nf' 

// Filter modules
include SelectVariants as GATK_SelectVariants from './NextflowModules/GATK/4.1.3.0/SelectVariants.nf' params(genome:"${params.genome}", optional: "--intervals ${params.intervals_of_interest}")

// Retrieve input data files
def idat_files = extractIdatPairFromDir(params.idat_path) // [sample_id, array_id, grn_path, red_path]

def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Raw idat to Genotypes (VCF format)
    AutoCall(idat_files) 
    PICARD_GtcToVcf(AutoCall.out.map{sample_id, array_id, gtc_file -> [sample_id, gtc_file]})
    
    // Contamination
    BafRegress(PICARD_GtcToVcf.out)
    PICARD_VcfToAdpc(PICARD_GtcToVcf.out) // sample_id, vcf, vcf_index
    VerifyIDIntensity(PICARD_VcfToAdpc.out) 
    PICARD_VerifyIDToMetrics(VerifyIDIntensity.out)
    
    // Select sites of interest
    GATK_SelectVariants(PICARD_GtcToVcf.out.map{ sample_id, vcf_file, vcf_idx_file -> ["", vcf_file, vcf_idx_file, sample_id] }) 

    // Repository versions
    VersionLog()
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
        sendMail(to: params.email, subject: subject, body: email_html)
    } else {
        def subject = "PG Workflow Failed: ${analysis_id}"
        sendMail(to: params.email, subject: subject, body: email_html)
    }
}


process AutoCall {
    // Raw idat to Genotypes
    tag {"AutoCall ${sample_id}"}
    label 'AutoCall'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = true 

    input:
        tuple(val(sample_id), val(array_id), path(grn_idat), path(red_idat))

    output:
        tuple(val(sample_id), val(array_id), path("${sample_id}.gtc"))

    script:
        """
        mkdir ${array_id}
        cp ${grn_idat} ${array_id}
        cp ${red_idat} ${array_id}

        ${params.iaap_path} gencall \
        ${params.bead_pool_manifest_file} \
        ${params.cluster_file} \
        . \
        --idat-folder ${array_id} \
        ${params.gender_autocall_option} \
        --output-gtc
        """
}


process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowPG' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
