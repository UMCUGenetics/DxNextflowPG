#!/usr/bin/env nextflow
nextflow.preview.dsl=2

include GtcToVcf as PICARD_GtcToVcf from './NextflowModules/Picard/2.25.5/GtcToVcf.nf' params(
    bead_pool_manifest_file: "${params.bead_pool_manifest_file}",
    cluster_file: "${params.cluster_file}", 
    extended_chip_manifest_file: "${params.extended_chip_manifest_file}",
    genome: "${params.genome}",
    optional: ""
    )

def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Raw idat to Genotypes (VCF format)
    AutoCall(analysis_id, params.chip_well_barcode, params.green_idat_path, params.red_idat_path) // TODO: change analysis_id to assay_id
    PICARD_GtcToVcf(params.chip_well_barcode, AutoCall.out)

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
        sendMail(to: params.email, subject: subject, body: email_html, attach: "${params.outdir}/QC/${analysis_id}_multiqc_report.html")
    } else {
        def subject = "PG Workflow Failed: ${analysis_id}"
        sendMail(to: params.email, subject: subject, body: email_html)
    }
}

process AutoCall {
    tag {"AutoCall ${assay_id}"}
    label 'AutoCall'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = true 

    output:
        path("${assay_id}.gtc")

    input:
        val(assay_id)
        val(array_id)
        path(green_idat)
        path(red_idat)
    
    script:
        """
        rm -rf ${array_id}
        mkdir ${array_id}
        cp ${green_idat} ${array_id}
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
