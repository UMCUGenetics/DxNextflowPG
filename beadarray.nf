#!/usr/bin/env nextflow
nextflow.preview.dsl=2

def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Raw idat to Genotypes
    AutoCall(analysis_id, params.chip_well_barcode, params.green_idat_path, params.red_idat_path)

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
    tag {"AutoCall ${analysis_id}"}
    label 'AutoCall'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = true 

    input:
        tuple(analysis_id, assay_id, path(green_idat), path(red_idat))
    
    script:
        """
        ${params.iaap_path} gencall \
        ${params.bead_pool_manifest_file} \
        ${params.cluster_file} \
        . \
        --idat-folder ${assay_id} \
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
        echo 'DxNextflowWes' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'Dx_tracks' >> repository_version.log
        git --git-dir=${params.dxtracks_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'ExonCov' >> repository_version.log
        git --git-dir=${params.exoncov_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'ExomeDepth' >> repository_version.log
        git --git-dir=${params.exomedepth_path}/../.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'TrendAnalysis' >> repository_version.log
        git --git-dir=${params.trend_analysis_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}
