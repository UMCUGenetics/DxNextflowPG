#!/usr/bin/env nextflow
nextflow.preview.dsl=2

def analysis_id = params.outdir.split('/')[-1]

workflow {
    // Raw idat to Genotypes
    AutoCall(analysis_id, params.chip_well_barcode, params.green_idat_path, params.red_idat_path)
    GtcToVcf(analysis_id, params.chip_well_barcode, AutoCall.out)

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

    output:
        path("${assay_id}.gtc")

    input:
        val(analysis_id)
        val(assay_id)
        path(green_idat)
        path(red_idat)
    
    script:
        """
        rm -rf ${assay_id}
        mkdir ${assay_id}
        cp ${green_idat} ${assay_id}
        cp ${red_idat} ${assay_id}

        ${params.iaap_path} gencall \
        ${params.bead_pool_manifest_file} \
        ${params.cluster_file} \
        . \
        --idat-folder ${assay_id} \
        ${params.gender_autocall_option} \
        --output-gtc
        """
}


process GtcToVcf {
    label 'PICARD_2_25_5'
    label 'PICARD_2_25_5_GtcToVcf'
    container = 'quay.io/biocontainers/picard:2.25.5--hdfd78af_0'
    shell = ['/bin/bash', '-euo', 'pipefail']
    
    input:
        val(analysis_id)
        val(assay_id)
        path(gtc_file)    
    
    output:
        tuple (assay_id, path("${assay_id}.vcf"),path("${assay_id}.vcf.tbi"), emit : genotyped_vcfs)

    script:
        """
        picard -Xmx${task.memory.toGiga()-4}G \
        GtcToVcf \
        TMP_DIR=\$TMPDIR \
        INPUT=${gtc_file} \
        OUTPUT=${assay_id}.vcf \
        --CLUSTER_FILE ${params.cluster_file} \
        --BPM_FILE ${params.bead_pool_manifest_file} \
        --MANIFEST ${params.chip_manifest_file} \
        --SAMPLE_ALIAS "${assay_id}" \
        --DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS true \
        --REFERENCE_SEQUENCE ${params.genome} \
        --MAX_RECORDS_IN_RAM 100000 \
        --CREATE_INDEX true
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
