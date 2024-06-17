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

include { validateParameters; paramsSummaryLog } from 'plugin/nf-validation'
log.info paramsSummaryLog(workflow)
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { extractFastqPairFromDir } from './modules/local/utils/fastq.nf'

include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS } from './modules/nf-core/gatk4/genotypegvcfs/main'
include { MOSDEPTH } from './modules/nf-core/mosdepth/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { SAMBAMBA_MARKDUP } from './modules/nf-core/sambamba/markdup/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { SEQKIT_SPLIT2 } from './modules/nf-core/seqkit/split2/main'
include { VCF2GLIMS } from './modules/local/vcf2glims/main'
include { VERIFYBAMID_VERIFYBAMID2 } from './modules/nf-core/verifybamid/verifybamid2/main'

include { CONTROLFREEC_FREEC } from './modules/nf-core/controlfreec/freec/main'
include { MANTA_GERMLINE } from './modules/nf-core/manta/germline/main'
include { DELLY_CALL } from './modules/nf-core/delly/call/main'
include { SV2 } from './modules/local/SV2/main'


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

    ch_bams = Channel.fromPath("${params.bam_path}/*.bam")

    ch_genome_fasta
        .map{name, path -> [id: name]}
        .set{ch_genome_meta}

    ch_bams_meta = ch_bams.map{ data -> [[id: data.getBaseName()], data] }


    // CONTROLFREEC_FREEC(
    //     ch_bams_meta,
    //     params.genome_fasta,
    //     "${params.genome_fasta}.fai",
    //     params.genome_chrfiles
    // )




    /*
     Manta
    */
    ch_bams_idx = Channel.fromPath("${params.bam_path}/*.bai")
    ch_manta_target = Channel.fromPath("$projectDir/assets/manta_target.bed.gz")
    ch_manta_target_index = Channel.fromPath("$projectDir/assets/manta_target.bed.gz.tbi")

    manta_input = ch_bams_meta
        .merge(ch_bams_idx)
        .combine(ch_manta_target)
        .combine(ch_manta_target_index)
        .combine(Channel.fromPath(params.manta_config))

    MANTA_GERMLINE(
        manta_input,
        ch_genome_fasta,
        ch_genome_fasta_index
    )

    delly_exclude = Channel.fromPath("$projectDir/assets/human.hg19.excl.tsv")

    delly_input = ch_bams_meta
        .merge(ch_bams_idx)
        .combine(delly_exclude)

    DELLY_CALL(
        delly_input,
        ch_genome_fasta,
        ch_genome_fasta_index
    )


    DELLY_CALL.out.bcf.combine(MANTA_GERMLINE.out.candidate_sv_vcf, by:0).view()
    SV2 (
        ch_bams_meta,
        ch_genome_fasta,
        DELLY_CALL.out.bcf.combine(MANTA_GERMLINE.out.candidate_sv_vcf, by:0)
    )


    // Variant calling
    // GATK4_HAPLOTYPECALLER(
    //     ch_bam_bai.combine(ch_intervals).map{ meta, bam, bai, intervals -> [meta, bam, bai, intervals, []] },
    //     ch_genome_fasta, ch_genome_fasta_index, ch_genome_dict, ch_dbsnp, ch_dbsnp_index
    // )
    // GATK4_GENOTYPEGVCFS(
    //     GATK4_HAPLOTYPECALLER.out.vcf.join(GATK4_HAPLOTYPECALLER.out.tbi).combine(ch_intervals).map{
    //         meta, vcf, tbi , intervals -> [meta, vcf, tbi, intervals, []]
    //     },
    //     ch_genome_fasta.map{ meta, file -> [file] },
    //     ch_genome_fasta_index.map{ meta, file -> [file] },
    //     ch_genome_dict.map{ meta, file -> [file] },
    //     ch_dbsnp.map{ meta, file -> [file] },
    //     ch_dbsnp_index.map{ meta, file -> [file] }
    // )

    // GLIMS output
    // VCF2GLIMS(GATK4_GENOTYPEGVCFS.out.vcf)

    // QC
    // FASTQC(ch_fastq)
    // MOSDEPTH(
    //     ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []] },
    //     ch_genome_fasta
    // )
    // VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_svd, Channel.empty().toList(), ch_genome_fasta.map{ meta, file -> [file] })

//     // Softare versions
//     ch_versions = channel.empty()
//     ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
//     ch_versions = ch_versions.mix(FASTQC.out.versions)
//     ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
//     ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
//     ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
//     ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)
//     ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
//     ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
//     ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)
//     CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

//     // MultiQC
//     ch_multiqc_files = Channel.empty()
//     ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(SAMBAMBA_MARKDUP.out.txt.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.collect{it[1]})
//     ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//     ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//     MULTIQC(
//         ch_multiqc_files.collect(),
//         ch_multiqc_config.toList(),
//         Channel.empty().toList(),
//         Channel.empty().toList()
//     )
 }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     def analysis_id = params.outdir.split('/')[-1]
//     // HTML Template
//     def template = new File("$baseDir/assets/workflow_complete.html")
//     def binding = [
//         runName: analysis_id,
//         workflow: workflow
//     ]
//     def engine = new groovy.text.GStringTemplateEngine()
//     def email_html = engine.createTemplate(template).make(binding).toString()

//     // Send email
//     if (workflow.success) {
//         def subject = "PG Workflow Successful: ${analysis_id}"
//         sendMail(to: params.email.trim(), subject: subject, body: email_html, attach: "${params.outdir}/QC/multiqc_report.html")
//     } else {
//         def subject = "PG Workflow Failed: ${analysis_id}"
//         sendMail(to: params.email.trim(), subject: subject, body: email_html)
//     }
// }
