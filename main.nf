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
    Validate parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; } from 'plugin/nf-validation'
validateParameters()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import modules/subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { SAMBAMBA_MARKDUP } from './modules/nf-core/sambamba/markdup/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Reference file channels
    // ch_genome = Channel.fromPath(params.genome).map {genome -> [genome.getSimpleName(), genome] }

    // Input channel
    ch_fastq = extractFastqPairFromDir(params.input)
    ch_fastq.view()
    ch_fastq

    // Mapping
    BWAMEM2_MEM(ch_fastq, ch_genome, true)
    SAMBAMBA_MARKDUP(BWAMEM2_MEM.out.bam.map{ meta, bam -> [ meta - meta.subMap('read_group'), bam ] }.groupTuple())

    // Variant calling
    // GATK_HaplotypeCallerGVCF
    // GATK_GenotypeGVCFs

    // GLIMS output
    // VCF2GLIMS


    // QC
    FASTQC(ch_fastq)
    // Mosdepth
    // VerifyBamID2

    // MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        Channel.empty().toList(),
        Channel.empty().toList()
    )

}

def flowcellLaneFromFastq(path) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab

    // parse first line of a FASTQ file (optionally gzip-compressed)
    // and return the flowcell id and lane number.
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven or eight fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String machine
    int run_nr
    String fcid
    int lane

    machine = fields[0]
    run_nr = fields[1].toInteger()
    fcid = fields[2]
    lane = fields[3].toInteger()

    [fcid, lane, machine, run_nr]
}

def extractFastqPairFromDir(dir) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    dir = dir.tokenize().collect{"$it/**_R1_*.fastq.gz"}
    Channel
    .fromPath(dir, type:'file')
    .ifEmpty { error "No R1 fastq.gz files found in ${dir}." }
    .filter { !(it =~ /.*Undetermined.*/) }
    .map { r1_path ->
        def fastq_files = [r1_path]
        def sample_id = r1_path.getSimpleName().split('_')[0]
        def r2_path = file(r1_path.toString().replace('_R1_', '_R2_'))
        if (r2_path.exists()) {
            fastq_files.add(r2_path)
        } else {
            exit 1, "R2 fastq.gz file not found: ${r2_path}."
        }
        def (flowcell, lane) = flowcellLaneFromFastq(r1_path)
        def rg_id = "${sample_id}_${flowcell}_${lane}"

        [['id': sample_id, 'read_group': rg_id], fastq_files]
    }
}
