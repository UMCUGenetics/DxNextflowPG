process SV2 {
    tag "$meta.id"
    label 'error_retry'

    container 'ghcr.io/umcugenetics/sv2:1.5'


    input:
    tuple val(meta),
        path(bam), path(bai),
        path(delly_vcf), path(delly_vcf_tbi),
        path(manta_vcf), path(manta_vcf_tbi),
        path(freec_vcf), path(freec_vcf_tbi),
        path(SNP_calls), path(SNP_calls_idx)

    output:
    tuple val(meta), path("sv2_preprocessing/", type: "dir")
    tuple val(meta), path("sv2_features/", type: "dir")
    tuple val(meta), path("sv2_genotypes/", type: "dir")

    script:
    def args = task.ext.args ?: ""
    def vcf  = "-v ${delly_vcf} ${manta_vcf} ${freec_vcf}"
    def bam  = bam ? "-i ${bam}" : ""
    def snp  = "-snv $SNP_calls"
    def sample_id = meta?["id"]

    """
    touch sample.ped
    #echo "#Family ID\tIndividual ID\tPaternal ID\tMaternal ID\tSex\tPhenotype" >> sample.ped
    # Gender is set to 1 (male) for sv2 to function, but is irrelevant for CYP2D6/7 CNVs
    echo "0\t${sample_id}\t0\t0\t1\tNA" >> sample.ped

    sv2 \\
    ${bam} \\
    ${vcf} \\
    ${snp} \\
    -ped sample.ped
    """

}
