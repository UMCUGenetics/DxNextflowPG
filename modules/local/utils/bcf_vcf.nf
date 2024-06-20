process BCF2VCF {
    tag "$meta.id"

    container "biocontainers/bcftools:1.20--h8b25389_0"

    input:
    tuple val(meta), path(bcf)
    val(tool_name)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_tbi


    script:
    """
    bcftools \\
        view \\
        -O v \\
        ${bcf} \\
    | sed 's/INS/DUP/g' > ${meta.id}_${tool_name}_noins.vcf

    bgzip ${meta.id}_${tool_name}_noins.vcf
    tabix -p vcf ${meta.id}_${tool_name}_noins.vcf.gz
    """
}


process CLEAN_VCF {
    tag "$meta.id"

    container "biocontainers/bcftools:1.20--h8b25389_0"


    input:
    tuple val(meta), path(vcf)
    val(tool_name)

    output:
    tuple val(meta), path("*_noins.vcf.gz"), emit: vcf
    tuple val(meta), path("*_noins.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    gunzip -c ${vcf} \\
    | sed 's/INS/DUP/g' > ${meta.id}_${tool_name}_noins.vcf

    bgzip ${meta.id}_${tool_name}_noins.vcf
    tabix -p vcf ${meta.id}_${tool_name}_noins.vcf.gz
    """

}


process FREEC2VCF {
    tag "$meta.id"

    //container "python:3.10.14-alpine3.20"

    input:
    tuple val(meta), path(freec_CNV)

    output:
    tuple val(meta), path("*_noins.vcf.gz"), emit: vcf
    tuple val(meta), path("*_noins.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    python ${projectDir}/modules/local/utils/freecToVCF.py ${freec_CNV} > ${meta.id}_FREEC_noins.vcf

    bgzip ${meta.id}_FREEC_noins.vcf
    tabix -p vcf ${meta.id}_FREEC_noins.vcf.gz
    """
}
