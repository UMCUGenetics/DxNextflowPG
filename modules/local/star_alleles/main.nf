process CALL_STARALLELES {
    debug true
    errorStrategy 'terminate'
    tag "$meta.id"
    label "process_single"


    container "docker://ghcr.io/umcugenetics/star_caller:v1.0.1"
    
    input:
    tuple val(meta), path(vcf), path(tbi), val(pgx_gene)
    path(excel)
    path(dbSNP)

    output:
    tuple val(pgx_gene), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = "${meta.id}_${pgx_gene}_alleles.csv"

    """

    echo ${meta}
    echo ${vcf}
    echo ${tbi}
    echo ${pgx_gene}

    echo ${excel}
    echo ${dbSNP}

    ${projectDir}/modules/local/star_alleles/star_calling/call_star_alleles_excel.py \\
        --vcf ${vcf} \\
        --dbSNP ${dbSNP} \\
        --translation_table ${excel} \\
        --pgx_gene ${pgx_gene} \\
        -o ${outfile}
    """
}
