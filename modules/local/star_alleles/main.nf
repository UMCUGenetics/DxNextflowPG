process CALL_STARALLELES {

    tag "$meta.id"
    label "process_single"


    container "docker://ghcr.io/umcugenetics/star_caller:v1.0.1"
    
    input:
    tuple val(meta), path(vcf), path(tbi), val(pgx_gene)
    tuple val(meta2), path(excel)
    path(dbSNP)

    output:
    tuple val(pgx_gene), path("*.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = "${meta.id}_${pgx_gene}_alleles.csv"

    """
    ${projectDir}/modules/local/star_alleles/star_calling/call_star_alleles_excel.py \\
        --vcf ${vcf} \\
        --dbSNP ${dbSNP} \\
        --translation_table ${excel} \\
        --pgx_gene ${pgx_gene} \\
        -o ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CALL_STAR_ALLELES: \$(${projectDir}/modules/local/star_alleles/star_calling/call_star_alleles_excel.py -v)
    END_VERSIONS
    """
}
