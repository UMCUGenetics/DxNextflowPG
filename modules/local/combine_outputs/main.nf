process COMBINERESULTS {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(pgx_gene), path(pypgx_dirs), path(excel_caller_csvs)

    output:
    path("*.csv"), emit: csv
    path("*.zip"), emit: zip

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${moduleDir}/combine_outputs.py \\
        --pypgx_dirs ${pypgx_dirs} \\
        --gene ${pgx_gene} \\
        --excel_csvs ${excel_caller_csvs}
    """

    stub:
    """
    touch ${pgx_gene}_results.zip
    touch ${pgx_gene}.csv
    """
}
