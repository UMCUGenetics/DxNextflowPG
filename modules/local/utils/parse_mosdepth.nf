process PARSE_MOSDEPTH {
    tag "${meta.id}"

    input:
    tuple val(meta), path(summary_txt)

    output:
    tuple val(meta), path("*.cov"), emit: average_pg_coverage

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${summary_txt} | grep "total_region" | cut -f 4 > ${prefix}.cov
    """
}
