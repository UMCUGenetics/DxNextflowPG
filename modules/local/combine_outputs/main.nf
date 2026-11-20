process COMBINERESULTS {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(pypgx_dirs)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.zip"), emit: zip

    when:
    task.ext.when == null || task.ext.when

    script:
    def suffix = task.ext.suffix ?: "--output_suffix _mqc.csv"
    """
    combine_outputs.py \\
        --pypgx_dirs ${pypgx_dirs} \\
        --gene ${meta.id} \\
        ${suffix}
    """

    stub:
    """
    touch ${meta.id}_results.zip
    touch ${meta.id}.csv
    """
}
