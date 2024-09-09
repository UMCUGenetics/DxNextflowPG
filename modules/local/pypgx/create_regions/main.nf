process PYPGX_CREATEREGIONS {
    tag "pypgx_createregions"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    output:
    path("pypgx_regions.bed"), emit: bed
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def assembly = task.ext.assembly_version ?: "GRCh38"
    def genes = "--genes ${task.ext.pgx_genes.join(' ')}" ?: ''
    """
    pypgx create-regions-bed \\
        ${genes} \\
        --assembly ${assembly} \\
        ${args} \\
        > pypgx_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
