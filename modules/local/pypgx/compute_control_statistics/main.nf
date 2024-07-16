process PYPGX_COMPUTECONTROLSTATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    val(control_gene)
    val(assembly_version)

    output:
    tuple val(meta), path('*control_statistics*.zip'), emit: control_stats
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def control = control_gene ?: "VDR"
    def assembly = assembly_version ?: "GRCh38"
    def prefix = "${meta.id}"

    """
    pypgx compute-control-statistics \\
        --assembly ${assembly} \\
        ${control} \\
        ${prefix}_control_statistics_${control}.zip \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

}
