
process pypgx_prepare {
    tag "$meta.id"
    container ""

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biocontainers/pypgx:0.25.0--pyh7e72e81_0' :
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam_file), path(bam_bai)
    tuple val(meta2), path(genome_fasta)

    val(control_gene)


    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path('*coverage.zip'), emit: coverage
    tuple val(meta), path('*control_statistics*.zip'), emit: control_stats
    path("versions.yml"), emit: versions


    script:
    def control = control_gene ?: "VDR"
    def prefix = "${meta.id}"
    """
    pypgx create-input-vcf \\
        ${prefix}_variants.vcf.gz \\
        ${genome_fasta} \\
        $bam_file

    pypgx prepare-depth-of-coverage \\
        ${prefix}_coverage.zip \\
        $bam_file

    pypgx compute-control-statistics \\
        ${control} \\
        ${prefix}_control_statistics_${control}.zip \\
        $bam_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}

process pypgx_run_ngs {
    tag "$meta.id"

    container "biocontainers/pypgx:0.25.0--pyh7e72e81_0"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(coverage), path(control_stats), val(pgx_gene)
    path(resource_bundle)


    output:
    tuple val(pgx_gene), path("*pypgx_output"), emit: outdir
    path("versions.yml"), emit: versions


    script:
    def prefix = "${meta.id}_${pgx_gene}"
    """

    export MPLCONFIGDIR="/tmp/"
    export PYPGX_BUNDLE=${resource_bundle}/

    pypgx run-ngs-pipeline \\
        ${pgx_gene} \\
        ${prefix}_pypgx_output/ \\
        --variants ${vcf} \\
        --depth-of-coverage ${coverage} \\
        --control-statistics ${control_stats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}

process create_pypgx_output_table {

    container "biocontainers/pypgx:0.25.0--pyh7e72e81_0"

    input:
    tuple val(pgx_gene), path(output_dirs)

    output:
    path("*.csv"), emit: csv

    script:
    """
    #!/usr/bin/env python

    from pypgx.sdk import utils as sdk
    import pandas as pd

    pypgx_archives = [sdk.Archive.from_file(pypgx_output+'/results.zip').data
            for pypgx_output in '${output_dirs}'.split()]

    pd.concat(pypgx_archives).to_csv(open('${pgx_gene}.csv', 'w'), sep='\t')
    """
}
