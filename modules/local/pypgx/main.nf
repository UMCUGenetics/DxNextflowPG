
process pypgx_prepare {
    tag "$meta.id"
    container ""

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam_file), path(bam_bai)
    tuple val(meta2), path(genome_fasta)
    val(control_gene)
    val(pypgx_gene_list)
    val(assembly_version)


    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path('*coverage.zip'), emit: coverage
    tuple val(meta), path('*control_statistics*.zip'), emit: control_stats
    path("versions.yml"), emit: versions


    script:
    def control = control_gene ?: "VDR"
    def assembly = assembly_version ?: "GRCh38"
    def pypgx_genes = pypgx_gene_list.join(' ')
    def prefix = "${meta.id}"
    """
    pypgx create-input-vcf \\
        --genes ${pypgx_genes} \\
        --assembly ${assembly} \\
        ${prefix}_variants.vcf.gz \\
        ${genome_fasta} \\
        $bam_file

    pypgx prepare-depth-of-coverage \\
        --genes ${pypgx_genes} \\
        --assembly ${assembly} \\
        ${prefix}_coverage.zip \\
        $bam_file

    pypgx compute-control-statistics \\
        --assembly ${assembly} \\
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

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(coverage), path(control_stats), val(pgx_gene)
    path(resource_bundle)
    val(assembly_version)


    output:
    tuple val(pgx_gene), path("*pypgx_output"), emit: outdir
    path("versions.yml"), emit: versions


    script:
    def prefix = "${meta.id}_${pgx_gene}"
    def assembly = assembly_version ?: "GRCh38"
    """

    export MPLCONFIGDIR="/tmp/"
    export PYPGX_BUNDLE=${resource_bundle}/

    pypgx run-ngs-pipeline \\
        --assembly ${assembly} \\
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

process combine_results {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(pgx_gene), path(pypgx_dirs), path(excel_caller_csvs)

    output:
    path("*.csv"), emit: csv
    path("*.zip"), emit: zip

    script:
    """
    #!/usr/bin/env python

    from pypgx.sdk import utils as sdk
    import pandas as pd

    pypgx_archives = [sdk.Archive.from_file(pypgx_output+'/results.zip').data
            for pypgx_output in '${pypgx_dirs}'.split()]
    combined_sample_data = pd.concat(pypgx_archives)

    # Pypygx specific outputting

    # Meta data is the same for all samples across the same pgx_gene
    metadata = sdk.Archive.from_file('${pypgx_dirs}'.split()[0]+'/results.zip').metadata
    merged_output = sdk.Archive(metadata, combined_sample_data)
    merged_output.to_file("${pgx_gene}_results.zip")


    excel_alleles = pd.concat(
    [pd.read_csv(excel_output, sep='\t', index_col=0)
        for excel_output in '${excel_caller_csvs}'.split()])


    
    merged = combined_sample_data.join(excel_alleles)

    # For human readable output
    merged.to_csv(open('${pgx_gene}.csv', 'w'), sep='\t')
    """

}
