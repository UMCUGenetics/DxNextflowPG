process UNZIP_PYPGX_OUTPUT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:16.02' :
        'biocontainers/p7zip:16.02' }"

    input:
    tuple val(meta), path(pypgx_outdir)

    output:
    tuple val(meta), path("${meta.id}_pypgx_multiqc.tsv"), emit: tsv

    script:
    def header_file = "${projectDir}/assets/pypgx_multiqc_header.txt"
    """
    7za \\
        x \\
        ${pypgx_outdir}/results.zip

    mv tmp*/data.tsv > ${meta.id}_pypgx_multiqc.tsv

    """

}
