process PYPGX_CREATEINPUTVCF {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypgx:0.25.0--pyh7e72e81_0':
        'biocontainers/pypgx:0.25.0--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    val(pypgx_gene_list)
    val(assembly_version)


    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def assembly = assembly_version ?: "GRCh38"
    def prefix = "${meta.id}"
    def pypgx_genes = pypgx_gene_list.join(' ')
    """
    pypgx create-input-vcf \\
        --genes ${pypgx_genes}
        --assembly ${assembly} \\
        ${prefix}_variants.vcf.gz \\
        ${fasta} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgx: \$(echo \$(pypgx -v 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
