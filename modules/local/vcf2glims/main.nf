process VCF2GLIMS {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/umcugenetics/vcf2glims:1.0.0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.analysis_id}_${meta.id}"

    """
    vcf2glims.py \\
        $args \\
        $meta.analysis_id \\
        $vcf \\
        > ${prefix}.csv
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
