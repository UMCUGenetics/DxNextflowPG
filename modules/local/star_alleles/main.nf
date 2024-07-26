process CALL_STARALLELES {
    tag "$meta.id"
    label "process_single"

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-344874846f44224e5f0b7b741eacdddffe895d1e:d3fff24ee1297b4c3bcef48354c2a30f0c82007a-0':
    //     'biocontainers/mulled-v2-344874846f44224e5f0b7b741eacdddffe895d1e:d3fff24ee1297b4c3bcef48354c2a30f0c82007a-0' }"

    container "docker://quay.io/biocontainers/mulled-v2-344874846f44224e5f0b7b741eacdddffe895d1e:d3fff24ee1297b4c3bcef48354c2a30f0c82007a-0"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(excel)
    tuple val(meta3), path(dbSNP)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = "${meta}_alleles.csv"
    """
    pip install openpyxl # for testing purposes
    
    call_star_alleles_excel.py \\
        --vcf ${vcf} \\
        --dbSNP ${dbSNP} \\
        --translation_table ${excel} \\
        -o ${outfile}
    """
}
