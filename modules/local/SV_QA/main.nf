process SV_QA {
    tag "SV_QA"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-125c15aa3b8c154feff40f5b1165da9a24ec59fe:93482dc239e0f6857ec3db3e8f1d4c2abc36d676-0':
        'biocontainers/mulled-v2-125c15aa3b8c154feff40f5b1165da9a24ec59fe:93482dc239e0f6857ec3db3e8f1d4c2abc36d676-0' }"

    input:
    path(csv)

    output:
    path("*_SV_QA.csv"), emit: csv

    script:
    """
    frequency_dir='${projectDir}/assets/Dx_tracks/pharmacogenetics/frequencies'
    bn=\$(basename ${csv} .csv)

    sv_qa.py \\
        \$frequency_dir \\
        "${csv}" \\
        "\$bn"_SV_QA.csv
    """

    stub:
    """
    bn=\$(basename ${csv} .csv)
    touch "\$bn"_SV_QA.csv

    """
}
