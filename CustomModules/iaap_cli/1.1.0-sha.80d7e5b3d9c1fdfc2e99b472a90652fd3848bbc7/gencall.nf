process GenCall {
    // Raw idat to Genotypes
    tag {"GenCall ${identifier}"}
    label 'IAAP_CLI_1_1_0'
    label 'IAAP_CLI_1_1_0_GenCall'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = true 

    input:
        tuple(val(identifier), val(array_id), path(grn_idat), path(red_idat))

    output:
        tuple(val(identifier), val(array_id), path("${identifier}.gtc"))

    script:
        gender_estimate_threshold = params.estimate_gender? -0.1 : 0.9
        cluster_file = params.estimate_gender? params.gender_estimate_file : params.cluster_file
        """
        mkdir ${array_id}
        cp ${grn_idat} ${array_id}
        cp ${red_idat} ${array_id}

        ${params.iaap_path} gencall \
            --idat-folder ${array_id} \
            --gender-estimate-call-rate-threshold ${gender_estimate_threshold} \
            --output-gtc \
            ${params.optional} \
            ${params.bead_pool_manifest_file} \
            ${cluster_file} \
            .
        """
}