process GenCall {
    // Raw idat to Genotypes
    tag {"GenCall ${identifier}"}
    label 'GenCall'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = true 

    input:
        tuple(val(identifier), val(array_id), path(grn_idat), path(red_idat))

    output:
        tuple(val(identifier), val(array_id), path("${identifier}.gtc"))

    script:
        """
        mkdir ${array_id}
        cp ${grn_idat} ${array_id}
        cp ${red_idat} ${array_id}

        ${params.iaap_path} gencall \
        ${params.bead_pool_manifest_file} \
        ${params.cluster_file} \
        . \
        --idat-folder ${array_id} \
        ${params.gender_gencall_option} \
        --output-gtc
        """
}