process GtcToVcf {
    tag {"GtcToVcf ${identifier}"}
    label 'IlluminaGtcToVcf_1_2_1'
    label 'IlluminaGtcToVcf_1_2_1_GtcToVcf'
    container = 'umcugenbioinf/illumina_gtctovcf:1.2.1'
    shell = ['/bin/bash', '-euo', 'pipefail']
    
    input:
        tuple(val(identifier), path(gtc_file))
    
    output:
        tuple(val(identifier), path("${identifier}.vcf.gz"), path("${identifier}.vcf.gz.tbi"), emit : genotyped_vcfs)
    
    script:
        """
        ${params.gtctovcf_path}/GTCtoVCF \
        --gtc-paths ${gtc_file} \
        --manifest-file ${params.bead_pool_manifest_file} \
        --genome-fasta-file ${params.genome} \
        --output-vcf-path ${identifier}.vcf.gz \
        --include-attributes GT GQ BAF FT IGC LRR NORMX NORMY R THETA X Y \
        --log-file ${identifier}.log \
         ${params.optional}
        """
}