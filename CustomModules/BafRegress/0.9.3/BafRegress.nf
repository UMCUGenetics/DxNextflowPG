process BafRegress {
    // Contamination estimate
    // B allele frequency regression models
    tag {"BafRegress ${sample_id}"}
    label 'BafRegress_0_9_3'
    label 'BafRegress_0_9_3_BafRegress'
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = 'umcugenbioinf/bafregress:0.9.3'

    input:
        tuple (val(sample_id), path(input_vcf), path(input_vcf_index))

    output:
        tuple (val(sample_id), path("${sample_id}_bafregress_report.txt"), emit: baffregress)

    script:
        """
        bcftools view ${input_vcf} | \
            python ${params.bafregress_path}/VCFtoFinalReportForBafRegress/parseVcfToBAFRegress.py > tmp_report.txt

        python ${params.bafregress_path}/bafRegress.py estimate \
            --freqfile ${params.maf_file} \
            tmp_report.txt \
            > ${sample_id}_bafregress_report.txt
        """
}