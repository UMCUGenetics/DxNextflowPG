process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  // Disable cache to force a new version log when restarting the workflow.
    container = 'umcugenbioinf/illumina_gtctovcf:1.2.1'  // Required to retrieve version of illumina gtctovcf.

    output:
        path('repository_version.log')

    script:
        """
        echo 'DxNextflowPG' >> repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo "IAAP_CLI" >> repository_version.log
        ${params.iaap_path} --version >> repository_version.log

        echo "BafRegress Estimate" >> repository_version.log
        python ${params.bafregress_path}/bafRegress.py estimate  --version >> repository_version.log
        echo "VCFtoFinalReportForBafRegress" >> repository_version.log
        git --git-dir=${params.bafregress_path}/VCFtoFinalReportForBafRegress/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo "Illumina GtcToVcf" >> repository_version.log
        ${params.gtctovcf_path}/gtc_to_vcf.py --version >> repository_version.log
        """
}