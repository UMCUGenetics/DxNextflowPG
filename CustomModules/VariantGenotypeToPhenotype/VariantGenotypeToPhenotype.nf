process VariantGenotypeToPhenotype {
    // Custom process to translate (sets of) variant genotype to a phenotype.
    tag {"VariantGenotypeToPhenotype ${identifier}"}
    label 'VariantGenotypeToPhenotype'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(identifier), path(vcf_file), path(vcf_idx_file)) // should be compressed VCF with tabix index.

    output:
        path("${identifier}_genotypes.txt")

    script:
        """
        source ${baseDir}/assets/venv/bin/activate
        python ${baseDir}/assets/variant_genotype_to_phenotype.py \
        ${vcf_file} \
        ${identifier} \
        ${params.translation_table} \
        --output_prefix ${identifier}_genotypes
        """
}
