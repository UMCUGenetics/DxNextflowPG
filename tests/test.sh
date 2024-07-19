#!/usr/bin/env sh


kgp_vcf=$1
dbsnp_vcf=$2
regions_file=$3
genotype_dir=$4
kgp_subset=./variants_regions.vcf
dbsnp_subset=./dbSNP_regions.vcf

if [ ! -f "$regions_file" ]; then
    singularity run \
        -B /hpc:/hpc \
        -B $TMPDIR:$TMPDIR \
        /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img \
            pypgx create-regions-bed \
                --assembly GRCh38 \
                --add-chr-prefix \
            > $regions_file

fi

if [ ! -f "$kgp_subset" ]; then
    singularity run \
        -B /hpc:/hpc \
        -B $TMPDIR:$TMPDIR \
        /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
            bcftools view \
                -R $regions_file \
                -O v \
                "$kgp_vcf" \
            > $kgp_subset
fi

if [ ! -f "$dbsnp_subset" ]; then
        singularity run \
        -B /hpc:/hpc \
        -B $TMPDIR:$TMPDIR \
        /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
            bcftools view \
                -R $regions_file \
                -O v \
                $dbsnp_vcf \
            > $dbsnp_subset
fi

if [ ! -d $genotype_dir ]; then
    echo "Sampling Genotypes"
    source testenv/bin/activate
    python random_sampling.py $kgp_subset $genotype_dir

fi
echo "Compressing VCFs"
cd $genotype_dir

for vcf in ./*.vcf; do
   singularity run \
        -B /hpc:/hpc \
        -B $TMPDIR:$TMPDIR \
        /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
            bgzip $vcf

   singularity run \
        -B /hpc:/hpc \
        -B $TMPDIR:$TMPDIR \
        /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
            tabix "$vcf".gz

   echo "Importing VCF"
   singularity run \
       -B /hpc:/hpc \
       -B $TMPDIR:$TMPDIR \
       /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img \
            pypgx import-variants \
                --assembly GRCh38 \
                CYP2C19 "$vcf".gz $(basename $vcf .vcf)_imported_variants.zip
done

# PyPGX star alleles
for vcf in ./*.vcf.gz; do
   singularity run \
       -B /hpc:/hpc \
       -B $TMPDIR:$TMPDIR \
       /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img \
            pypgx predict-alleles \
                $(basename $vcf .vcf.gz)_imported_variants.zip $(basename $vcf .vcf.gz)_star_alleles.zip
done


