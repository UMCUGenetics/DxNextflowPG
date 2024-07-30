#!/usr/bin/env sh


# kgp_vcf=$1
# dbsnp_vcf=$2
# regions_file=$3
# genotype_dir=$4
# kgp_subset=./variants_regions.vcf
# dbsnp_subset=./dbSNP_regions.vcf

# if [ ! -f "$regions_file" ]; then
#     singularity run \
#         -B /hpc:/hpc \
#         -B $TMPDIR:$TMPDIR \
#         /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img \
#             pypgx create-regions-bed \
#                 --assembly GRCh38 \
#                 --add-chr-prefix \
#             > $regions_file

# fi

# if [ ! -f "$kgp_subset" ]; then
#     singularity run \
#         -B /hpc:/hpc \
#         -B $TMPDIR:$TMPDIR \
#         /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
#             bcftools view \
#                 -R $regions_file \
#                 -O v \
#                 "$kgp_vcf" \
#             > $kgp_subset
# fi

# if [ ! -f "$dbsnp_subset" ]; then
#         singularity run \
#         -B /hpc:/hpc \
#         -B $TMPDIR:$TMPDIR \
#         /hpc/diaggen/software/singularity_cache/quay.io-biocontainers-bcftools-1.20--h8b25389_0.img \
#             bcftools view \
#                 -R $regions_file \
#                 -O v \
#                 $dbsnp_vcf \
#             > $dbsnp_subset
# fi
p

source testenv/bin/activate

dummy_vcf_dir="dummy_vcfs"

mkdir -p output/
echo "Creating genotype dummy files"
python CreateGenotypeSimulationFiles.py \
    -i input/aug2023_translation_table_20230810.xlsx \
    -v input/sites_of_interest_GRCh38_dbsnp_146.vcf

mkdir -p $dummy_vcf_dir
echo "Converting dummy files to vcf"
for dummy_file in output/*.csv; do
   echo $dummy_file
   python dummyToVCF.py \
       $dummy_file > dummy_vcfs/$(basename $dummy_file).vcf

done

cd dummy_vcfs
for vcf in ./*.vcf; do
   echo "$vcf"
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

   pg_gene=$(basename $vcf | cut -f 1 -d '_')
   echo $pg_gene
   echo "Importing VCF"
   singularity run \
       -B /hpc:/hpc \
       -B $TMPDIR:$TMPDIR \
       /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img \
            pypgx import-variants \
                --assembly GRCh38 \
                $pg_gene "$vcf".gz $(basename $vcf .vcf)_"$pg_gene"_imported_variants.zip
done


# PyPGX star alleles

   singularity shell \
       -B /hpc:/hpc \
       -B $TMPDIR:$TMPDIR \
       /hpc/diaggen/software/singularity_cache/depot.galaxyproject.org-singularity-pypgx-0.25.0--pyh7e72e81_0.img

   for imported_vars in ./*imported_variants.zip; do
       pypgx predict-alleles \
           $imported_vars $(basename $imported_vars _imported_variants.zip)_star_alleles.zip
   done

   exit
   # rm $imported_vars # Cleanup tmp file
