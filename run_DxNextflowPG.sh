#!/bin/bash
set -euo pipefail

workflow_path='/home/cog/jvansteenbrugge/workdir/pypgx/DxNextflowPG/'

# Set input and output dirs
input=`realpath -e $1`
output=`realpath $2`
email=$3
bam_path=$4
genome=$5
optional_params=( "${@:6}" )


sbatch <<EOT
#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem 5G
#SBATCH --gres=tmpspace:10G
#SBATCH --job-name DxNextflowPG
#SBATCH -o log/slurm_DxNextflowPG.%j.out
#SBATCH -e log/slurm_DxNextflowPG.%j.err
#SBATCH --mail-user $email
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=diaggen

export NXF_JAVA_HOME='/hpc/diaggen/users/joris/software/tools/java/jdk'

/hpc/diaggen/users/joris/software/tools/nextflow/nextflow run $workflow_path \
    --input $input \
    --outdir $output \
    --email $email \
    --bam_path $bam_path \
    --genome_fasta $genome \
    -resume \
    -ansi-log false \
    -profile slurm \
    ${optional_params[@]:-""}


EOT
