#!/bin/bash
set -euo pipefail

workflow_path='/home/cog/jvansteenbrugge/workdir/pypgx/DxNextflowPG/'

# Set input and output dirs
input=`realpath -e $1` # is currently ignored since already mapped bam files are used
output=`realpath $2`
email=`realpath $3`
bam_path=`realpath $4`
genome=`realpath $5`
pypgx_resources=`realpath $6`
optional_params=( "${@:7}" )
workdir=`pwd`

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
cd $workdir
/hpc/diaggen/users/joris/software/tools/nextflow/nextflow run $workflow_path \
    --input $input \
    --outdir $output \
    --email $email \
    --bam_path $bam_path \
    --genome_fasta $genome \
    --pypgx_resource_bundle $pypgx_resources \
    -resume \
    -ansi-log false \
    -profile slurm \
    ${optional_params[@]:-""}


EOT
