#!/bin/bash
set -euo pipefail

workflow_path='/home/cog/edejong2/repos/DxNextflowPG' # change to production

# Set input and output dirs
input=$(realpath -e $1)
data_file_type = $2
output=$(realpath $3)
email=$4
mkdir -p $output && cd $output
mkdir -p log

if { [ ${data_file_type@L} -e "idat" ] }; then
    data_path = "--idat_path ${data_file_type}"
elif { [ ${data_file_type@L} -e "gtc" ] }; then
    data_path = "--gtc_path ${data_file_type}"
else
    echo "Data file type ${data_file_type} not supported."
    exit 0
fi

if { [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; }; then
    echo "Workflow job not submitted, please check $output for 'workflow.status' files."
else
    touch workflow.running

    sbatch <<-EOT
    #!/bin/bash
    #SBATCH --time=08:00:00
    #SBATCH --nodes=1
    #SBATCH --mem 5G
    #SBATCH --gres=tmpspace:10G
    #SBATCH --job-name Nextflow_PG
    #SBATCH -o log/slurm_nextflow_pg.%j.out
    #SBATCH -e log/slurm_nextflow_pg.%j.err
    #SBATCH --mail-user $email
    #SBATCH --mail-type FAIL
    #SBATCH --export=NONE
    #SBATCH --account=diaggen

    module load Java/1.8.0_60

    /hpc/diaggen/users/ellen/software/nextflow run $workflow_path/beadarray.nf \
    -c $workflow_path/beadarray.config \
    $data_path \
    --outdir $output \
    --email $email \
    -profile slurm \
    -resume -ansi-log false

    if [ \$? -eq 0 ]; then
        echo "Nextflow done."

        # echo "Zip work directory"
        # find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

        # echo "Remove work directory"
        # rm -r work

        echo "Creating md5sum"
        find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

        echo "PG workflow completed successfully."
        rm workflow.running
        touch workflow.done

        exit 0
    else
        echo "Nextflow failed"
        rm workflow.running
        touch workflow.failed
        exit 1
    fi
EOT    
fi
