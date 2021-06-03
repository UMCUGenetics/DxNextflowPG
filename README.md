# DxNextflowPG
Repo with scripts to process pharmacogenetic data.

## Workflow 'blocks'
This workflow contains the following blocks:
- Retrieving and reading raw intensity data files (idat).
- Genotyping

## Settings / Technical details
### Git branching model
In this repository, a branching model is applied.[a-successful-git-branching-model](https://nvie.com/posts/a-successful-git-branching-model/)

### Get Nextflow Modules
```bash
git submodule update --init --recursive
```

### Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

### Running PG workflow
```bash
nextflow run beadarray.nf -c beadarray.config --idat_path <idat_dir_path> --outdir <output_dir_path> --email <email> [-profile slurm|mac]
```