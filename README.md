[![GitHub Actions CI Status](https://github.com/UMCUGenetics/DxNextflowPG/workflows/nf-core%20CI/badge.svg)](https://github.com/UMCUGenetics/DxNextflowPG/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/UMCUGenetics/DxNextflowPG/workflows/nf-core%20linting/badge.svg)](https://github.com/UMCUGenetics/DxNextflowPG/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# DxNextflowPG

Nextflow workflow to analyse pharmacogenetics NGS data.

## UMCU HPC

```
export NXF_JAVA_HOME='/hpc/diaggen/software/tools/jdk-18.0.2.1/'
tools/nextflow run DxNextflowPG --input <input_path> --outdir <output_path> --email <email> -resume
```

## Development

### Install nf-core and Nextflow

```
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install nf-core

curl -s https://get.nextflow.io | bash
mv nextflow venv/bin/
```

## Testing

To run nf-test:

``` sh
git submodule update --init --recursive
nf-test test . --tag=local #optionally set a profile with --profile <docker|singularity|etc.>
```
