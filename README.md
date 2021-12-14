- [1. DxNextflowPG](#1-dxnextflowpg)
  - [1.1. Git branching model](#11-git-branching-model)
  - [1.2. Workflow 'blocks'](#12-workflow-blocks)
  - [1.3. Running PG workflow](#13-running-pg-workflow)
  - [1.4. Genome builds / versions.](#14-genome-builds--versions)
  - [1.5. Assets and pytest](#15-assets-and-pytest)
  - [1.6. Limitations](#16-limitations)
- [2. Software dependencies](#2-software-dependencies)
  - [2.1. python environment](#21-python-environment)
  - [2.2. IAAP](#22-iaap)
  - [2.3. BafRegress](#23-bafregress)
  - [Docker files](#docker-files)
  - [2.4. Nextflow Modules](#24-nextflow-modules)
  - [2.5. Install Nextflow](#25-install-nextflow)
  - [2.6. References](#26-references)
    - [2.6.1. Workflow references](#261-workflow-references)
 
# 1. DxNextflowPG
Repo with scripts to process pharmacogenetic data.
## 1.1. Git branching model
In this repository, a branching model is applied: [a-successful-git-branching-model](https://nvie.com/posts/a-successful-git-branching-model/)

## 1.2. Workflow 'blocks'
This workflow contains the following blocks:
- Retrieving and reading raw intensity data files (idat).
- Genotyping

## 1.3. Running PG workflow
```bash
nextflow run beadarray.nf -c beadarray.config --idat_path idat_dir_path --outdir output_dir_path --email email [-profile slurm|mac]
```
## 1.4. Genome builds / versions.
GRCh38 is used.

## 1.5. Assets and pytest
Several custom python scripts used within the workflow are added to the assets directory.
Pytests are included for some of these custom python scripts.

I choose to set-up the [tests-outside-the-application-code](https://docs.pytest.org/en/6.2.x/goodpractices.html#choosing-a-test-layout-import-rules).
Therefore, I used the package set up.
- Install package in repos directory: 
    ```bash 
    python3 -m venv venv
    . venv/bin/activate
    pip install -e .
    ```

To run pytest for all tests with verbosity and print statements shown in console, use:
```bash
pytest --flake8 -v -s tests/
``` 
To run pytest for a specific class with verbosity, use for example:
```bash
pytest --flake8 -v tests/test_genotype_to_phenotype.py::TestInput
```
The scope of tests you want to run, can be selected. [specifying-tests-or-selecting-tests](https://docs.pytest.org/en/6.2.x/usage.html#specifying-tests-selecting-tests)

Some pytests require testdata and is included in this repository alongside the actual test code.

To test the total pytest coverage, use [pytest-cov](https://pytest-cov.readthedocs.io/en/latest/readme.html)
## 1.6. Limitations
Not supported yet:
- Structural variants, such as CNV, INDELs (> 1 bp) 
- INDELs of 1 bp with more than 2 alleles.

Other remarks:
- ENSEMBL REST API
  - Since ensembl REST Api is only available for a single GRCh38, version cannot be specified. Assumed that patches will not effect the sites of interests.

TODO: retrieve MAF for GRCh38 or remove BAFREGRESS.
# 2. Software dependencies
##  2.1. python environment
- Login HPC
- Go to repository ./assets/
- Run:
    ```bash
    python3 -m venv venv
    . venv/bin/activate
    pip install --upgrade pip
    python3 -m pip install -r requirements.txt
    ```

## 2.2. IAAP
- Retrieve tar.gz from https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html
- Move tar.gz to HPC /hpc/diaggen/software/tools/
- Login HPC
- Run:
    ```bash
    cd /hpc/diaggen/software/tools/
    tar -xzvf iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz
    ```

Docker image (using Docker Desktop locally). Replace `version` to match with Dockerfile argument.
```bash
cd ./CustomModules/IlluminaGtcToVcf/version/
docker build -t umcugenbioinf/illumina_gtctovcf:version -f Dockerfile .
docker push umcugenbioinf/illumina_gtctovcf:version
```
## 2.3. BafRegress
HPC installation
```bash
cd /hpc/diaggen/software/tools/
wget https://genome.sph.umich.edu/w/images/d/d5/BafRegress.tar.gz
tar -xzvf BafRegress.tar.gz -C bafregress_1_0_0
cd bafregress_1_0_0
git clone https://github.com/jpdna/VCFtoFinalReportForBafRegress.git
```

Docker image (using Docker Desktop locally). Replace `version` to match with Dockerfile argument.
```bash
cd ./tools/IlluminaGtcToVcf/version/
docker build -t umcugenbioinf/illumina_gtctovcf:version -f Dockerfile .
docker push umcugenbioinf/illumina_gtctovcf:version
```
## Docker files
Build docker image for software dependencies. 
- [Install Docker Desktop](https://docs.docker.com/desktop/mac/apple-silicon/)
    ```bash
    docker build -t organization_or_username/toolname:version -f path_to_dockerfile
    docker push organization_or_username/toolname:version
    ```
## 2.4. Nextflow Modules
Get Nextflow Modules
```bash
git submodule add git@github.com:UMCUGenetics/NextflowModules.git
```
Update Nextflow Modules
```bash
git submodule update --init --recursive
```

## 2.5. Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

## 2.6. References
### 2.6.1. Workflow references
Use assets/create_workflow_reference_files.py to create required .bed and .yaml files.
The .bed file is required to create an intervallist.

Note: it is possible to use docker and/or Nextflow as well (except BedToIntervalList), see scripts in:
- NextflowModules/Picard/2.26.4--hdfd78af_0/CreateSequenceDictionary.nf
- NextflowModules/Picard/2.26.4--hdfd78af_0/CreateExtendedIlluminaManifest.nf

To generate the reference files, run following commands at HPC computenode. Replace the parameters with actual strings. 

***Generate .dict file***

Make sure the --GENOME_ASSEMBLY matches with the reference fasta. for example: 'GRCh38'.

```bash
java -jar -Xmx4G /hpc/local/CentOS7/cog_bioinf/picard-tools-2.5.0/picard.jar \
CreateSequenceDictionary \
REFERENCE=path_to_ref_fasta \
OUTPUT=path_and_filename_of_output_dict \
GENOME_ASSEMBLY=genome_assembly
```

***Generate intervallist***

Make sure the sequence_dictionary matches with the one used in the workflow.
```bash
./assets/migrate_translation_table_to_input.R --blabla

java -jar -Xmx4G /hpc/local/CentOS7/cog_bioinf/picard-tools-2.5.0/picard.jar \
BedToIntervalList \
I=path_to_generated_bed_file  \
O=path_and_filename_of_output_intervallist \
SEQUENCE_DICTIONARY=path_to_reference_genome_dict \
UNIQUE=true
```

***Generate extended manifest file***
```bash
set +u; env - PATH="$PATH" SINGULARITYENV_TMP="$TMP" SINGULARITYENV_TMPDIR="$TMPDIR" singularity exec \
-B "$PWD" -B /hpc:/hpc -B $TMPDIR:$TMPDIR \
/hpc/diaggen/software/singularity_cache/quay.io-biocontainers-picard-2.26.4--hdfd78af_0.img /bin/bash -c "\
picard -Xmx4G \
CreateExtendedIlluminaManifest \
--TMP_DIR $TMPDIR \
--INPUT path_and_filename_of_manifest_csv \
--OUTPUT path_and_filename_of_output_extended_manifest_csv>\
--REFERENCE_SEQUENCE path_to_ref_fasta \
--REPORT_FILE path_and_filename_of_output_report \
--CLUSTER_FILE path_and_filename_of_clusterfile \
--MAX_RECORDS_IN_RAM 100000"
```