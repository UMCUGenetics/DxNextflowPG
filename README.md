

# 1. Table of Contents
- [1. Table of Contents](#1-table-of-contents)
- [2. DxNextflowPG](#2-dxnextflowpg)
  - [2.1. Git branching model](#21-git-branching-model)
  - [2.2. Workflow 'blocks'](#22-workflow-blocks)
  - [2.3. Running PG workflow](#23-running-pg-workflow)
  - [2.4. Assets and pytest](#24-assets-and-pytest)
- [3. Software dependencies](#3-software-dependencies)
  - [3.1. python environment](#31-python-environment)
  - [3.2. IAAP](#32-iaap)
  - [3.3. BafRegress](#33-bafregress)
  - [3.4. Nextflow Modules](#34-nextflow-modules)
  - [3.5. Install Nextflow](#35-install-nextflow)
  - [3.6. References](#36-references)
    - [3.6.1. Workflow references](#361-workflow-references)
    - [3.6.2. Pytest references](#362-pytest-references)
 
# 2. DxNextflowPG
Repo with scripts to process pharmacogenetic data.
## 2.1. Git branching model
In this repository, a branching model is applied: [a-successful-git-branching-model](https://nvie.com/posts/a-successful-git-branching-model/)

## 2.2. Workflow 'blocks'
This workflow contains the following blocks:
- Retrieving and reading raw intensity data files (idat).
- Genotyping

## 2.3. Running PG workflow
```bash
nextflow run beadarray.nf -c beadarray.config --idat_path <idat_dir_path> --outdir <output_dir_path> --email <email> [-profile slurm|mac]
```

## 2.4. Assets and pytest
Several custom python scripts used within the workflow are added to the assets directory.
Pytests are included for some of these custom python scripts.

I choose to set-up the [tests-outside-the-application-code](https://docs.pytest.org/en/6.2.x/goodpractices.html#choosing-a-test-layout-import-rules).
Therefore, I used the package set up.
- Install package in repos directory: `pip install -e .`

To run pytest for all tests with verbosity and print statements shown in console, use:
```bash
pytest -v -s tests/
``` 

The scope of tests you want to run, can be selected. [specifying-tests-or-selecting-tests](https://docs.pytest.org/en/6.2.x/usage.html#specifying-tests-selecting-tests)
# 3. Software dependencies
##  3.1. python environment
- Login HPC
- Go to repository ./assets/
- Run:
```bash
python3 -m venv venv
. venv/bin/activate
python3 -m pip install -r requirements.txt
```

## 3.2. IAAP
- Retrieve tar.gz from https://emea.support.illumina.com/downloads/iaap-genotyping-cli.html
- Move tar.gz to HPC /hpc/diaggen/software/tools/
- Login HPC
- Run:
```bash
cd /hpc/diaggen/software/tools/
tar -xzvf iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz
```

## 3.3. BafRegress
```bash
cd /hpc/diaggen/software/tools/
wget https://genome.sph.umich.edu/w/images/d/d5/BafRegress.tar.gz
tar -xzvf BafRegress.tar.gz -C bafregress_1_0_0
cd bafregress_1_0_0
git clone https://github.com/jpdna/VCFtoFinalReportForBafRegress.git
```

## 3.4. Nextflow Modules
Get Nextflow Modules
```bash
git submodule add git@github.com:UMCUGenetics/NextflowModules.git
```
Update Nextflow Modules
```bash
git submodule update --init --recursive
```

## 3.5. Install Nextflow
```bash
mkdir tools && cd tools
curl -s https://get.nextflow.io | bash
```

## 3.6. References
### 3.6.1. Workflow references
To generate correct references, run following commands at HPC computenode. 
Replace all <> with strings. Make sure the sequence_dictionary matches with the one used in the workflow.
```bash
./assets/migrate_translation_table_to_input.R --blabla

java -jar -Xmx4G /hpc/local/CentOS7/cog_bioinf/picard-tools-2.5.0/picard.jar \
BedToIntervalList \
I=<path_to_generated_bed_file>  \
O=<path_and_filename_of_output_intervallist> \
SEQUENCE_DICTIONARY=<path_to_reference_genome_dict> \
UNIQUE=true
```

### 3.6.2. Pytest references
bgzip -c ./tests/test_data/fake_000000000000_R00C00.vcf > ./tests/test_data/fake_000000000000_R00C00.vcf.gz
tabix -p vcf ./tests/test_data/fake_000000000000_R00C00.vcf.gz