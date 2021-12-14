################## BASE IMAGE ######################

FROM --platform=linux/amd64 python:3.9.7

################## METADATA ######################
LABEL base_image="python:3.9.7"
LABEL version="3.9.7"
LABEL software="Illumina_GTCtoVCF"
LABEL about.summary="GTC to VCF converter."
LABEL about.home="https://github.com/Illumina/GTCtoVCF"
LABEL about.documentation="https://github.com/Illumina/GTCtoVCF#readme"
LABEL about.tags="Microarray"

################## INSTALLATIONS ######################
RUN python -m pip install pip==21.2.4
RUN pip install --no-cache-dir 'setuptools<58'
RUN pip install --no-cache-dir numpy==1.21.4 pyvcf==0.6.8 pysam==0.18.0