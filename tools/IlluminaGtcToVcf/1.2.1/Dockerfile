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
RUN pip install --no-cache-dir 'setuptools<58'
RUN pip install --no-cache-dir numpy pyvcf pysam

WORKDIR /tools/
ARG GTCTOVCF_VERSION="1.2.1"
RUN wget -q "https://github.com/Illumina/GTCtoVCF/archive/refs/tags/${GTCTOVCF_VERSION}.tar.gz" -O "GTCtoVCF-${GTCTOVCF_VERSION}.tar.gz" \
    && tar -xzf "GTCtoVCF-${GTCTOVCF_VERSION}.tar.gz" \
    && rm "GTCtoVCF-${GTCTOVCF_VERSION}.tar.gz"
ENV GTCTOVCF="GTCtoVCF-${GTCTOVCF_VERSION}/gtc_to_vcf.py"
