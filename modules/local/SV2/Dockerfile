################## BASE IMAGE ######################
FROM  conda/miniconda2

LABEL version="1.5"
ADD Homo_sapiens.GRCh37.GATK.illumina.fasta /data/Hg37.fasta
ADD Homo_sapiens.GRCh37.GATK.illumina.fasta.fai /data/Hg37.fasta.fai

RUN sed -i 's/deb.debian.org/archive.debian.org/g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list && \
    sed -i '/security.debian.org/d' /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid-until && \
    apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

RUN conda install sv2 -c bioconda -c conda-forge
RUN yes | sv2 -download
RUN sv2 -hg19 /data/Hg37.fasta

#RUN chmod -R +rw /usr/local/lib/python2.7/site-packages/sv2/config
