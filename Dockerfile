# Specify the parent image
FROM ubuntu:focal as app
#Creating variable that is only accessible in build
ARG SYLENS_VER="0.1.0"
# Metadata
LABEL base.image="Ubuntu Focal 20.04"
LABEL dockerfile.version="1"
LABEL software="Sylens"
LABEL software.version="${SYLENS_VER}"
LABEL description="FASTQ downsampler and ASCII converter"
LABEL website="https://github.com/evagunawan/SYLENS"
LABEL license="GNU GPL3"
LABEL license.url="https://github.com/evagunawan/SYLENS/blob/main/LICENSE"
LABEL maintainer="Eva Gunawan"
LABEL maintainer.email="eva.gunawan@slh.wisc.edu"
# Getting dependencies 
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    python3 \
    python3-pip \
    gzip && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/* 
RUN pip3 install argparse \
    biopython>=1.81 && \
    mkdir data
# Grabbing Sylens from github
RUN wget https://github.com/evagunawan/SYLENS/archive/main.tar.gz && \
    tar -xvf main.tar.gz && \
    rm main.tar.gz
# Add SYLENS directory to PATH
ENV PATH="/SYLENS-main:$PATH"
# Setting the work directory. Staphb uses /data
WORKDIR /data
# Running the image, default command if nothing is entered
CMD ["Sylens_main.py", "--help"]
#Testing the image
FROM app as test
WORKDIR /test
RUN wget -O sample1_R1.fastq.gz "https://github.com/nf-core/test-datasets/raw/viralrecon/illumina/amplicon/sample1_R1.fastq.gz"
RUN Sylens_main.py sample1_R1.fastq.gz -s 10