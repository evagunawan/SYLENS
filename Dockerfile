# Specify the parent image
FROM ubuntu:focal as app

#Creating variable that is only accessible in build
ARG SYLENS_VER="0.1.0"

# Metadata
LABEL base.image="Ubuntu Focal 20.04"
LABEL dockerfile.version="1"
LABEL software="Sylens"
LABEL software.version="${SYLENS_VER}"
LABEL description="FASTQ downsampler and ASC II converter"
LABEL website="https://github.com/evagunawan/SYLENS"
LABEL license="GNU GPL3"
LABEL license.url="https://github.com/evagunawan/SYLENS/blob/main/LICENSE"
LABEL maintainer="Eva Gunawan"
LABEL maintainer.email="eva.gunawan@slh.wisc.edu"

# Getting dependencies 
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    git \
    python3 \
    python3-pip \
    gzip && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/* && \
    pip3 install argparse && \
    mkdir data 

# Cloning Sylens
RUN wget https://github.com/evagunawan/SYLENS/archive/main.tar.gz -O SYLENS.tgz && \
    tar -xvf SYLENS.tgz && \
    rm SYLENS.tgz && \
    mkdir SYLENS && \
    cd SYLENS && \
    chmod +x SYLENS && \
    cp SYLENS .

# Setting the work directory. Staphb uses /data always
WORKDIR /data

# Running the image
CMD ["Sylens_main.py", "--help"]

#Testing the image
FROM app as test

WORKDIR /test

RUN wget -O SRR026674.fastq.gz "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRR026674&format=fastq"

CMD ["sylens_main.py", "SRR026674.fastq.gz", "-s", "10"]

