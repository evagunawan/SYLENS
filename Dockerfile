# Specify the parent image
FROM ubuntu:focal as app

#Creating variable that is only accessible in build
ARG SYLENS_VER="0.1.0"

# Metadata
LABEL base.image="Ubuntu Focal 20.04"
LABEL dockerfile.version="1"
LABEL software="Sylens"
LABEL software.version="${SYLENS_VER}}"
LABEL description="FASTQ downsampler and ASC II converter"
LABEL website="https://github.com/evagunawan/SYLENS"
LABEL license="GNU GPL3"
LABEL license.url="https://github.com/evagunawan/SYLENS/blob/main/LICENSE"
LABEL maintainer="Eva Gunawan"
LABEL maintainer.email="eva.gunawan@slh.wisc.edu"

# Getting dependencies 
RUN apt-get update && apt-get -y install git \
    wget \
    pip \
    gzip \
    python3 && \
    apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Cloning git repository 
RUN git clone https://github.com/evagunawan/SYLENS.git

RUN pip install -r requirements.txt

# Setting the work directory. Staphb uses /data always
WORKDIR /data

# Running the image
CMD ["Sylens_main.py", "--help"]

#Testing the image
FROM app as test

WORKDIR /test