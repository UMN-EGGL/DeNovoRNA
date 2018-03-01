FROM ubuntu:16.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "Tools to remap RNASeq data for gene annotation"

RUN apt-get -y update && apt-get install -y \
    curl wget git gcc build-essential zlib1g-dev

RUN mkdir -p /src/ && \
    cd /src && \
    git clone https://github.com/alexdobin/STAR.git && \ 
    cd STAR/source && \
    make STAR

ENV PATH=/src/STAR/source/:${PATH}

WORKDIR /root

VOLUME ["/scratch/sshfs/data_release/"]

ENTRYPOINT [ "STAR" ]


