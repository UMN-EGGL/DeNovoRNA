FROM ubuntu:16.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "Tools to remap RNASeq data for gene annotation"

# Install some base packages
RUN apt-get -y update && apt-get install -y \
    curl wget git gcc build-essential zlib1g-dev \
    vim parallel


# Install STAR
RUN mkdir -p /src/ && \
    cd /src && \
    git clone https://github.com/alexdobin/STAR.git && \ 
    cd STAR/source && \
    make STAR

# Install conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f
ENV PATH=/src/STAR/source:/root/miniconda3/bin:${PATH}
RUN conda update -n base conda

# Create a python environment
RUN conda create -y -n default python=3 
RUN /bin/bash -c "source activate default"


# Install deps
RUN conda install Cython numpy scipy pandas
RUN pip install --upgrade pip
RUN pip install git+https://github.com/LinkageIO/Minus80.git
RUN pip install git+https://github.com/LinkageIO/LocusPocus.git
RUN pip install git+https://github.com/schae234/htseq.git
RUN pip install ipython
RUN pip install ipdb

# Install Adapter Removal
RUN conda install -c maxibor adapterremoval2

WORKDIR /root

COPY data data
COPY scripts scripts

ENTRYPOINT [ "/bin/bash" ]


