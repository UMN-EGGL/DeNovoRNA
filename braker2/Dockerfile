FROM ubuntu:16.04

WORKDIR /workspace

ARG BRAKER2_PATH="/opt/braker2/"
ARG BRAKER2_DL_URL="http://bioinf.uni-greifswald.de/augustus/binaries/BRAKER2.tar.gz" 

ADD $BRAKER2_DL_URL /tmp
ENV PATH="$BRAKER2_PATH:${PATH}"

RUN mkdir $BRAKER2_PATH
RUN tar xzvf /tmp/BRAKER2.tar.gz -C $BRAKER2_PATH
RUN mv $BRAKER2_PATH/BRAKER_v2*/* $BRAKER2_PATH
RUN apt-get update -qy
RUN apt-get -qy install perl-modules build-essential
RUN cpan install File::Compare Scalar::Util::Numeric Parallel::ForkManager File::Which
