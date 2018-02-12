FROM ubuntu:16.04

WORKDIR /workspace

ARG BRAKER2_PATH="/opt/braker2/"
ADD ./BRAKER_v2.0.5 $BRAKER2_PATH
ENV PATH="$BRAKER2_PATH:${PATH}"

RUN apt-get update -qy
RUN apt-get -qy install perl-modules build-essential
RUN cpan install File::Compare Scalar::Util::Numeric Parallel::ForkManager File::Which
