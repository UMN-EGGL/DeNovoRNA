FROM ubuntu:16.04

WORKDIR /braker2

ADD ./BRAKER_v2.0.5 /braker2

RUN apt-get update -qy
RUN apt-get -qy install perl-modules build-essential
RUN cpan install File::Compare Scalar::Util::Numeric Parallel::ForkManager File::Which
