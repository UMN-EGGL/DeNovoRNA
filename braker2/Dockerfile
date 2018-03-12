FROM ubuntu:16.04

WORKDIR /workspace

ARG BRAKER2_PATH="/opt/braker2/"
ARG BRAKER2_DL_URL="http://bioinf.uni-greifswald.de/augustus/binaries/BRAKER2.tar.gz" 
ARG AUGUSTUS_TMP_PATH="/tmp/augustus-build"
ARG AUGUSTUS_DL_URL="http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz"

# Add braker2 and augustus sources 
ADD $BRAKER2_DL_URL /tmp
RUN mkdir -p $AUGUSTUS_TMP_PATH
ADD $AUGUSTUS_DL_URL $AUGUSTUS_TMP_PATH
ENV PATH="$BRAKER2_PATH:$AUGUSTUS_PATH:${PATH}"

# Set up build environment and compile
RUN apt-get update -qy
RUN mkdir -p $BRAKER2_PATH
RUN mkdir -p $AUGUSTUS_TMP_PATH
RUN tar xzvf /tmp/BRAKER2.tar.gz -C $BRAKER2_PATH
RUN tar xzvf $AUGUSTUS_TMP_PATH/augustus.current.tar.gz -C $AUGUSTUS_TMP_PATH
RUN apt-get -qy install perl-modules build-essential libboost-iostreams-dev zlib1g-dev libgsl-dev libsqlite3-dev libboost-graph-dev libsuitesparse-dev liblpsolve55-dev libbamtools-dev 
RUN cd $AUGUSTUS_TMP_PATH/augustus/ && make && make install
RUN cpan install File::Compare Scalar::Util::Numeric Parallel::ForkManager File::Which
