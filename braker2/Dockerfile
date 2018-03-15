FROM ubuntu:16.04

WORKDIR /workspace

# Temporary build paths, as well as install paths
ARG TMP_PATH="/tmp"
ARG BRAKER2_PATH="/opt/braker2/"
ARG BRAKER2_TMP_PATH="/tmp/braker2"
ARG BRAKER2_DL_URL="http://bioinf.uni-greifswald.de/augustus/binaries/BRAKER2.tar.gz" 
ARG AUGUSTUS_TMP_PATH="/tmp/augustus-build"
ARG AUGUSTUS_DL_URL="http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz"
ARG GENEMARK_TMP_PATH="/tmp/genemark"
ARG GENEMARK_DL_URL="http://topaz.gatech.edu/GeneMark/tmp/GMtool_vmkL3/gm_et_linux_64.tar.gz"
ARG GENEMARK_PATH="/opt/genemark"
ENV PATH "$GENEMARK_PATH:$BRAKER2_PATH:${PATH}"


# Add braker2, GeneMark, and augustus sources
RUN mkdir -p $BRAKER2_PATH $BRAKER2_TMP_PATH $AUGUSTUS_TMP_PATH $GENEMARK_PATH $GENEMARK_TMP_PATH
ADD $BRAKER2_DL_URL  $BRAKER2_TMP_PATH
ADD $AUGUSTUS_DL_URL $AUGUSTUS_TMP_PATH
ADD $GENEMARK_DL_URL $GENEMARK_TMP_PATH

# Do a few things here.
# 1. Extract archives, and move BRAKER2 scripts to their path.
# 2. Install dependencies for BRAKER2 and AUGUSTUS.
# 3. Set up build environment and compile AUGUSTUS.
RUN apt-get update -qy
RUN tar xzvf $BRAKER2_TMP_PATH/BRAKER2.tar.gz -C $BRAKER2_TMP_PATH
RUN cp -r $BRAKER2_TMP_PATH/BRAKER*/* $BRAKER2_PATH
RUN tar xzvf $AUGUSTUS_TMP_PATH/augustus.current.tar.gz -C $AUGUSTUS_TMP_PATH
RUN tar xzvf $GENEMARK_TMP_PATH/gm_et_linux_64.tar.gz -C $GENEMARK_TMP_PATH
RUN cp $GENEMARK_TMP_PATH/gm_et_linux_64/gmes_petap/gm_key ~/.gm_key
RUN cp -r $GENEMARK_TMP_PATH/gm_et_linux_64/gmes_petap/* $GENEMARK_PATH/
RUN apt-get -qy install perl-modules build-essential libboost-iostreams-dev zlib1g-dev libgsl-dev libsqlite3-dev libboost-graph-dev libsuitesparse-dev liblpsolve55-dev libbamtools-dev bamtools samtools ncbi-blast+
RUN cd $AUGUSTUS_TMP_PATH/augustus/ && make && make install
RUN cpan install File::Compare Scalar::Util::Numeric Parallel::ForkManager File::Which YAML Hash::Merge Logger::Simple
