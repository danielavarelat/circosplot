FROM ubuntu:14.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    php5-mcrypt \
    python-pip
# configure locale, see https://github.com/rocker-org/rocker/issues/19
RUN \
    # Install linux dependencies
    echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections && \
    apt-get update && \
    apt-get install -yqq python \
    r-base r-base-core r-recommended && \
    \
    # Install Circos
    Rscript -e "withCallingHandlers(install.packages('RCircos', repos=c('http://cran.us.r-project.org', 'http://cran.mtu.edu')), warning = function(w) stop(w))" && \
    \
    # Install python3 dependencies
    pip install pysam pandas importlib requests subprocess argparse