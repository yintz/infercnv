# Docker file for inferCNV

FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
# ENV R_VERSION 3.3.1-1xenial0
ENV R_VERSION 3.2.5-1xenial

### Locations
ENV SRC /usr/local/src

######################
## Dependencies and Tools
######################
##############
## Helper tools
RUN apt-get update && \
    apt-get install -y git libxml2-dev

## Install R
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | apt-key add - && \
    apt-get update && \
    apt-get install -y r-base-core=${R_VERSION} && \
    apt-get install -y r-recommended=${R_VERSION} && \
    apt-get install -y r-base=${R_VERSION} && \
    apt-get install -y r-base-dev=${R_VERSION}

RUN echo "source(\"https://bioconductor.org/biocLite.R\")" > ${SRC}/install_pkgs.R && \
    echo "biocLite(c(\"GMD\", \"ape\", \"RColorBrewer\", \"optparse\", \"logging\"), dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    Rscript ${SRC}/install_pkgs.R

# Pull in inferCNV source
WORKDIR ${SRC}
RUN git clone https://github.com/broadinstitute/inferCNV
