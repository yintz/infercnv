# Docker file for inferCNV

FROM ubuntu:trusty
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### Locations
ENV SRC /usr/local/src

######################
## Dependencies and Tools
######################
##############
## Helper tools
RUN apt-get update && \
    apt-get install -y git libxml2-dev r-cran-xml libcurl4-openssl-dev

# Install dependencies for r-base source build
RUN apt-get build-dep -y r-base

## Install R 3.2.1
WORKDIR ${SRC}
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz && \
    apt-get install -y tcl tcl-dev tk-dev tk && \
    tar xf R-3.2.1.tar.gz

WORKDIR ${SRC}/R-3.2.1
RUN ./configure --prefix=/usr/local --enable-R-shlib --with-tcltk && \
    make && \
    make install

# Install R libraries
RUN echo "options(repos=\"http://cran.rstudio.com/\")" > ${SRC}/install_pkgs.R && \
    echo "install.packages(c(\"GMD\", \"ape\", \"RColorBrewer\", \"optparse\", \"logging\"), repos=\"http://cran.us.r-project.org\", dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    Rscript ${SRC}/install_pkgs.R

# Pull in inferCNV source
WORKDIR ${SRC}
RUN git clone https://github.com/broadinstitute/inferCNV
