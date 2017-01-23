# Docker file for inferCNV

# FROM ubuntu:xenial
FROM ubuntu:trusty
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
# ENV R_VERSION 3.3.1-1xenial0
# ENV R_VERSION 3.0.2-1ubuntu1
# ENV R_VERSION 3.2.5-1xenial

### Locations
ENV SRC /usr/local/src

######################
## Dependencies and Tools
######################
##############
## Helper tools
## software-properties-common needed for add-apt-repository
RUN apt-get update && \
    apt-get install -y git libxml2-dev r-cran-xml libcurl4-openssl-dev
#software-properties-common

# Install dependencies for r-base source build
RUN apt-get build-dep -y r-base

## Java needed for 3.2.1 R build
## Need toupdate this so you don;t have to click
#RUN add-apt-repository ppa:webupd8team/java && \
#    apt-get update && \
#    apt-get install oracle-java7-installer

## Install R 3.2.1
#RUN apt-get install -y r-base=${R_VERSION}

WORKDIR ${SRC}
# gfortran
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz && \
    apt-get install -y tcl tcl-dev tk-dev tk && \
    tar xf R-3.2.1.tar.gz

WORKDIR ${SRC}/R-3.2.1
RUN ./configure --prefix=/usr/local --enable-R-shlib --with-tcltk && \
    make && \
    make install

## Install R
#RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
#    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
#    gpg -a --export E084DAB9 | apt-key add - && \
#    apt-get update && \
#    apt-get install -y r-base-core=${R_VERSION} && \
#    apt-get install -y r-recommended=${R_VERSION} && \
#    apt-get install -y r-base=${R_VERSION} && \
#    apt-get install -y r-base-dev=${R_VERSION}

# options(repos='http://cran.rstudio.com/')
RUN echo "options(repos=\"http://cran.rstudio.com/\")" > ${SRC}/install_pkgs.R && \
    echo "install.packages(c(\"GMD\", \"ape\", \"RColorBrewer\", \"optparse\", \"logging\"), repos=\"http://cran.us.r-project.org\", dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    Rscript ${SRC}/install_pkgs.R

# Pull in inferCNV source
WORKDIR ${SRC}
RUN git clone https://github.com/broadinstitute/inferCNV
