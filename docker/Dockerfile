# Docker file for inferCNV
FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
gpg -a --export E084DAB9 | apt-key add -

RUN apt-get update && \
apt-get -y install curl libssl-dev libcurl4-openssl-dev libxml2-dev r-base r-base-dev git

RUN echo "options(repos = c(CRAN = \"https://cran.rstudio.com\"))" >.Rprofile && \
echo "install.packages(c(\"devtools\",\"ape\",\"GMD\",\"RColorBrewer\",\"optparse\",\"logging\"), dependencies =TRUE)" > install_devtools.r && \
echo "library(\"devtools\")" >> install_devtools.r && \
echo "install_github(\"broadinstitute/inferCNV\")" >> install_devtools.r && \
R --no-save < install_devtools.r

RUN apt-get update && apt-get install -y python3

RUN git clone https://github.com/broadinstitute/inferCNV.git
ENV PATH=${PATH}:/inferCNV/scripts

CMD inferCNV.R --help
