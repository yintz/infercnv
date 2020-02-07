# Docker file for inferCNV
FROM bioconductor/bioconductor_docker:devel

LABEL org.label-schema.license="BSD-3-Clause" \
      org.label-schema.vendor="Broad Institute" \
      maintainer="Christophe Georgescu <cgeorges@broadinstitute.org>"

RUN apt-get update && apt-get -y install curl libssl-dev libcurl4-openssl-dev \
                                        libxml2-dev git python3 jags \
                                        r-cran-rjags time && \
                      apt-get clean && rm -rf /var/tmp/* \
                                          /tmp/* /var/lib/apt/lists/*

# Install R and Bioconductor packages
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com'))" >.Rprofile
RUN R -e "install.packages(c('devtools','KernSmooth', 'lattice', 'Matrix', \
                             'survival', 'MASS', 'TH.data', 'nlme', 'ape', \
                             'fitdistrplus', 'multcomp', 'coin', 'binhf', \
                             'caTools', 'coda', 'dplyr', 'doparallel', \
                             'fastcluster', 'foreach', 'futile.logger', \
                             'future', 'gplots', 'ggplot2', 'HiddenMarkov', \
                             'reshape', 'rjags', 'RColorBrew', 'doParallel', \
                             'tidyr', 'gridExtra', 'argparse', 'knitr', \
                             'rmarkdown', 'testthat', 'optparse', 'logging', \
                             'data.table', 'BiocManager'), repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages(c('cluster', 'Seurat'), repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install(c('BiocGenerics', 'edgeR', 'SingleCellExperiment', \
            'SummarizedExperiment', 'BiocStyle', 'BiocCheck'))"

# Checkout and install infercnv
# update to version bump commit
RUN git clone https://github.com/broadinstitute/infercnv && cd infercnv && \
      git checkout master && git checkout cc19b7efaf09e3eff329a93df1ada475c73d2b40 && \
      R CMD INSTALL . 

ENV PATH=${PATH}:/infercnv/scripts

CMD inferCNV.R --help

