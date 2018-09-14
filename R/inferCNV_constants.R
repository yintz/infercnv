#!/usr/bin/env Rscript


options(error = function() { traceback(2); stop("Error encountered") } )

C_CHR <- "chr"
C_START <- "start"
C_STOP <- "stop"
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
C_OUTPUT_FORMAT <- c("pdf", "png")

