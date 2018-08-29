#!/usr/bin/env Rscript


options(error = function() traceback(2))

C_CHR <- "chr"
C_START <- "start"
C_STOP <- "stop"
C_VIS_OUTLIER_CHOICES <- c("average_bound")
C_REF_SUBTRACT_METHODS <- c("by_mean", "by_quantiles")
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
C_OUTPUT_FORMAT <- c("pdf", "png")

