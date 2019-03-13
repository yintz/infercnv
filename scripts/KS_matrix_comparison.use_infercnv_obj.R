#!/usr/bin/env Rscript

set.seed(1234)

suppressPackageStartupMessages(library("argparse"))

library(tidyverse)

parser = ArgumentParser()
parser$add_argument("--matrix1", required=T, nargs=1)
parser$add_argument("--infercnv_obj", required=T, nargs=1)
parser$add_argument("--log", required=F, default=FALSE, action="store_true")
parser$add_argument("--output", required=T, nargs=1, help="output filename pdf")

args = parser$parse_args()


#' learn distribution parameters:
data1 = as.matrix(read.table(args$matrix1, header=T, row.names=1))



infercnv_obj_file = args$infercnv_obj
infercnv_obj = readRDS(infercnv_obj_file)
data2 = as.matrix(infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)])


png(args$output)
if (args$log) {
    data1 = log(data1+1)
    data2 = log(data2+1)
}


## plotting ideas borrowed from
## https://stackoverflow.com/questions/39162178/kolmogorov-smirnov-plot-in-r-ggplot


m1_ecdf = ecdf(data1)
m2_ecdf = ecdf(data2)
val_range = range(data1, data2)
step = (val_range[2] - val_range[1])/100
vals = seq(val_range[1], val_range[2], step)


m1_cdf = m1_ecdf(vals)
m2_cdf = m2_ecdf(vals)

cdfs = data.frame(vals,
                  m1_cdf,
                  m2_cdf)

ks_point = which.max(abs(cdfs$m1_cdf - cdfs$m2_cdf))
ks_point_info = cdfs[ks_point,]
##message("KS point info: ", paste(ks_point_info, collapse=', '))

cdfs = cdfs %>% gather('m1_cdf', 'm2_cdf', key='type', value='cdf')


ggplot(cdfs, aes(x=vals, y=cdf)) +
    geom_line(aes(color=type, linetype=type)) +
    geom_segment(aes(x=ks_point_info$vals,
                     y=ks_point_info$m1_cdf,
                     xend=ks_point_info$vals,
                     yend=ks_point_info$m2_cdf), color='magenta', size=2) +
    ggtitle(sprintf("%s vs. %s KS", args$matrix1, args$matrix2)) + xlab("number") + ylab("cdf")





