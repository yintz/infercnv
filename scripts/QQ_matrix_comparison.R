#!/usr/bin/env Rscript

set.seed(1234)

suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--matrix1", required=T, nargs=1)
parser$add_argument("--matrix2", required=T, nargs=1)
parser$add_argument("--log", required=F, default=FALSE, action="store_true")
parser$add_argument("--output", required=T, nargs=1, help="output filename png")

args = parser$parse_args()


#' learn distribution parameters:
data1 = as.matrix(read.table(args$matrix1, header=T, row.names=1))
data2 = as.matrix(read.table(args$matrix2, header=T, row.names=1))


png(args$output)
if (args$log) {
    data1 = log(data1+1)
    data2 = log(data2+1)
}
qqplot(data1, data2)
abline(a=0,b=1, col='red')



