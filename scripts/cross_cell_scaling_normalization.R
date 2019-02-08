#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--log", help="log transform expr", action='store_true', default=FALSE)

args = parser$parse_args()

library(infercnv)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

expr.data = infercnv_obj@expr.data


## do upper quartile normalization
upper_quart = apply(expr.data, 2, quantile, probs=0.75)
mean_upper_quart = mean(upper_quart)
revised.expr.data = sweep(expr.data, 2, mean_upper_quart/upper_quart, "*")

new_upper_quart = apply(revised.expr.data, 2, quantile, probs=0.75) 

print(new_upper_quart)

infercnv_obj@expr.data = revised.expr.data

saveRDS(infercnv_obj, 'rescaled.obj')

