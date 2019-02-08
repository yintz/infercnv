#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--window_size", help="window size", required=FALSE, type='integer', default=11)
args = parser$parse_args()

library(infercnv)
library(ggplot2)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

infercnv_obj = infercnv:::.subcluster_tumors_general(infercnv_obj)

mf_infercnv_obj = infercnv:::.apply_heatmap_median_filtering(infercnv_obj, window_size=args$window_size)

saveRDS(mf_infercnv_obj, file=sprintf("%s-median_filtered.W%d.obj", infercnv_obj_file, args$window_size) )

plot_cnv(mf_infercnv_obj, output_filename=paste0(infercnv_obj_file, sprintf(".mf.W%d", args$window_size)))


