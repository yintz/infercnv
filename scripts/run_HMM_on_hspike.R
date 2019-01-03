#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(futile.logger)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

hspike = infercnv_obj@.hspike

hspike.hmm = infercnv:::predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj=hspike,
                                                                 cnv_mean_sd=infercnv:::get_spike_dists(hspike),
                                                                 cnv_level_to_mean_sd_fit=infercnv:::get_hspike_cnv_mean_sd_trend_by_num_cells_fit(hspike),
                                                                 p_val=0.05,
                                                                 hclust_method='ward.D2'
                                                                 )

plot_cnv(hspike.hmm, x.center=3, x.range=c(1,5), output_filename="hspike.hmm")

saveRDS(hspike.hmm, file=sprintf("%s-HMM.obj", infercnv_obj_file))

