#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(futile.logger)
library(dplyr)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

ref_group_cell_indices = infercnv:::get_reference_grouped_cell_indices(infercnv_obj)
pdf_filename = paste0(infercnv_obj_file, ".chr_expr_densities.pdf")

cnv_mean_sd = infercnv:::get_spike_dists(infercnv_obj@.hspike)

pdf(pdf_filename)

chrs = unique(infercnv_obj@gene_order$chr)


for (chr in chrs) {
        
    gene_idx = which(infercnv_obj@gene_order$chr == chr)
    
    ref_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx,ref_group_cell_indices])
    
    df = data.frame(class='normal', vals=ref_data_pts)
    
    for (tumor in names(infercnv_obj@observation_grouped_cell_indices) ) {
        
        tumor_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ tumor ]]
        tumor_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx, tumor_cell_idx])
        
        df = rbind(df, data.frame(class=tumor, vals=tumor_data_pts))
    }

    flog.info(sprintf("Plotting data for chr: %s", chr))
    
    p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + scale_y_continuous(trans='log10', limits=c(1,NA)) + ggtitle(chr)

    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd)) 
    


    plot(p)
}

