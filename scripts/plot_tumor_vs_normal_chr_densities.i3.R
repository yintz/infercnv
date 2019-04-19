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
pdf_filename = paste0(infercnv_obj_file, ".i3.chr_expr_densities.pdf")

normal_sd_trend = infercnv:::.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj)

mu = normal_sd_trend$mu
sigma = normal_sd_trend$sigma



pdf(pdf_filename)

chrs = unique(infercnv_obj@gene_order$chr)

delta = infercnv:::get_HoneyBADGER_setGexpDev(gexp.sd=sigma, alpha=0.05, k_cells=7)

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
    
    p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + ggtitle(chr) # + scale_y_continuous(trans='log10', limits=c(1,NA))
    
    
    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=mu,'sd'=sigma)) +
        stat_function(fun=dnorm, color='blue', args=list('mean'=mu-delta,'sd'=sigma)) +
        stat_function(fun=dnorm, color='blue', args=list('mean'=mu+delta,'sd'=sigma)) 
    

    plot(p)
}

