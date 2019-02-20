#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)

args = parser$parse_args()

library(infercnv)
library(futile.logger)
library(tidyverse)


infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

gene_order = infercnv_obj@gene_order
gene_order = cbind(gene_order, gene=rownames(gene_order))

cnv_to_expr_vals = list()

expr.data <- infercnv_obj@expr.data

cnv_mean_sd = infercnv:::get_spike_dists(infercnv_obj@.hspike)

chrs = unique(infercnv_obj@gene_order$chr)

groups = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)

samples = names(groups)


for (sample in samples) {
    pdf_name = sprintf("%s-%s.cnv_expr_densities_each_chr.pdf", infercnv_obj_file, sub("[^A-Za-z0-9]", "_", sample, perl=TRUE))
    pdf(pdf_name)
    
    message(sprintf("plotting sample: %s", sample))

    sample_cells = groups[[ sample ]]
    
    sample_expr = expr.data[, sample_cells]

    for (chr in chrs) {
        chr_gene_idx = which(infercnv_obj@gene_order$chr == chr)

        sample_gene_expr = sample_expr[chr_gene_idx,]

        normal_gene_expr = expr.data[chr_gene_idx, unlist(infercnv_obj@reference_grouped_cell_indices)]

        df = rbind(data.frame(class='allnormal', vals=as.numeric(normal_gene_expr) ),
                   data.frame(class='sample', vals=as.numeric(sample_gene_expr)) )
        
        message(sprintf("plotting sample: %s, %s", sample, chr))

        p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + ggtitle(sprintf("%s, %s", sample, chr))

        p = p +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd))



        plot(p)

    }
    dev.off()
}

