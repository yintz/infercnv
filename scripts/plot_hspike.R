#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(ggplot2)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

if (! is.null(infercnv_obj@.hspike)) {
    out_prefix = paste0(infercnv_obj_file, '.hspike')
    plot_cnv(infercnv_obj@.hspike,
             output_filename=basename(out_prefix))

    hspike_obj = infercnv_obj@.hspike
    hspike_gene_expr_by_cnv <- infercnv:::.get_gene_expr_by_cnv(hspike_obj)
    hspike_cnv_mean_sd <- infercnv:::.get_gene_expr_mean_sd_by_cnv(hspike_gene_expr_by_cnv)
    p = infercnv:::.plot_gene_expr_by_cnv(hspike_gene_expr_by_cnv, hspike_cnv_mean_sd)
    pdf(paste0(infercnv_obj_file, '.hspike.dist.pdf'))
    plot(p)
    dev.off()
    
} else {
    message("no hspike to plot")
}



