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
    pdfname = paste0(infercnv_obj_file, '.hspike.diff_normal_tumor.pdf')

    pdf(pdfname)
    hspike = infercnv_obj@.hspike

    normal_matrix = hspike@expr.data[,unlist(hspike@reference_grouped_cell_indices)]
    tumor_matrix = hspike@expr.data[,unlist(hspike@observation_grouped_cell_indices)]

    normal.means = rowMeans(normal_matrix)
    tumor.means = rowMeans(tumor_matrix)

    plot(normal.means, ylim=range(normal.means, tumor.means))
    points(tumor.means, col='green')

    plot(tumor.means - normal.means)
    abline(h=0, col='red')
    
    sm = caTools::runmean(tumor.means - normal.means, k=31)
    points(sm, col='magenta')
    
} else {
    message("no hspike to plot")
}



