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

pdf('ladeda.pdf')

normal_groups = infercnv_obj@reference_grouped_cell_indices
tumor_groups = infercnv_obj@observation_grouped_cell_indices

expr.data = infercnv_obj@expr.data

num_tumor_groups = length(tumor_groups)

par(mfrow=c(num_tumor_groups, 1))

library(tidyverse)


plotme <- function(normal_pts, tumor_pts, windowsize) {

    all_pts = c(normal_pts, tumor_pts)

    all_pts_names = names(all_pts)

    my.colors = rainbow(length(all_pts))
    
    yrange = range(unlist(all_pts))

    text.adj = 0.7
    for (i in 1:length(all_pts)) {
        if (i == 1) {
            plot(all_pts[[i]], t='l', col=my.colors[i], main=sprintf("windowsize: %g, tumor: %s", windowsize, all_pts_names[length(all_pts_names)]), ylim=yrange,
                 cex.lab=text.adj, cex.main=text.adj, cex.axis=text.adj)
        } else {
            points(all_pts[[i]], t='l', col=my.colors[i])
        }
    }
    abline(h=0)
    legend('top', legend=all_pts_names, col=my.colors, pch=1, horiz=T, bty='n', cex=text.adj)
    
}



get_smoothed <- function(cell_idx, windowsize) {
    group_expr_data = expr.data[, cell_idx]
    smoothed = apply(group_expr_data, 2, caTools::runmean, k=windowsize)
    smoothed_mean = rowMeans(smoothed)

    ## center it:
    smoothed_mean = smoothed_mean - median(smoothed_mean)
    
    return(smoothed_mean)
}

plot_chr_smooths <- function(tumor_type) {

    
    tumor_pts = tumor_groups[[tumor_type]]
    
    windowsizes = c(25,50,75,100)
    for (windowsize in windowsizes) {
        message(sprintf("\t-plotting %s", tumor_type))

        normal_pts = list()
        for (normal_type in names(normal_groups)) {
            normal_pts[[ normal_type ]] <- get_smoothed(normal_groups[[normal_type]], windowsize)
        }
        
        tumor_pts = list()
        tumor_pts[[ tumor_type ]] = get_smoothed(tumor_groups[[tumor_type]], windowsize)
        plotme(normal_pts, tumor_pts, windowsize)
    }
}




for (tumor_type in names(tumor_groups)) {
    message(sprintf("plotting for %s", tumor_type))
    plot_chr_smooths(tumor_type)
}
