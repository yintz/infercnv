#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--counts_matrix", help="raw counts matrix file", required=TRUE, nargs=1)
parser$add_argument("--cnv", help='cnv level', required=TRUE, type='double', nargs=1)

args = parser$parse_args()

library(infercnv)
library(tidyverse)

data = read.table(args$counts_matrix)
data = as.matrix(data)

params <- infercnv:::.estimateSingleCellParamsSplatterScrape(data)

params[['nCells']] = 100

#' normalize first:
cs = colSums(data)
median_cs = median(cs)
data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
data <- data * median_cs

## sim using specified gene means
gene_means = rowMeans(data)
gene_means = gene_means[gene_means>0]


params[['nGenes']] = length(gene_means)
sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
sim_matrix <- counts(sim_matrix)


num_genes_change=200
fold_change = args$cnv

pdf(sprintf("sim_cnv-%g-delta.pdf", fold_change))

orig_gene_means = gene_means
gene_means[1:num_genes_change] = gene_means[1:num_genes_change] * fold_change
alt_matrix =  infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
alt_matrix = counts(alt_matrix)

# normalize each to total sum counts
target_median = median(c(colSums(sim_matrix), colSums(alt_matrix)))

normalize_counts <- function(matrix, target_sum) {
    matrix = sweep(matrix, 2, colSums(matrix), "/")
    matrix = matrix * target_sum
    return(matrix)
}   

sim_matrix = normalize_counts(sim_matrix, target_median)
alt_matrix = normalize_counts(alt_matrix, target_median)



## restrict to changed subset
sim_matrix = sim_matrix[1:num_genes_change,]
alt_matrix = alt_matrix[1:num_genes_change,]


alt_genemeans = rowMeans(alt_matrix)



plot(orig_gene_means[1:num_genes_change], alt_genemeans, pch='.', log='xy', main='orig_target_gene_means vs. sim w/ cnv')
abline(a=0,b=1,col='red')


plot(rowMeans(sim_matrix), alt_genemeans, pch='.', log='xy', main='mean(gene) sim w/o vs. sim w/ cnv')
abline(a=0,b=1,col='red')


density_plot <- function(title, normal_expr, tumor_expr) {

    tumor_expr <- as.numeric(tumor_expr)
    hspike_expr <- as.numeric(normal_expr)
    
    df = data.frame(class='tumor', expr=as.numeric(tumor_expr))
    df = rbind(df, data.frame(class='norm', expr=as.numeric(normal_expr)))

    p = df %>% ggplot(aes(expr, fill=class)) + geom_density(alpha=0.3) + ggtitle(title)
    
    plot(p)
}


sim_matrix.log = log(sim_matrix+1)
alt_matrix.log = log(alt_matrix+1)

density_plot(sprintf("cnv: %g", fold_change), sim_matrix.log, alt_matrix.log)


sim_matrix.log.means = rowMeans(sim_matrix.log)
sim_matrix.log = sweep(sim_matrix.log, 1, sim_matrix.log.means, '-')
alt_matrix.log = sweep(alt_matrix.log, 1, sim_matrix.log.means, '-')

density_plot(sprintf("delta cnv: %g", fold_change), sim_matrix.log, alt_matrix.log)

smoothScatter(sim_matrix.log.means, rowMeans(alt_matrix.log), main='gene means, normal subtracted')
abline(h=log(fold_change), col='magenta')

sim_matrix.log.smoothed = apply(sim_matrix.log, 1, caTools::runmean, k=101)
alt_matrix.log.smoothed = apply(alt_matrix.log, 1, caTools::runmean, k=101)

density_plot(sprintf("cnv: %g, smoothed", fold_change), sim_matrix.log.smoothed, alt_matrix.log.smoothed)

all_alt_smoothed_vals = c(as.numeric(alt_matrix.log.smoothed))

## smoothing results depends on gene chr order.
for (round in seq(5)) {
    alt_matrix.log.reordered = alt_matrix.log[ sample(seq(ncol(alt_matrix.log)), replace=T), ]
    alt_matrix.log.reordered.smoothed = apply(alt_matrix.log.reordered, 1, caTools::runmean, k=101)
    density_plot(sprintf("smoothR: %g, cnv: %g, smoothed", round, fold_change), sim_matrix.log.smoothed, alt_matrix.log.reordered.smoothed)
    
    all_alt_smoothed_vals = c(all_alt_smoothed_vals, as.numeric(alt_matrix.log.reordered.smoothed))
}

density_plot(sprintf("all smoothR cnv: %g, smoothed", fold_change), sim_matrix.log.smoothed, all_alt_smoothed_vals)
