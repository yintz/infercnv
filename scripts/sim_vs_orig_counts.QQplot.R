#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--counts_matrix", help="raw counts matrix file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)

data = read.table(args$counts_matrix)
data = as.matrix(data)

params <- infercnv:::.estimateSingleCellParamsSplatterScrape(data)


#' normalize first:
cs = colSums(data)
median_cs = median(cs)
data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
data <- data * median_cs

## get sim params
sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params)    
sim_matrix <- counts(sim_matrix)

pdf("sim_vs_orig_counts.qqplots.pdf")
qqplot(log(as.numeric(data)+1), log(as.numeric(sim_matrix)+1), main='orig vs. full sim')
abline(a=0,b=1,col='red')


## sim using specified gene means
gene_means = rowMeans(data)
gene_means = gene_means[gene_means>0]

params[['nGenes']] = length(gene_means)
sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
sim_matrix <- counts(sim_matrix)

qqplot(log(as.numeric(data)+1), log(as.numeric(sim_matrix)+1), main='orig vs. sim w/ means')
abline(a=0,b=1,col='red')


