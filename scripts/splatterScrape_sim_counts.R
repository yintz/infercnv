#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
options(error = function() {traceback(2);quit(save = "no", status = 0, runLast = FALSE)})

parser = ArgumentParser()

parser$add_argument("--counts_matrix", help="raw counts matrix file", required=TRUE, default=NULL, nargs=1)
parser$add_argument("--ncells", help="number of cells to simulate", required=TRUE, type='integer', nargs=1)
parser$add_argument("--ngenes", help="number of genes to simulate", required=TRUE, type='integer', nargs=1)
parser$add_argument("--output", help='name of output matrix file', required=TRUE, nargs=1)

args = parser$parse_args()

library(infercnv)
library(SingleCellExperiment)
library("methods") 
library(splatter)


counts_matrix = read.table(args$counts_matrix)
params_file = sprintf("%s.params_obj", args$counts_matrix)
if (file.exists(params_file)) {
    message("-note, reusing stored params")
    params = readRDS(params_file)
} else {
    params <- infercnv:::.estimateSingleCellParamsSplatterScrape(counts_matrix)
    saveRDS(params, file=sprintf("%s.params_obj", args$counts_matrix))
}

ncells = args$ncells
ngenes = args$ngenes
output_filename = args$output

data = as.matrix(counts_matrix)

#' normalize first:
cs = colSums(counts_matrix)
median_cs = median(cs)
data <- sweep(counts_matrix, STATS=cs, MARGIN=2, FUN="/")
data <- data * median_cs

## sim using specified gene means
gene_means = rowMeans(data)
gene_means = gene_means[gene_means>0]

gene_means = sample(x=gene_means, size=ngenes, replace=T)

newnames = paste0('gene', 1:ngenes)

names(gene_means) = newnames


params[['nGenes']] = ngenes
params[['nCells']] = ncells


sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
sim_matrix <- counts(sim_matrix)

write.table(sim_matrix, file=output_filename, quote=F, sep='\t')






