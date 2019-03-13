#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
options(error = function() {traceback(2);quit(save = "no", status = 0, runLast = FALSE)})

parser = ArgumentParser()

parser$add_argument("--infercnv_obj", help="total sum normalized infercnv obj", required=TRUE, default=NULL, nargs=1)
parser$add_argument("--ncells", help="number of cells to simulate", required=FALSE, type='integer', nargs=1, default=-1)
parser$add_argument("--ngenes", help="number of genes to simulate", required=FALSE, type='integer', nargs=1, default=-1)
parser$add_argument("--output_prefix", help='prefix for output matrix file', required=TRUE, nargs=1)

args = parser$parse_args()

library(infercnv)
library(SingleCellExperiment)
library("methods")
library(tidyverse)


infercnv_obj_file = args$infercnv_obj

ncells = args$ncells
ngenes = args$ngenes
output_prefix = args$output_prefix

infercnv_obj = readRDS(infercnv_obj_file)

expr.data = infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)]

if (ncells < 0) {
    ncells = ncol(expr.data)
}
if (ngenes < 0) {
    ngenes = nrow(expr.data)
}

## sim using specified gene means
gene_means = rowMeans(expr.data)
gene_means = gene_means[gene_means>0]

gene_means = sample(x=gene_means, size=ngenes, replace=T)

newnames = paste0('gene', 1:ngenes)

names(gene_means) = newnames


sim_matrix <- infercnv:::.get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj, gene_means, ncells, TRUE)


output_filename = paste0(output_prefix, ".counts.matrix")
write.table(sim_matrix, file=output_filename, quote=F, sep='\t')

pdf(paste0(output_prefix, ".KS.pdf"))
infercnv:::KS_plot("meanVarSim", as.numeric(log(expr.data+1)), as.numeric(log(sim_matrix+1)))






