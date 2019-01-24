#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
options(error = function() {traceback(2);quit(save = "no", status = 0, runLast = FALSE)})

parser = ArgumentParser()
parser$add_argument("--counts_matrix", help="raw counts matrix file", required=TRUE, nargs=1)
parser$add_argument("--ncells", help="number of cells to simulate", required=TRUE, type='integer', nargs=1)
parser$add_argument("--ngenes", help="number of genes to simulate", required=TRUE, type='integer', nargs=1)
parser$add_argument("--output", help='name of output matrix file', required=TRUE, nargs=1)

args = parser$parse_args()

library(infercnv)
library(SingleCellExperiment)
library("methods") 
library(splatter)

data = read.table(args$counts_matrix)
ncells = args$ncells
ngenes = args$ngenes
output_filename = args$output

data = as.matrix(data)


write_sim_matrix <- function(counts_matrix, ncells, ngenes, output_filename) {

    params <- infercnv:::.estimateSingleCellParamsSplatterScrape(counts_matrix)
    
    #' normalize first:
    cs = colSums(counts_matrix)
    median_cs = median(cs)
    data <- sweep(counts_matrix, STATS=cs, MARGIN=2, FUN="/")
    data <- data * median_cs
    
    ## sim using specified gene means
    gene_means = rowMeans(data)
    gene_means = gene_means[gene_means>1]
        
    gene_means = sample(x=gene_means, size=ngenes, replace=T)
    
    newnames = paste0('gene', 1:ngenes)
    
    names(gene_means) = newnames
    
    
    params[['nGenes']] = ngenes
    params[['nCells']] = ncells

        
    sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
    sim_matrix <- counts(sim_matrix)
    
    write.table(sim_matrix, file=output_filename, quote=F, sep='\t')
    
}

write_sim_matrix(data, ncells, ngenes, output_filename)







