#!/usr/bin/env Rscript


hclust_method='ward.D2'

num_rand_iters = 100
MAX_PVAL=0.05

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(futile.logger)
library(pheatmap)

obj = readRDS(args$infercnv_obj)

tumor.expr.data = obj@expr.data[, unlist(obj@observation_grouped_cell_indices)]

gene_order = obj@gene_order
chrs = unique(gene_order$chr)


pdf("test.recursive_trees.pdf")


ALL_CLUSTERS = list()
MIN_CLUSTER_SIZE=3

library(sigclust2)

recursive_cluster_cutting <- function(expr.matrix) {

    message("recursive_cluster_cutting()")
    print(dim(expr.matrix))

    if (dim(expr.matrix)[2] < MIN_CLUSTER_SIZE) {
        message("cluster size too small. Storing cluster")
        ALL_CLUSTERS[[length(ALL_CLUSTERS)+1]] <<- colnames(expr.matrix)

        print("Returning")
        return(NULL)
        print("Didn't actually return...")
    }

    print("Onward")
    print(dim(expr.matrix))
    
    t_tumor.expr.data = t(expr.matrix) # cells as rows, genes as cols

    shc_result = shc(t_tumor.expr.data, metric='euclidean', linkage='ward.D2')
    plot(shc_result)
    
    for(chr in chrs) {
        chr_genes = which(gene_order$chr == chr)
        
        message(sprintf("plotting %s", chr))
        
        shc_result = shc(t_tumor.expr.data[,chr_genes], metric='euclidean', linkage='ward.D2')
        plot(shc_result)
    }
    
    
            
}

recursive_cluster_cutting(tumor.expr.data)

dev.off()

print(ALL_CLUSTERS)








