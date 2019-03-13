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

tumor.expr.data[tumor.expr.data>3] <- 4
tumor.expr.data[tumor.expr.data<3] <- 2


pdf("test.recursive_trees.pdf")


ALL_CLUSTERS = list()
MIN_CLUSTER_SIZE=3


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
    d = dist(t_tumor.expr.data)

    h_obs = hclust(d, method=hclust_method)

    # permute by chromosomes
    
    permute_chr_col_vals <- function(df) {

        num_cells = nrow(df)

        for(chr in chrs) {
            chr_genes = which(gene_order$chr == chr)

            df[, chr_genes] = df[sample(x=1:num_cells, size=num_cells, replace=F), chr_genes]
        }

        df
    }

    permute_col_vals <- function(df) {

        num_cells = nrow(df)
        for (i in 1:ncol(df)) {
            df[,i] = df[sample(x=1:num_cells, size=num_cells, replace=F), i]
        }
        
        df
    }

    
    example_rand_matrix <- NULL
    max_rand_heights = c()
    for (i in 1:num_rand_iters) {
        
        ##rand.tumor.expr.data = permute_chr_col_vals(t_tumor.expr.data)
        rand.tumor.expr.data = permute_col_vals(t_tumor.expr.data)
        example_rand_matrix <- rand.tumor.expr.data
        rand.dist = dist(rand.tumor.expr.data)
        h_rand <- hclust(rand.dist, method=hclust_method)

        max_rand_heights = c(max_rand_heights, max(h_rand$height))
    }
        
    h = h_obs$height

    max_height = max(h)
    
    message(sprintf("Max Rand Heights(h): %s", paste(max_rand_heights, sep=",", collapse=",")))
    
    max_rand_height_dens = density(max_rand_heights)
    plot(max_rand_height_dens, xlim=range(max_rand_height_dens$x, max_height))
    
    e = ecdf(max_rand_heights)
    message(sprintf("pvals(Lengths(h)): %s", paste(1-e(h), sep=",", collapse=",")))
    
    pval = 1- e(max_height)
    message(sprintf("pval for max obs height: %g = %g", max_height, pval))
    
    abline(v=max_height, col='red')
    
    pheatmap(t(expr.matrix), cluster_cols=F)
    pheatmap(example_rand_matrix, cluster_cols=F)


    #stop("stopping")
    
    if (max_height > 0 & pval <= MAX_PVAL) {
        ## keep on cutting.
        cut_height = mean(c(h[length(h)-1], h[length(h)]))
        message(sprintf("cutting at height: %g",  cut_height))
        grps = cutree(h_obs, h=cut_height)
        print(grps)
        uniqgrps = unique(grps)
        for (grp in uniqgrps) {
            grp_idx = which(grps==grp)
            
            message(sprintf("grp: %s  contains idx: %s", grp, paste(grp_idx,sep=",", collapse=","))) 
            df = expr.matrix[,grp_idx,drop=F]
            recursive_cluster_cutting(df)
        }
    } else {
        message("No cluster pruning")
        ALL_CLUSTERS[[length(ALL_CLUSTERS)+1]] <<- colnames(expr.matrix)
    }
        
}

recursive_cluster_cutting(tumor.expr.data)

dev.off()

print(ALL_CLUSTERS)








