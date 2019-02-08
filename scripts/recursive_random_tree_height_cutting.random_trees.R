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

infercnv_obj = readRDS(args$infercnv_obj)


pdf("test.recursive_trees.pdf")

adj.obj = infercnv:::define_signif_tumor_subclusters(infercnv_obj, p_val=0.05, hclust_method='ward.D2', partition_method='random_trees')





