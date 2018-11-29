#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(dplyr)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

infercnv_obj = remove_spike(infercnv_obj)

normal_idx = unlist(infercnv_obj@reference_grouped_cell_indices)
normal_expr_mtx = infercnv_obj@expr.data[,normal_idx]

simnormal_idx = unlist(unlist(infercnv_obj@observation_grouped_cell_indices))
simnormal_expr_mtx = infercnv_obj@expr.data[,simnormal_idx]

pdf(paste0(infercnv_obj_file, ".qq_plots.ALL.pdf"))

#par(mfrow=c(3,1))


plot_dists <- function(normal_mtx, simnormal_mtx, title) {

    df = data.frame(class='normal', vals=as.numeric(normal_mtx))
    df = rbind(df, data.frame(class='simnormal', vals=as.numeric(simnormal_mtx)))

    p = df %>% filter(vals>0) %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + scale_y_continuous(limits=c(0,NA)) + ggtitle(title)
    plot(p)

}

get_pct_zeros <- function(mtx) {
    vals = as.numeric(mtx)
    n_vals = length(vals)
    n_zero = sum(vals==0)
    pct_zero = n_zero/n_vals * 100
    return(pct_zero)
}

plot_pct_zeros <- function(normal_mtx, simnormal_mtx) {
    normal_pct_zeros = get_pct_zeros(normal_mtx)
    simnormal_pct_zeros = get_pct_zeros(simnormal_mtx)

    barplot(c(normal_pct_zeros, simnormal_pct_zeros), names.arg=c('norm',  'simnorm'), main='pct zeros', ylim=c(0,100))
}

plot_dists(normal_expr_mtx, simnormal_expr_mtx, 'all')
qqplot(normal_expr_mtx, simnormal_expr_mtx, main="normal vs. simnormal")
abline(a=0,b=1, col='red')

plot_pct_zeros(normal_expr_mtx, simnormal_expr_mtx)

dev.off()

# do per chr
chrs = unique(infercnv_obj@gene_order$chr)
for (chr in chrs) {

    gene_idx = which(infercnv_obj@gene_order$chr == chr)

    message(sprintf("Processing chr: %s, w/ %d  genes", chr, length(gene_idx)))
    
    chr_normal_mtx = normal_expr_mtx[gene_idx, ]

    chr_simnormal_mtx = simnormal_expr_mtx[gene_idx, ]
    
    pdf(paste0(infercnv_obj_file, sprintf(".qq_plots.%s.pdf", chr)))
    qqplot(chr_normal_mtx, chr_simnormal_mtx, main=chr)
    abline(a=0,b=1, col='red')
    
    plot_dists(chr_normal_mtx, chr_simnormal_mtx, chr)

    plot_pct_zeros(chr_normal_mtx, chr_simnormal_mtx)

    dev.off()
}



