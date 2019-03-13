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

if (! is.null(infercnv_obj@.hspike)) {
    hspike_obj = infercnv_obj@.hspike


    pdf(paste0(infercnv_obj_file, '.hspike.dist_by_numcells.pdf'))



    
    gene_expr_by_cnv <- infercnv:::.get_gene_expr_by_cnv(hspike_obj)
    cnv_level_to_mean_sd = list()

    for (ncells in c(1,2,3,4,5,10,20,50,100)) {
        
        cnv_to_means = list()
        cnv_mean_sd = list()
        
        for (cnv_level in names(gene_expr_by_cnv) ) {
            expr_vals = gene_expr_by_cnv[[ cnv_level ]]
            nrounds = 100
            
            means = c()
                        
            for(i in 1:nrounds) {
                vals = sample(expr_vals, size=ncells, replace=T)
                m_val = mean(vals)
                means = c(means,  m_val)
            }
            cnv_to_means[[ cnv_level ]] = means
            cnv_mean_sd[[ cnv_level ]] = list(sd=sd(means), mean=mean(means))
        }
        
        ## plot

        df  = do.call(rbind, lapply(names(cnv_to_means), function(x) { data.frame(cnv=x, expr=cnv_to_means[[x]]) }))
        
        p = df %>% ggplot(aes(expr,  fill=cnv, colour=cnv))  +  geom_density(alpha=0.1)
        
        p = p +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd)) +
            stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd)) 
        
        p = p + ggtitle(sprintf("num cells: %g", ncells))
        
        plot(p)

        
    }
    
        
    dev.off()
    
} else {
    message("no hspike to plot")
}



