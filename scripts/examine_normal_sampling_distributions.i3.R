#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(tidyverse)
library(futile.logger)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)


expr_vals <- infercnv_obj@expr.data


sd_trend_info = infercnv:::.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj)


mu = sd_trend_info$mu
sigma = sd_trend_info$sigma

sds = c()
ngenes = nrow(expr_vals)

tumor_samples = infercnv_obj@observation_grouped_cell_indices

print(tumor_samples)

num_tumor_samples = length(tumor_samples)

print(num_tumor_samples)

mean_vals_df = NULL;
z_p_val = 0.05


num_cells_to_empirical_sd = list()

nrounds=100

ncells_partitions = seq (1,100,5)
for (ncells in ncells_partitions) {
    means = c()
    
    message(sprintf("num cells: %g", ncells))

    cells_counted = 0;
    
    for(i in 1:nrounds) {
        ## pick a random gene
        rand.gene = sample(1:ngenes, size=1)
        
        ## pick a random normal cell type
        rand.sample = sample(1:num_tumor_samples, size=1)
                                        #rand.sample=1
        #print(rand.sample)
        
        vals = sample(expr_vals[rand.gene, tumor_samples[[rand.sample]] ], size=ncells, replace=T)
        m_val = mean(vals)
        means = c(means,  m_val)
        
        cells_counted = cells_counted + length(vals)
                
    }
    means.sd = sd(means)
    means.mean = mean(means)
    
    num_cells_to_empirical_sd[[ ncells ]] = means.sd
    
    df = data.frame(num_cells=ncells, vals=means)

    message(sprintf("plotting ncells distribution: %g", ncells))
    
    data.want = df
    
    
    p = data.want %>% ggplot(aes(vals, fill=num_cells)) +
        geom_density(alpha=0.3)   +
        ggtitle(sprintf("num_cells: %g", ncells))

    ## draw parameterized distribution
    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=means.mean,'sd'=means.sd))
    

    alpha=0.05
    ks_delta = infercnv:::get_HoneyBADGER_setGexpDev(gexp.sd=sd_trend_info$sigma, k_cells=ncells, alpha=alpha, plot=T)
    
    left_mean = means.mean - ks_delta
    message("left_mean: ", left_mean)
    p = p +
        stat_function(fun=dnorm, color='blue', args=list('mean'=left_mean,'sd'=means.sd))


    right_mean = means.mean + ks_delta
    message("right_mean: ", right_mean)        
        p = p +
            stat_function(fun=dnorm, color='blue', args=list('mean'=right_mean,'sd'=means.sd))  
    
        
    plot(p)
}

