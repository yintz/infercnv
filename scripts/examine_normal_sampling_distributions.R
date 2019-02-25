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

if (! infercnv:::has_reference_cells(infercnv_obj)) {
    stop("Error, cannot tune parameters without reference 'normal' cells defined")
}

expr_vals <- infercnv_obj@expr.data
mu = mean(expr_vals)
sigma = sd(expr_vals)
nrounds = 1000
sds = c()
ngenes = nrow(expr_vals)

normal_samples = infercnv_obj@reference_grouped_cell_indices

num_normal_samples = length(normal_samples)

mean_vals_df = NULL;
z_p_val = 0.05

num_cells_to_empirical_sd = list()

ncells_partitions = seq (1,100,5)
for (ncells in ncells_partitions) {
    means = c()
    
    message(sprintf("num cells: %g", ncells))

    cells_counted = 0;
    
    for(i in 1:nrounds) {
        ## pick a random gene
        rand.gene = sample(1:ngenes)
        
        ## pick a random normal cell type
        rand.sample = sample(num_normal_samples)
        #rand.sample=1
        
        vals = sample(expr_vals[rand.gene, normal_samples[[rand.sample]] ], size=ncells, replace=T)
        m_val = mean(vals)
        means = c(means,  m_val)

        cells_counted = cells_counted + length(vals)

        
    }
    my.sd = sd(means)
    sds = c(sds, my.sd)

    num_cells_to_empirical_sd[[ ncells ]] = my.sd
    
    df = data.frame(num_cells=ncells, vals=means)
                                        #print(df)
    if(is.null(mean_vals_df)) {
        mean_vals_df = df
    } else {
        mean_vals_df = rbind(mean_vals_df, df)
    }
    
}

## fit linear model
num_cells = ncells_partitions

write.table(data.frame(num_cells=num_cells, sds=sds), file='num_cells_vs_sds.table.dat', quote=F, sep="\t")


fit = lm(log(sds) ~ log(num_cells)) #note, hbadger does something similar, but not for the hmm cnv state levels

my.spline = smooth.spline(log(num_cells), log(sds)) 

message("plotting log(sd) vs. log(num_cells)")

plot(log(num_cells), log(sds), main='log(sd) vs. log(num_cells)')

plot(num_cells, sds, main='sd vs. num_cells')

my.spline2 = smooth.spline(num_cells, sds) 

## store mean_delta for the single gene for convenience sake
mean_delta = qnorm(p=1-z_p_val, sd=sigma, mean=0)

normal_sd_trend = list(mu=mu,
                       sigma=sigma,
                       fit=fit,
                       spline=my.spline,
                       mean_delta=mean_delta)



### do some plotting


for (ncells in ncells_partitions) {

    message(sprintf("plotting ncells distribution: %g", ncells))
    
    data.want = mean_vals_df %>% filter(num_cells == ncells)
    
    
    p = data.want %>% ggplot(aes(vals, fill=num_cells)) +
        geom_density(alpha=0.3)

    sigma <- exp(predict(normal_sd_trend$fit,
                         newdata=data.frame(num_cells=ncells))[[1]]) 

    message("ncells:", ncells, " sigma: ", sigma)
    
    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=1,'sd'=sigma))  +
        ggtitle(sprintf("num_cells: %g, sd: %g", ncells, sigma))

    p = p +
        stat_function(fun=dnorm, color='magenta', args=list('mean'=1,'sd'=num_cells_to_empirical_sd[[ ncells]] )) 


    pval=0.01
    
    left_mean = 1 - 2 * (1-qnorm(p=pval, mean=1, sd=sigma))
    message("left_mean: ", left_mean)
    p = p +
        stat_function(fun=dnorm, color='blue', args=list('mean'=left_mean,'sd'=sigma))


    right_mean = 1 + 2 * (qnorm(p=1-pval, mean=1, sd=sigma)-1)
    message("right_mean: ", right_mean)        
        p = p +
            stat_function(fun=dnorm, color='blue', args=list('mean'=right_mean,'sd'=sigma))  
    
    

    

    if (FALSE) {
    
        spline.sd = exp(predict(my.spline, x=log(ncells))$y)
        
        
        p = p +
            stat_function(fun=dnorm, color='green', args=list('mean'=1,'sd'=spline.sd)) 
        
        spline2.sd = predict(my.spline2, x=ncells)$y
        
        message(spline2.sd)
        
        p = p +
            stat_function(fun=dnorm, color='orange', args=list('mean'=1,'sd'=spline2.sd)) 
    }
    
    plot(p)
}

