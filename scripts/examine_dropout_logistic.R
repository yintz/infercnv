#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args) == 0) {
    stop("Error, require params: infercnv.obj");
}

infercnv_obj_file = args[1]

pdf(paste0(infercnv_obj_file, '.dropout.pdf'))

infercnv_obj = readRDS(infercnv_obj_file)


library(edgeR)
library(fitdistrplus)
library(infercnv)

# borrowing some code from splatter

get_parameters <- function(group_name, expr.matrix) {

    params = list()
    params[['group_name']] = group_name
    
    # estimate gamma for  genes
    lib.sizes <- colSums(expr.matrix)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(expr.matrix) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    
    # estimate dropout params
    mean_vs_p0_table = infercnv:::.get_mean_vs_p0_from_matrix(expr.matrix)
    logistic_params = infercnv:::.get_logistic_params(mean_vs_p0_table)

    params[['dropout.logistic.midpt']] = logistic_params$midpt
    params[['dropout.logistic.slope']] = logistic_params$slope
        
    

    mean_vs_p0_table = cbind(mean_vs_p0_table, logm=log(mean_vs_p0_table$m + 1))
    smoothScatter(mean_vs_p0_table$logm, mean_vs_p0_table$p0, main=group_name)
    points(mean_vs_p0_table$logm,
           infercnv:::.logistic(mean_vs_p0_table$logm, logistic_params$midpt, logistic_params$slope), col='red')


    midpt_use = mean(mean_vs_p0_table$logm[mean_vs_p0_table$p0>0.48 & mean_vs_p0_table$p0<0.52])
    
    points(mean_vs_p0_table$logm,
           infercnv:::.logistic(mean_vs_p0_table$logm, midpt_use, logistic_params$slope), col='magenta')
    

    s = smooth.spline(mean_vs_p0_table$logm, mean_vs_p0_table$p0)
    r = range(mean_vs_p0_table$logm)
    x=seq(r[1], r[2], 0.1)
    points(x, predict(s, x)$y, col='orange')
    

    return(params)

}




# examine each group
all_groups = c(infercnv_obj@observation_grouped_cell_indices,  infercnv_obj@reference_grouped_cell_indices)
all_groups[['combined_normal']] <- unlist(infercnv_obj@reference_grouped_cell_indices)

for (group in names(all_groups)) {

    group_idxs = all_groups[[ group ]]
    expr.data = infercnv_obj@expr.data[,  group_idxs]

    params = get_parameters(group, expr.data)
    params = t(as.data.frame(params))
    
    print(params)
    
}

