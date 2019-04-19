#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args) == 0) {
    stop("Error, require params: infercnv.obj");
}

infercnv_obj_file = args[1]

infercnv_obj = readRDS(infercnv_obj_file)


library(edgeR)
library(fitdistrplus)
library(infercnv)
library(Matrix)

# borrowing some code from splatter

get_parameters <- function(group_name, expr.matrix) {

    message(sprintf("getting params for: %s", group_name))
    params = list()
    params[['group_name']] = group_name
    
    
    # estimate common dispersion
    design <- matrix(1, ncol(expr.matrix), 1)
    disps <- edgeR::estimateDisp(expr.matrix, design = design)

    params[[ 'common.dispersion' ]] = disps$common.dispersion
    

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

