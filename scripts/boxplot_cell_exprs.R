#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--log", help="log transform expr", action='store_true', default=FALSE)

args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(tidyverse)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

expr.data = infercnv_obj@expr.data

## build df of expr values.
cell_groups = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

cell_group_names = names(cell_groups)


pngname = sprintf("%s-boxplot.png", infercnv_obj_file)
png(pngname)

expr_tibble = do.call(rbind, lapply(cell_group_names, function(cell_group_name) {
    cell_group_expr = expr.data[, cell_groups[[ cell_group_name ]] ]
    
    cell_group_expr = as.tibble(cell_group_expr)
    
    cell_group_expr = cell_group_expr %>% gather(key='cellname', value='expr')
    
    cell_group_expr = cell_group_expr %>% mutate(group_name=cell_group_name)
}))



p = expr_tibble %>% ggplot(aes(y=expr, x=cellname, color=group_name)) + geom_boxplot(outlier.shape=NA) + facet_wrap(~group_name, scales='free_x')

plot(p)

saveRDS(expr_tibble, 'my.tibble.obj')

