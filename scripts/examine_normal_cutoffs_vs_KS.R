#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--scale", help="scale", action='store_true', default=FALSE)
parser$add_argument("--subtract", help="subtract", action='store_true', default=FALSE)
parser$add_argument("--smooth", help="smooth", action='store_true', default=TRUE)
parser$add_argument("--show_tumor", help="show tumor instead of normal", action='store_true', default=FALSE)
parser$add_argument("--output", help="name of output png file", required=TRUE)

args = parser$parse_args()

library(infercnv)
library(tidyverse)
library(futile.logger)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

if (! infercnv:::has_reference_cells(infercnv_obj)) {
    stop("Error, cannot tune parameters without reference 'normal' cells defined")
}

if (args$scale) {
    infercnv_obj <- infercnv:::scale_infercnv_expr(infercnv_obj)
}

if (args$subtract) {
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE)
}


if (args$smooth) {
    infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=101, smooth_ends=TRUE)
}

if (args$show_tumor) {
    expr_vals <- infercnv_obj@expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)]
} else {
    expr_vals <- infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)]
}


mu = mean(expr_vals)
sigma = sd(expr_vals)

data.want = data.frame(vals=as.numeric(expr_vals))

mean_delta = infercnv:::determine_mean_delta_via_Z(sigma, p=0.05)
KS_delta = infercnv:::get_HoneyBADGER_setGexpDev(gexp.sd=sigma, alpha=0.05)


png(args$output)

message("plotting ncells distribution")

message("mean delta: ", mean_delta)
message("KS_delta: ", KS_delta)

p = data.want %>% ggplot(aes(vals)) +
    geom_density(alpha=0.3)

p = p +
    stat_function(fun=dnorm, color='black', args=list('mean'=mu,'sd'=sigma))


## add Z-based

p = p +
    stat_function(fun=dnorm, color='blue', args=list('mean'=mu-mean_delta,'sd'=sigma)) +
    stat_function(fun=dnorm, color='blue', args=list('mean'=mu+mean_delta,'sd'=sigma))

## add KS-based

p = p +
    stat_function(fun=dnorm, color='magenta', args=list('mean'=mu-KS_delta,'sd'=sigma)) +
    stat_function(fun=dnorm, color='magenta', args=list('mean'=mu+KS_delta,'sd'=sigma))

plot(p)
