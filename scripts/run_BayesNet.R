#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--prelim_infercnv_obj", help="preliminary infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--i6HMM_infercnv_obj", help="i6HMM infercnv_obj file", required=TRUE, nargs=1)

parser$add_argument("--BayesMaxPNormal", help="BayesMaxPNormal", required=TRUE, nargs=1, type='double')
parser$add_argument("--out_dir", help="output directory", required=TRUE, nargs=1)

args = parser$parse_args()

library(infercnv)
library(futile.logger)

infercnv_obj_prelim = readRDS(args$prelim_infercnv_obj)

hmm.infercnv_obj = readRDS(args$i6HMM_infercnv_obj)


flog.info("Running Bayesian Network Model on HMM predicted CNV's\n")

hmm.infercnv_obj <- infercnv::inferCNVBayesNet(infercnv_obj    = infercnv_obj_prelim,
                                               HMM_obj         = hmm.infercnv_obj,
                                               BayesMaxPNormal = args$BayesMaxPNormal,
                                               file_dir        = args$out_dir,
                                               postMcmcMethod  = "removeCNV",
                                               out_dir         = file.path(args$out_dir, "BayesNetOutput"),
                                               quietly = TRUE)



