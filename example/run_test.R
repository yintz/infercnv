#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
    stop()
}

as.numeric(args[1])


options(error = function() traceback(2))

library("infercnv")

# create the infercnv object
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
#                                     annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
#                                     delim="\t",
#                                     gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
#                                     ref_group_names=NULL) 

infercnv_obj <- readRDS("default_input_infercnv_object.rds")

out_dir="output_dir_memory_test"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE,
                             plot_steps=FALSE,
                             no_plot=TRUE,
                             denoise=TRUE,
                             debug=TRUE,
                             HMM=TRUE,
                             up_to_step=as.numeric(args[1]))

