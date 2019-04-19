#!/usr/bin/env Rscript

options(error = function() traceback(2))

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../inst/extdata/oligodendroglioma_expression_downsampled.counts.matrix.gz",
                                    annotations_file="../inst/extdata/oligodendroglioma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="../inst/extdata/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

out_dir="output_dir"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=TRUE
                             )

