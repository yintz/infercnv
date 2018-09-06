#!/usr/bin/env Rscript


options(error = function() traceback(2))

library("infercnv")

infercnv_obj = CreateInfercnvObject(data_file="oligodendrocytoma_expression_downsampled.txt",
                                    annotations_file="oligodendrocytoma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="gencode_downsampled.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

ret = infercnv::run(infercnv_obj,
              cutoff=1, 
              out_path="output_dir", 
              cluster_by_groups=T, 
              plot_steps=F,
              use_zscores=T,
              )


