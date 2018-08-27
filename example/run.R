#!/usr/bin/env Rscript

library("infercnv")

infercnv::infercnv(x="oligodendrocytoma_expression_downsampled.txt", 
                   gene_order="gencode_downsampled.txt", 
                   annotations="oligodendrocytoma_annotations_downsampled.txt",
                   cutoff=1, 
                   #noise_filter=0.2,
                   #max_centered_expression=3,
                   output_dir="outdir", 
                   hclust_method="ward.D",
                   name_ref_groups="Microglia/Macrophage,Oligodendrocytes (non-malignant)", 
                   cluster_by_groups=T, 
                   plot_steps=F,
                   use_zscores=T,
                   log_level="DEBUG"
                   )

