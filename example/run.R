#!/usr/bin/env Rscript


options(error = function() traceback(2))

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(data_file="oligodendroglioma_expression_downsampled.txt",
                                    annotations_file="oligodendroglioma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="gencode_downsampled.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, 
                             out_path="output_dir", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             use_zscores=T,
                             )

# generate final plot
plot_cnv(infercnv_obj,
         out_dir="output_dir", 
         cluster_by_groups=T,
         color_safe_pal=FALSE,
         x.center=0,
         title="inferCNV",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         output_filename="infercnv")


