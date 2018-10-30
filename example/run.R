#!/usr/bin/env Rscript

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="oligodendroglioma_expression_downsampled.counts.matrix",
                                    annotations_file="oligodendroglioma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="gencode_downsampled.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

out_dir="output_dir"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             include.spike=T
                             )


# save the final object in case we want to experiment with it later, replotting w/ diff thresholds, etc.
save('infercnv_obj', file=file.path(out_dir, 'infercnv.final.obj'))

# generate final plot
plot_cnv(infercnv_obj,
         out_dir=out_dir, 
         cluster_by_groups=T,
         color_safe_pal=FALSE,
         x.center=1,
         x.range=c(0,2),
         title="inferCNV",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         output_filename="infercnv")


