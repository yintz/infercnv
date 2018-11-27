#!/usr/bin/env Rscript

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="merged.matrix",
                                    annotations_file="my.cell.annots",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("normal"))

out_dir="output_dir"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             plot_steps=T,
                             include.spike=T,  # used for final scaling to fit range (0,2) centered at 1.
                             HMM=T
                             )

