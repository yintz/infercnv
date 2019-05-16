#!/usr/bin/env Rscript

options(error = function() traceback(2))
options("warning.length" = 8000)

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 
                                    
out_dir="../example/output_dir_sampling_testscript"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             plot_steps=FALSE,
                             denoise=TRUE,
                             HMM=FALSE,
			                 no_prelim_plot=TRUE
                             )
t_out_dir = paste0(out_dir, "/subplots_1")
if(t_out_dir != "." & !file.exists(t_out_dir)){
    dir.create(t_out_dir)
}
infercnv:::plot_per_group(infercnv_obj, out_dir=t_out_dir, png_res=100, sample=TRUE, n_cells=100)


t_out_dir = paste0(out_dir, "/subsamples_1")
if(t_out_dir != "." & !file.exists(t_out_dir)){
    dir.create(t_out_dir)
}
sample_obj <- infercnv:::sample_object(infercnv_obj)

subsample_obj <- infercnv:::sample_object(infercnv_obj, n_cells=10)

upsubsample_obj <- infercnv:::sample_object(subsample_obj, n_cells=100)

every_2_object2 <- infercnv:::sample_object(infercnv_obj, every_n=2, above_m=2)

only_1_per_object <- infercnv:::sample_object(infercnv_obj, every_n=1000, above_m=2)

only_1_10times_per_object <- infercnv:::sample_object(only_1_per_object, n_cells=10)

infercnv_obj_filtered <- infercnv::apply_median_filtering(infercnv_obj, window_size=5)


infercnv::plot_cnv(sample_obj,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_sampled",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)

infercnv::plot_cnv(subsample_obj,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_subsampled",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)

infercnv::plot_cnv(upsubsample_obj,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_subsampled_then_upsampled",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)


infercnv::plot_cnv(every_2_object2,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_sample_every_2",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)



infercnv::plot_cnv(only_1_per_object,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_sample_only_1",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)


infercnv::plot_cnv(only_1_10times_per_object,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=t_out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_sample_only_1_10_times",
                   png_res=300,
                   output_format="png",
                   write_expr_matrix=TRUE)

infercnv::plot_cnv(infercnv_obj_filtered,
                   k_obs_groups=2,
                   cluster_by_groups=TRUE,
                   out_dir=out_dir,
                   x.center=1,
                   x.range="auto",
                   title="infercnv",
                   output_filename="infercnv_sampled_median_filtered",
                   png_res=300,
                   output_format = NA,
                   write_expr_matrix=TRUE)
