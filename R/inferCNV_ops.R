#!/usr/bin/env Rscript




#' Function doing the actual analysis before calling the plotting functions.
#'
#' @title Infer CNV changes given a matrix of RNASeq counts. Output a pdf and matrix of final values.
#'
#' @param data Expression matrix (genes X samples),
#'                 assumed to be log2(TPM+1) .
#' @param gene_order Ordering of the genes (data's rows)
#'                       according to their genomic location
#'                       To include all genes use 0.
#' @param cutoff Cut-off for the average expression of genes to be
#'                   used for CNV inference.
#' @param reference_obs Column names of the subset of samples (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all samples will
#'                          be the reference.
#' @param transform_data Indicator to log2 + 1 transform
#' @param window_length Length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#' @param max_centered_threshold The maximum value a a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value.
#' @param noise_filter The minimum difference a value can be from the
#'                            average reference in order for it not to be
#'                            cleared (zeroed out as noise).
#' @param noise_quantiles quantile range within the residual reference
#'                            distribution to be cleared (zeroed out as noise). Alternative
#'                            to param noise_filter.
#' @param ref_group_names Names of groups from the "annotations" table whose cells
#' are to be used as reference groups.
#' @param num_ref_groups The number of reference groups or a list of
#'                           indices for each group of reference indices in
#'                           relation to reference_obs.
#' @param out_path The path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
#' @param obs_annotations_groups Vector with group index of observations cells,
#' based on the annotation ot the cells.
#' @param k_obs_groups Number of groups in which to break the observations.
#' @param plot_steps If true turns on plotting intermediate steps.
#' @param contig_tail Length of the tail removed from the ends of contigs.
#' @param method_bound_vis Method to use for bounding values in the visualization.
#' @param lower_bound_vis Lower bound to normalize data to for visualization.
#' @param upper_bound_vis Upper bound to normalize data to for visualization.
#' @param ref_subtract_method Method used to subtract the reference values from the observations.
#' Valid choices are: "by_mean", "by_quantiles".
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' @param use_zscores If true, converts log(expression) data to zscores
#'
#' @return
#' Returns a list including:
#'     CNV matrix before visualization.
#'     CNV matrix after outlier removal for visualization.
#'     Contig order
#'     Column names of the subset of samples that should be used as references.
#'     Names of samples in reference groups.
#'
process_data <- function(infercnv_obj,
                      cutoff,
                      transform_data,
                      window_length,
                      max_centered_threshold,
                      noise_filter=NA,
                      noise_quantiles=c(0.025, 0.975),
                      num_ref_groups,
                      out_path,
                      grouping_key_coln,
                      cluster_by_groups,
                      k_obs_groups=1,
                      plot_steps=FALSE,
                      contig_tail= (window_length - 1) / 2,
                      method_bound_vis=NA,
                      lower_bound_vis=NA,
                      upper_bound_vis=NA,
                      ref_subtract_method="by_mean",
                      hclust_method='complete',
                      min_cells_per_gene=3,
                      use_zscores=FALSE,
                      make_zero_NA=FALSE) {

    logging::loginfo(paste("::process_data:Start", sep=""))

    # Split the reference data into groups if requested
    if (!is.null(num_ref_groups)) {
        ##TODO: update to use infercnv_obj
        groups_ref <- split_references(average_data=data, #data_smoothed,
                                       ref_obs=reference_obs,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)


        logging::loginfo(paste("::process_data:split_reference. ",
                               "found ",length(groups_ref)," reference groups.",
                               sep=""))
        
    }
    
            
    # Plot incremental steps.
    if (plot_steps){
        
        logging::loginfo(paste("\n\tSTEP 01: incoming data\n\n"))

        save(list=ls(), file=file.path(out_path, "01_incoming_data.Rdata"))

        plot_step(data=data,
                  plot_name=file.path(out_path,
                                      "01_incoming_data.pdf"))

        plot_cnv(plot_data=data,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="01_incoming_data",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.01_incoming_data")
    }


    # Make sure data is log transformed + 1
    if (transform_data){
        data <- log2(data + 1)


                                        # Plot incremental steps.
        if (plot_steps){
            logging::loginfo(paste("\n\tSTEP 02: log transformation of data\n\n"))

            save(list=ls(), file=file.path(out_path, "02_transformed.Rdata"))

            plot_step(data=data,
                      plot_name=file.path(out_path,
                                          "02_transformed.pdf"))

            plot_cnv(plot_data=data,
                     contigs=chr_order_for_plotting,
                     k_obs_groups=k_obs_groups,
                     obs_annotations_groups=obs_annotations_groups,
                     obs_annotations_names=obs_annotations_names,
                     grouping_key_coln=grouping_key_coln,
                     cluster_by_groups=cluster_by_groups,
                     reference_idx=ret_list[["REF_OBS_IDX"]],
                     ref_contig=NULL,
                     contig_cex=1,
                     ref_groups=ret_list[["REF_GROUPS"]],
                     ref_group_names=ref_group_names,
                     out_dir=out_path,
                     color_safe_pal=FALSE,
                     x.center=0,
                     title="02_log_transformed_data",
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename="infercnv.02_log_transformed"
                     )
        }
    }


    if (make_zero_NA) {
        data[data==0] = NA
    }

    ###################################################
    ## Step 03: removing insufficiently expressed genes
    
    # Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    # Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    #  with mean values below cutoff.

    genes_min_expr_cutoff <- above_min_mean_expr_cutoff(data, cutoff)

    ## require each gene to be present in a min number of cells for both obs and ref sets

    ## note, changed to just using the reference cells and not the observed here. See method for more details.
    genes_min_cells_obs_and_ref <- above_min_cells_obs_and_ref(data, min_cells_per_gene=min_cells_per_gene,
                                                               obs_idx=ret_list[["REF_OBS_IDX"]], ref_idx=unlist(ret_list[["REF_GROUPS"]]))

    # require both critiera:  min expression and min cell occupancy
    keep_gene_indices <- intersect(genes_min_expr_cutoff, genes_min_cells_obs_and_ref)



    if (!is.null(keep_gene_indices)){
        data <- data[keep_gene_indices, , drop=FALSE]
        gene_order <- gene_order[keep_gene_indices, , drop=FALSE]
        logging::loginfo(paste("::process_data:Reduce by removing genes below mean threshold, ",
                      "new dimensions (r,c) = ",
                      paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
        logging::logdebug(paste("::process_data:Keeping gene indices.", sep=""))


        # Plot incremental steps.
        chr_order_for_plotting <- paste(as.vector(as.matrix(gene_order[1])))
        if (plot_steps){

            logging::loginfo(paste("\n\tSTEP 03: Removing lowly expressed genes\n\n"))

            save(list=ls(), file=file.path(out_path, "03_reduced_by_cutoff.Rdata"))

            plot_step(data=data,
                      plot_name=file.path(out_path,
                                          "03_reduced_by_cutoff.pdf"))

            plot_cnv(plot_data=data,
                     contigs=chr_order_for_plotting,
                     k_obs_groups=k_obs_groups,
                     obs_annotations_groups=obs_annotations_groups,
                     obs_annotations_names=obs_annotations_names,
                     grouping_key_coln=grouping_key_coln,
                     cluster_by_groups=cluster_by_groups,
                     reference_idx=ret_list[["REF_OBS_IDX"]],
                     ref_contig=NULL,
                     contig_cex=1,
                     ref_groups=ret_list[["REF_GROUPS"]],
                     ref_group_names=ref_group_names,
                     out_dir=out_path,
                     color_safe_pal=FALSE,
                     x.center=0,
                     title="03_reduced_by_cutoff",
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename="infercnv.03_reduced_by_cutoff")

        }

    } else {
        logging::loginfo(paste("::process_data:Reduce by cutoff.", sep=""))
        logging::logwarn(paste("::No indicies left to keep.",
                               " Stoping."))
        stop(998)
    }


    # Reduce contig info
    chr_order <- gene_order[1]  # resetting
    chr_order_for_plotting = paste(as.vector(as.matrix(chr_order)))
    gene_order <- NULL



    ##########################################################################################
    ## Step 4: Centering data (w/ or w/o z-score transform) and max-centered threshold applied
    
    if (use_zscores) {

        # center and convert to z-scores
        logging::loginfo(paste("::center_and_Zscore_conversion", sep=""))

        # remember, genes are rows, cells are cols

        # centering and z-scores based on the reference (normal) cells:
        
        # ref data represent the null distribution
        ref_idx=unlist(ret_list[["REF_GROUPS"]])
        ref_data = data[,ref_idx]

        gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=T)})
        gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=T)})

        # center all genes at the ref (normal) center:
        data = sweep(data, 1, gene_ref_mean, FUN="-")

        # convert to z-scores based on the ref (normal) distribution
        data = sweep(data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.
        
    }
    else {
        # just center
        logging::loginfo(paste("::centering", sep=""))
        data <- sweep(data, 1, rowMeans(data, na.rm=T), FUN="-")
    }

    logging::loginfo(paste("::process_data:Outlier removal, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data, na.rm=TRUE),
                           " Min=", min(data, na.rm=TRUE),
                           " Max=", max(data, na.rm=TRUE),
                           ".", sep=""))

    #######################################################
    ## Apply maximum centered expression thresholds to data
    threshold = max_centered_threshold
    if (is.na(max_centered_threshold)) {

        threshold = mean(abs(get_average_bounds(data)))

        logging::loginfo(paste("::process_data:setting max centered expr threshold using quantiles, set to: +/-: ", threshold))
    }
    else {
        logging::loginfo(paste("::process_data:setting max centered expr threshold using specified settings of +/-: ", threshold))
    }

    # Cap values between threshold and -threshold, retaining earlier center
    data[data > threshold] <- threshold
    data[data < (-1 * threshold)] <- -1 * threshold

    # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 04: centering with max expression threshold\n\n"))

        save(list=ls(), file=file.path(out_path, "04_center_with_threshold.Rdata"))

        plot_step(data=data,
                  plot_name=file.path(out_path,
                                      "04_center_with_threshold.pdf"))

        plot_cnv(plot_data=data,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="04_center_with_threshold",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.04_center_with_threshold")

    }


    ###########################################################################
    # Step 5: For each cell, smooth the data along chromosome with gene windows

    logging::loginfo(paste("::process_data:Step 5: Smoothing per cell data by chromosome.", sep=""))
    data_smoothed <- smooth_window(data, window_length)
    #data_smoothed = data
    
    data <- NULL

    # Plot incremental steps.
    if (plot_steps){
        logging::loginfo(paste("\n\tSTEP 05: Smoothing data by chromosome\n\n"))
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "05_smoothed.pdf"))

        save(list=ls(), file=file.path(out_path, "05_smoothed.Rdata"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="05_smoothed",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.05_smoothed")
    }



    ## 
    # Step 6: 
    # Center cells/observations after smoothing. This helps reduce the
                                        # effect of complexity.
    logging::loginfo(paste("::Step 06: centering smoothed data by cell", sep=""))
    data_smoothed <- center_smoothed(data_smoothed)

    #                                    # recenter by gene now, so all genes centered at zero.
    #logging::loginfo(paste("::Step 06: recentering by gene", sep=""))
    #data_smoothed <- sweep(data_smoothed, 1, rowMeans(data_smoothed, na.rm=T), FUN="-")
    
    
    # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 06: re-centering data after smoothing\n\n"))

        save(list=ls(), file=file.path(out_path, "06_recentered.Rdata"))
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "06_recentered.pdf"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="06_centering_of_smoothed",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.06_centering_of_smoothed")

    }


    # Remove average reference
    i_ref_obs <- which(colnames(data_smoothed) %in% reference_obs)
    data_smoothed <- subtract_ref(average_data=data_smoothed,
                                  ref_groups=groups_ref,
                                  ref_subtract_method=ref_subtract_method)
    logging::loginfo(paste("::process_data:Remove average, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    # Plot incremental steps.
    if (plot_steps){


        logging::loginfo(paste("\n\tSTEP 07: removing average of reference data\n\n"))

        save(list=ls(), file=file.path(out_path, "07_remove_average.Rdata"))
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "07_remove_average.pdf"))
        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="07_remove_average",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.07_remove_average")

    }


    ## Step 08:
    # Remove Ends
    #logging::logdebug(c("chr_order: ", chr_order))
    #logging::logdebug(chr_order)
    remove_indices <- c()
    for (chr in unlist(unique(chr_order))){
        #logging::loginfo(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- remove_tails(data_smoothed,
                                             which(chr_order == chr),
                                             contig_tail)
        #logging::logdebug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))
        remove_indices <- c(remove_indices, remove_chr)

    }
    if (length(remove_indices) > 0){
        chr_order <- chr_order[remove_indices,]
        chr_order_for_plotting = paste(as.vector(as.matrix(chr_order)))
        data_smoothed <- data_smoothed[remove_indices, , drop=FALSE]
    }

                                        # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 08: removing genes at chr ends\n\n"))

        save(list=ls(), file=file.path(out_path, "08_remove_ends.Rdata"))

        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "08_remove_ends.pdf"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="08_remove_Ends",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.08_remove_ends")

    }
    logging::loginfo(paste("::process_data:Remove ends, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed, na.rm=TRUE),
                           " Min=", min(data_smoothed, na.rm=TRUE),
                           " Max=", max(data_smoothed, na.rm=TRUE),
                           ".", sep=""))


    ########################################
    # Step 08B: return to regular space from log space
    
    data_smoothed = 2^data_smoothed - 1

    # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 08B: inverting log transform\n\n"))

        save(list=ls(), file=file.path(out_path, "08B_inv_log_transform.Rdata"))

        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "08B_inv_log_transform.pdf"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="08B_invert_log_transform",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.08B_invert_log_transform")

    }
    
    
    ################################
    # Step 09: de-noising 

    if (! is.na(noise_filter)) {

        if (noise_filter > 0) {
            logging::loginfo(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
            data_smoothed <- clear_noise(smooth_matrix=data_smoothed,
                                         threshold=noise_filter,
                                         adjust_towards_zero=TRUE)
        }
        else {
                                        # noise == 0 or negative...
                                        # don't remove noise.
        }
        
    }
    else {
        # default, use quantiles, if NA 
        logging::loginfo(paste("::process_data:Remove noise, noise threshold defined via quantiles: ", noise_quantiles))
        data_smoothed <- clear_noise_via_ref_quantiles(smooth_matrix=data_smoothed,
                                                        ref_idx=unlist(ret_list[["REF_GROUPS"]]),
                                                       quantiles=noise_quantiles,
                                                       adjust_towards_zero=TRUE)
    }
    

    logging::loginfo(paste("::process_data:Remove noise, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed, na.rm=TRUE),
                           " Min=", min(data_smoothed, na.rm=TRUE),
                           " Max=", max(data_smoothed, na.rm=TRUE),
                           ".", sep=""))


                                        # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 09: Denoising\n\n"))

        save(list=ls(), file=file.path(out_path, "09_denoise.Rdata"))

        plot_step(data=data_smoothed,
                  plot_name=file.path(out_path,
                                      "09_denoise.pdf"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="09_denoised",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.09_denoised")

    }


    # Output before viz outlier
    ret_list[["PREVIZ"]] = data_smoothed

    # Remove outliers for viz
    remove_outlier_viz_pdf <- NA
    if (plot_steps){

        remove_outlier_viz_pdf <- file.path(out_path,
                                            "10A_remove_outlier.pdf")
    }

    data_smoothed = remove_outliers_norm(data=data_smoothed,
                                         out_method=method_bound_vis,
                                         lower_bound=lower_bound_vis,
                                         upper_bound=upper_bound_vis,
                                         plot_step=remove_outlier_viz_pdf)
    ret_list[["VIZ"]] <- data_smoothed

    
    # Plot incremental steps.
    if (plot_steps){

        logging::loginfo(paste("\n\tSTEP 10: Removing outliers\n\n"))

        save(list=ls(), file=file.path(out_path, "10B_remove_outlier.Rdata"))

        plot_step(data=ret_list[["VIZ"]],
                  plot_name=file.path(out_path,
                                      "10B_remove_outlier.pdf"))

        plot_cnv(plot_data=data_smoothed,
                 contigs=chr_order_for_plotting,
                 k_obs_groups=k_obs_groups,
                 obs_annotations_groups=obs_annotations_groups,
                 obs_annotations_names=obs_annotations_names,
                 grouping_key_coln=grouping_key_coln,
                 cluster_by_groups=cluster_by_groups,
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 ref_group_names=ref_group_names,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="10_removed_outliers",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.10_removed_outliers")
    }

    logging::loginfo(paste("::process_data:remove outliers, ",
                           "new dimensions (r,c) = ",
                           paste(dim(ret_list[["VIZ"]]), collapse=","),
                           " Total=", sum(ret_list[["VIZ"]]),
                           " Min=", min(ret_list[["VIZ"]]),
                           " Max=", max(ret_list[["VIZ"]]),
                           ".", sep=""))

    ret_list[["CONTIGS"]] = paste(as.vector(as.matrix(chr_order)))


    return(ret_list)
}





make_ngchm <- function(infercnv_obj, out_dir=".", title="NGCHM", gene_symbol=NULL, path_to_shaidyMapGen=NULL) {
    
    if (!is.null(path_to_shaidyMapGen)) {
        shaidy.path <- unlist(strsplit(path_to_shaidyMapGen, split = .Platform$file.sep))
        if (!file.exists(path_to_shaidyMapGen) || tail(shaidy.path, n = 1L) != "ShaidyMapGen.jar"){
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using the parameter \"path_to_shaidyMapGen\".",
                                   "Check that the correct pathway is being used.")
            logging::logerror(error_message)
            stop(error_message)
        }
    } else {
        path_to_shaidyMapGen <- Sys.getenv("SHAIDYMAPGEN")
        if (!file.exists(path_to_shaidyMapGen)){ ## check if envionrmental variable is passed
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using SHAIDYMAPGEN.",
                                   "Check that the correct pathway is being used.")
            logging::logerror(error_message)
            stop(error_message)
        }
    }
    
    if (!requireNamespace("NGCHM", quietly=TRUE)) {
        stop("The \"NGCHM\" library is required to use \"-ngchm=TRUE\" but it is not available.", .call=FALSE)
    }
    
    logging::loginfo("Creating NGCHM as infercnv.ngchm")
    Create_NGCHM(infercnv_obj,
                 path_to_shaidyMapGen = path_to_shaidyMapGen,
                 out_dir = output_dir,
                 title = fig_title,
                 gene_symbol = gene_symbol)
}



# USE_MEANS_FLAG = FALSE

# Remove the average of the genes of the reference observations from all
# observations' expression. Normalization by column.
#
# Args:
# average_data Matrix containing the data to remove average
#               from (this includes the reference observations).
#               Row = Genes, Col = Cells.
# ref_observations Indices of reference observations.
#                   Only these are used in the average.
# ref_groups A list of vectors of indices refering to the
#             different groups of the reference indices.
#
# ref_subtract_method method used to subtract reference data from obs.
#                      options are: "by_mean", "by_quantiles"  (default: "by_mean")
#
# quantiles reference quantiles to use if ref_subtract_method == 'by_quantiles'
#
# Returns:
# Expression with the average gene expression in the reference
#          observations removed.
subtract_ref <- function(average_data,
                         ref_groups,
                         ref_subtract_method="by_mean",
                         quantiles=c(0.25, 0.75)
                         ) {

                                        # r = genes, c = cells
    logging::loginfo(paste("::subtract_ref:Start", sep=""))
    # Max and min mean gene expression within reference groups.
    average_max <- NULL
    average_min <- NULL
    # average_reference_obs <- average_data[,ref_observations, drop=FALSE]
    # Reference gene within reference groups
    # now reference indexes of ref_groups are relative to the full average_data matrix and not the average_reference_obs references submatrix
    for (ref_group in ref_groups) {

        if (ref_subtract_method == "by_mean") {

            grp_average <- rowMeans(average_data[,ref_group, drop=FALSE], na.rm=TRUE)
            if(is.null(average_max)){
                average_max <- grp_average
            }
            if(is.null(average_min)){
                average_min <- grp_average
            }
            average_max <- pmax(average_max, grp_average)
            average_min <- pmin(average_min, grp_average)

        } else if (ref_subtract_method == "by_quantiles") {

            grp_expression_data = average_data[,ref_group, drop=FALSE, na.rm=TRUE]
            quants = x = apply(grp_expression_data, 1, function(x) { quantile(x, quantiles, na.rm=TRUE);})
            quants = t(x)
            q_low_bound = quants[,1]
            q_high_bound = quants[,2]
            if (is.null(average_min)) {
                average_min = q_low_bound
            } else {
                average_min = pmin(average_min, q_low_bound)
            }

            if (is.null(average_max)) {
                average_max = q_high_bound
            } else {
                average_max = pmax(average_max, q_high_bound)
            }
        } else {
            stop(paste("Error, unsupported ref_subtract_method specified: ", ref_subtract_method, sep=""))
        }
    }

    # Remove the Max and min averages (or quantiles) of the reference groups for each gene

    # debugging
    #ref_gene_group_means = data.frame(avg_min=average_min, avg_max=average_max);
    #write.table(ref_gene_group_means, file.path(out_path, "ref_gene_group_means.dat"), quote=F, sep="\t")


    for(gene_i in 1:nrow(average_data)){
        current_col <- average_data[gene_i, ]

        # original code
        i_max <- which(current_col > average_max[gene_i])
        i_min <- which(current_col < average_min[gene_i])

        row_init <- rep(0, length(current_col))
        if(length(i_max) > 0){
            row_init[i_max] <- current_col[i_max] - average_max[gene_i]
        }
        if(length(i_min) > 0){
            row_init[i_min] <- current_col[i_min] - average_min[gene_i]
        }
        average_data[gene_i, ] <- row_init
    }

    return(average_data)
}


# Not testing, parameters ok.
# Helper function allowing greater control over the steps in a color palette.
# Source:http://menugget.blogspot.com/2011/11/define-color-steps-for-
#               colorramppalette.html#more

# Args:
# steps Vector of colors to change use in the palette
# between: Steps where gradients change
#
# Returns:
# Color palette
color.palette <- function(steps,
                          between=NULL, ...){

    if (is.null(between)){
        between <- rep(0, (length(steps) - 1))
    }
    if (length(between) != length(steps) - 1){
        stop("Must have one less \"between\" value than steps")
    }

    fill.steps <- cumsum(rep(1, length(steps)) + c(0, between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[, fill.steps] <- col2rgb(steps)
    for(i in which(between > 0)){
        col.start <- RGB[, fill.steps[i]]
        col.end <- RGB[, fill.steps[i + 1]]
        for(j in seq(3)){
            vals <- seq(col.start[j],
                        col.end[j],
                        length.out=between[i] + 2)
            vals <- vals[2:(2 + between[i] - 1)]
            RGB[j, (fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue=255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Create a sepList forthe heatmap.3 plotting function given integer vectors
# of rows and columns where speration should take place.
# The expected input to the heatmap function is a list of 2 lists.
# The first list are column based rectangles, and the second row.
# To define a rectagle the index of the row or column where the line of the rectagle
# should be placed is done with a vector of integers, left, bottom, right and top line.
# Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
#
# Args:
# row_count Total number of rows
# col_count Total number of columns
# row_seps Vector of integers indices for row breaks
# col_seps Vector of integer indices for column breaks
#
# Returns
# List of lists of vectors
create_sep_list <- function(row_count,
                            col_count,
                            row_seps=NULL,
                            col_seps=NULL){
    sepList <- list()
    # Heatmap.3 wants a list of boxes for seperating columns
    # Column data
    if(!is.null(col_seps) &&
       !is.na(col_seps) &&
       (length(col_seps)>0) &&
       col_count > 0){
        colList <- list()
        for(sep in 1:length(col_seps)){
            colList[[sep]] <- c(col_seps[sep],0,col_seps[sep],row_count)
        }
        sepList[[1]] <- colList
    } else {
        sepList[[1]] <- list()
        sepList[[1]][[1]] <- c(0,0,0,0)
    }

    # Row data
    # This is measured from bottom to top
    # So you have to adjust the values of the data
    row_seps <- row_count-row_seps
    if(!is.null(row_seps) &&
       !is.na(row_seps) &&
       (length(row_seps)>0) &&
       row_count > 0){
        rowList <- list()
        for(sep in 1:length(row_seps)){
            rowList[[sep]] <- c(0,row_seps[sep],col_count,row_seps[sep])
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}

# Split up reference observations in to k groups and return indices
# for the different groups.
#
# Args:
# average_data Matrix containing data. Row = Genes, Col = Cells.
# ref_obs Indices of reference obervations.
# num_groups The number of groups to partition nodes in or a list
#                       of already partitioned indices.
#
# Returns:
# Returns a list of grouped reference observations given as
#             vectors of groups. These are indices relative to the reference
#             observations only, so a return 1 indicates the first reference
#             row, not the first row.
split_references <- function(average_data,
                             ref_obs,
                             num_groups,
                             hclust_method='complete') {
    logging::loginfo(paste("::split_references:Start", sep=""))
    ret_groups <- list()
    split_groups <- NULL
    if(is.null(average_data)){
        return(ret_groups)
    }
    if(is.null(num_groups)){
        num_groups <- 1
    }
    if(is.null(ref_obs)){
        ref_obs <- 1:ncol(average_data)
    }

    # If a supervised grouping is given as a commandline argument
    # num_groups will be a list of length of more than 1.
    # Otherwise do this.
    if (!is.list(num_groups)){
        # num_groups <- unlist(num_groups)
        if (num_groups > ncol(average_data)){
            num_groups <- ncol(average_data)
        }

        # If only one group is needed, short-circuit and
        # return all indices immediately
        # If only one group is asked for or only one reference is
        # available to give.
        if ((num_groups < 2) ||
           (length(ref_obs) < 2)){
            ret_groups[[1]] <- ref_obs
            return(ret_groups)
        }

        # Get HCLUST
        # Get reference observations only.
        average_reference_obs <- t(average_data[, ref_obs, drop=FALSE])

        hc <- hclust(dist(average_reference_obs), method=hclust_method)

        split_groups <- cutree(hc, k=num_groups)
        split_groups <- split_groups[hc$order]
        # Keep the sort of the hclust
        for(cut_group in unique(split_groups)){
            group_idx <- which(split_groups == cut_group)
            ret_groups[[cut_group]] <- ref_obs[hc$order[group_idx]]
            # ret_groups[[cut_group]] <- which(colnames(average_data) %in% names(group_idx))
        }
    } else {
        ret_groups <- num_groups
    }
    return(ret_groups)
}

# Set outliers to some upper or lower bound. Then normalize values to
# approximately [-1, 1]. This is to prep the data for visualization.
#
# Args:
# data: data to remove outliers. Outliers removed within columns.
# out_method Method to remove outliers [(average_bound, NA (hard threshold))]
# lower_bound Lower bound which identifies a measurement
#                        as an outlier.
# upper_bound Upper bound which identifies a measurement
#                        as an outlier.
# plot_step: True will plot this analysis step.
#
# Returns:
# Return data matrix with outliers removed
remove_outliers_norm <- function(data,
                                 out_method=NA,
                                 lower_bound=NA,
                                 upper_bound=NA,
                                 plot_step=NA) {
    logging::loginfo(paste("::remove_outlier_norm:Start",
                           "out_method:", out_method,
                           "lower_bound:" , lower_bound,
                           "upper_bound:", upper_bound,
                           "plot_step:" , plot_step))
    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        logging::logerror("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }
    if (is.na(lower_bound) || is.na(upper_bound)){
                                        # using out_method instead of specified bounds.
        logging::loginfo(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))
        if(is.na(out_method)){
            logging::loginfo("::remove_outlier_norm:WARNING outlier removal was not performed.")
            return(data)
        }
        if (out_method == "average_bound"){

            bounds = get_average_bounds(data)
            lower_bound = bounds[1]
            upper_bound = bounds[2]

            # Plot bounds on data
            if(!is.na(plot_step)){
                pdf(plot_step, useDingbats=FALSE)
                boxplot(data)
                points(1:ncol(data), rep(lower_bound, ncol(data)),
                       pch=19, col="orange")
                points(1:ncol(data), rep(upper_bound, ncol(data)),
                       pch=19, col="orange")
                dev.off()
            }

        } else {
            logging::logerror(paste("::remove_outlier_norm:Error, please",
                                    "provide an approved method for outlier",
                                    "removal for visualization."))
            stop(991)
        }
    }
    else {
        # Hard threshold given bounds
        logging::loginfo(paste("::remove_outlier_norm:",
                               "lower_bound:" , lower_bound,
                               "upper_bound:", upper_bound) )

        if(!is.na(plot_step)){
            pdf(plot_step, useDingbats=FALSE)
            boxplot(data)
            points(1:ncol(data), rep(lower_bound, ncol(data)),
                   pch=19, col="orange")
            points(1:ncol(data), rep(upper_bound, ncol(data)),
                   pch=19, col="orange")
            dev.off()
        }

    }

    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound

    return(data)

}

# Center data after smoothing. Center with in cells using median.
#
# Args:
# data_smoothed Matrix to center.
#                          Row = Genes, Col = cells.
#
# Returns:
# Matrix that is median centered.
#             Row = Genes, Col = cells.
center_smoothed <- function(data_smoothed){

    logging::loginfo(paste("::center_smoothed:Start"))

    # Center within columns (cells)
    row_median <- apply(data_smoothed, 2, function(x) { median(x, na.rm=T) } )

    return(t(apply(data_smoothed, 1, "-", row_median)))
}



# Return the indices of the rows that average above the cut off
#
# Args:
# data Data to measure the average row and evaluate
#                 against the cutoff. Row = Genes, Col = Cells.
# cutoff Threshold to be above to be kept.
#
# Returns:
# Returns a vector of row indicies to keep (are above the cutoff).
above_min_mean_expr_cutoff <- function(data, cutoff){

    logging::loginfo(paste("::above_min_mean_expr_cutoff:Start", sep=""))
    average_gene <- log2(rowMeans( ( (2 ^ data) - 1), na.rm=TRUE) + 1 )
    logging::loginfo(paste("::process_data:Averages (counts).", sep=""))
    # Find averages above a certain threshold
    indicies <- which(average_gene > cutoff)
    if (length(indicies) > 0){
        return(indicies)
    } else {
        return(NULL)
    }
}


#' indicate which genes (rows) have at least specified min_cells_per_gene
#'
#' Args
#' @param data Data (expression) matrix
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#' @param obs_idx vector containing the column indices for the observed (tumor) cells
#' @param ref_idx vector containing the column indices for teh reference (normal) cells

above_min_cells_obs_and_ref = function(data, min_cells_per_gene, obs_idx, ref_idx) {

    ref_data = data[,ref_idx]

    ref_genes_passed = which(apply(ref_data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))


    #### chromosomes lost in the observed (tumor) may have no expression and don't want to lose those cells via filtering!

    ## require expression in reference, needed for determining gain / loss.

    #obs_data = data[,obs_idx]

    #obs_genes_passed = which(apply(obs_data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))

    #both_passed = intersect(ref_genes_passed, obs_genes_passed)

    #return(both_passed)

    return(ref_genes_passed)

}



#' Order the data and subset the data to data in the genomic position file.
#'
#' Args:
#' @param data Data (expression) matrix where the row names should be in
#'                 the row names of the genomic_position file.
#' @param genomic_position Data frame read in from the genomic position file
#'
#' @return Returns a matrix of expression in the order of the
#'            genomic_position file. NULL is returned if the genes in both
#'            data parameters do not match.
#'
order_reduce <- function(data, genomic_position){
    logging::loginfo(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL, chr_order=NULL)
    if (is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }

    # Drop pos_gen entries that are position 0
    remove_by_position <- -1 * which(genomic_position[2] + genomic_position[3] == 0)
    if (length(remove_by_position)){
        logging::logdebug(paste("::process_data:order_reduce: removing genes specified by pos == 0, count: ",
                                length(remove_by_position), sep=""))

        genomic_position <- genomic_position[remove_by_position, , drop=FALSE]
    }

    # Reduce to genes in pos file

    logging::logdebug(paste("::process_data:order_reduce: gene identifers in expression matrix: ",
                            row.names(data), collapse="\n", sep=""))
    logging::logdebug(paste("::process_data:order_reduce: gene identifers in genomic position table: ",
                            row.names(data), collapse="\n", sep=""))



    keep_genes <- row.names(data)[which(row.names(data)
                                  %in% row.names(genomic_position))]
    logging::logdebug(paste("::process_data:order_reduce: keep_genes size: ", length(keep_genes),
                            sep=""))

    # Keep genes found in position file
    if(length(keep_genes)){
        ret_results$expr <- data[keep_genes, , drop=FALSE]
        ret_results$order <- genomic_position[keep_genes, , drop=FALSE]
    } else {
        logging::loginfo(paste("::process_data:order_reduce:The position file ",
                               "and the expression file row (gene) names do not match."))
        return(list(expr=NULL, order=NULL, chr_order=NULL))
    }

    # Set the chr to factor so the order can be arbitrarily set and sorted.
    chr_levels <- unique(genomic_position[[CHR]])
    ret_results$order[[CHR]] <- factor(ret_results$order[[CHR]],
                                   levels=chr_levels)

    # Sort genomic position file and expression file to genomic position file
    # Order genes by genomic region
    order_names <- row.names(ret_results$order)[with(ret_results$order, order(chr,start,stop))]
    ret_results$expr <- ret_results$expr[order_names, , drop=FALSE]

    # This is the contig order, will be used in visualization.
    # Get the contig order in the same order as the genes.
    ret_results$order <- ret_results$order[order_names, , drop=FALSE]
    ret_results$chr_order <- ret_results$order[1]

    # Remove any gene without position information
    # Genes may be sorted correctly by not have position information
    # Here they are removed.
    logging::loginfo(paste("::process_data:order_reduce:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    logging::logdebug(paste("::process_data:order_reduce end."))
    return(ret_results)
}

# Remove values that are too close to the average and are considered noise.
#
# Args:
# smooth_matrix A matrix of values, smoothed, and with average
#                          reference removed. Row = Genes, Col = Cells.
# threshold The amount of difference a value must be from the
#                      reference before the value can be kept and not
#                      removed as noise.
# Returns:
# Denoised matrix
clear_noise <- function(smooth_matrix, threshold, adjust_towards_zero){
        
    logging::loginfo(paste("********* ::clear_noise:Start. threshold: ", threshold,  " adj_towards_zero: ", adjust_towards_zero, sep=""))
        
    if (threshold > 0){

        if (adjust_towards_zero) {
            
            upper_noise_flags = (smooth_matrix > 0 & smooth_matrix <= threshold)
            smooth_matrix[upper_noise_flags] = smooth_matrix[upper_noise_flags] - threshold
            smooth_matrix[ smooth_matrix[upper_noise_flags] < 0 ] = 0

            lower_noise_flags = (smooth_matrix < 0 & smooth_matrix >= -1*threshold)
            smooth_matrix[lower_noise_flags] = smooth_matrix[lower_noise_flags] + threshold
            smooth_matrix[ smooth_matrix[lower_noise_flags] > 0 ] = 0
                        
        }
        else {
            smooth_matrix[abs(smooth_matrix) < threshold] <- 0
        }
    }
    return(smooth_matrix)
}

# clear_noise_via_ref_quantiles: define noise levels based on quantiles within the ref (normal cell) distribution.
# Any data points within this defined quantile are set to zero.

clear_noise_via_ref_quantiles <- function(smooth_matrix, ref_idx, quantiles=c(0.025, 0.975), adjust_towards_zero=TRUE) {

    #save(list=c('smooth_matrix', 'ref_idx'), file='ladeda.rdata')
    
    vals = smooth_matrix[,ref_idx]

    vals[vals==0] = NA  # use remaining ref vals that weren't already turned to zeros

    lower_quantile = quantiles[1]
    upper_quantile = quantiles[2]

    #logging::loginfo(paste("::clear_noise_via_ref_quantiles: using noise quantiles set at: ",
    #                       lower_quantile, "-", upper_quantile, sep=""))

    #lower_bound <- mean(apply(vals, 2,
    #                          function(x) quantile(x, probs=lower_quantile, na.rm=TRUE)))

    #upper_bound <- mean(apply(vals, 2,
    #                          function(x) quantile(x, probs=upper_quantile, na.rm=TRUE)))

    upper_bound <- mean(apply(vals, 2, function(x) sd(x, na.rm=T))) * 1.5
    
    lower_bound <- -1 * upper_bound
    
    
    logging::loginfo(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                           lower_bound, "-", upper_bound, sep=" "))


    if (adjust_towards_zero) {
        if (upper_bound > 0) {
            upper_noise_flags = (smooth_matrix > 0 & smooth_matrix <= upper_bound)
            smooth_matrix[upper_noise_flags] = smooth_matrix[upper_noise_flags] - upper_bound
            smooth_matrix[ upper_noise_flags & smooth_matrix < 0 ] = 0 # dealing w/ over-correction
        }
        
        if (lower_bound < 0) {
            lower_noise_flags = (smooth_matrix < 0 & smooth_matrix > lower_bound)
            smooth_matrix[lower_noise_flags] = smooth_matrix[lower_noise_flags] - lower_bound # subtracting a negative val, so making more pos
            smooth_matrix[ lower_noise_flags & smooth_matrix > 0 ] = 0 # dealing w/ over-correction
        }
        
    }
    else {
        smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = 0
    }
    return(smooth_matrix)
}


# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr Indices of the chr in which the tails are to be removed.
# tail_length Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.
remove_tails <- function(smooth_matrix, chr, tail_length){

    #logging::loginfo(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
         tail_length <- floor(chr_length / 3)
    }
    remove_indices <- -1 * chr[1:tail_length]
    remove_indices <- c(remove_indices,
                        -1 * (chr[ ( (chr_length + 1) - tail_length):
                                   chr_length]))
    return(remove_indices)
}

# Smooth a matrix by column using a simple moving average.
# Tails of the averages use a window length that is truncated to
# available data.
#
# Args:
# data Data matrix to smooth. Row = Genes, Col = Cells.
# window_length Length of window to use for the moving average.
#        Should be a positive, odd integer.
#
# Returns:
# Matrix with columns smoothed with a simple moving average.
smooth_window <- function(data, window_length, smooth_ends=TRUE){

    logging::loginfo(paste("::smooth_window:Start.", sep=""))
    if (window_length < 2){
        return(data)
    }
    if (window_length > nrow(data)){
        return(data)
    }

    tail_length <- (window_length - 1) / 2
    num_genes <- nrow(data)
    data_sm <- apply(data,
                     2,
                     smooth_window_helper,
                     window_length=window_length)
    logging::logdebug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))

    if (smooth_ends) {
        # Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
        data_sm <- apply(data_sm,
                         2,
                         smooth_ends_helper,
                         tail_length=tail_length)

    }
    
    # Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    return(data_sm)
}

# Helper function for smoothing the ends of a moving average.
#
# Args:
# obs_data: Data to smooth
# tail_length:  Length of the tail to smooth.
#
# Returns:
# Data smoothed.
smooth_ends_helper <- function(obs_data, tail_length) {

    # strip NAs out and replace after smoothing
    orig_obs_data = obs_data

    nas = is.na(obs_data)

    obs_data = obs_data[!nas]

    obs_length <- length(obs_data)
    end_data <- obs_data

    # end_data will have the end positions replaced with mean values, smoothing just at the ends.

    obs_count <- length(obs_data)

    for (tail_end in 1:tail_length) {

        # algorithm generates smoothing windows from the end like so:
        # <|>
        # < | >
        # <  |  >
        # <   |   >
        # <    |    >
        # where | is the central position assigned the mean of observations in the window.

        bounds <- tail_end - 1
        end_tail <- obs_count - bounds

        show_debug_logging_here = FALSE

        if (show_debug_logging_here) {
            logging::logdebug(paste("::smooth_ends_helper: tail range <",
                                    tail_end - bounds,
                                    "|", tail_end, "|",
                                    tail_end + bounds,">", sep=" "))
        }

        end_data[tail_end] <- mean(obs_data[(tail_end - bounds):
                                            (tail_end + bounds)],
                                   na.rm=TRUE)


        if (show_debug_logging_here) {
            logging::logdebug(paste("::smooth_ends_helper: tail range <",
                                    end_tail - bounds,
                                    "|", end_tail, "|",
                                    end_tail + bounds, ">", sep=" "))
        }

        end_data[end_tail] <- mean(obs_data[(end_tail - bounds):
                                            (end_tail + bounds)],
                                   na.rm=TRUE)
    }

    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data

    return(orig_obs_data)
}

# Smooth vector of values over the given window length.
#
# Args:
# obs_data Vector of data to smooth with a moving average.
# window_length Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.
smooth_window_helper <- function(obs_data, window_length){

    nas = is.na(obs_data)
    vals = obs_data[! nas]

    smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)

    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]

    obs_data[! nas] = vals

    return(obs_data)
}




get_average_bounds = function (data) {

    lower_bound <- mean(apply(data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))

    return(c(lower_bound, upper_bound))

}
