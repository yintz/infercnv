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
run <- function(infercnv_obj,
                cutoff=1,
                out_path=".",
                transform_data=FALSE,
                window_length=101,
                num_ref_groups=NULL,
                max_centered_threshold=NA,
                noise_filter=NA,
                cluster_by_groups=FALSE,
                k_obs_groups=1,
                plot_steps=FALSE,
                method_bound_vis="average_bound",
                lower_bound_vis=NA,
                upper_bound_vis=NA,
                ref_subtract_method="by_mean",
                hclust_method='complete',
                min_cells_per_gene=3,
                sd_amplifier = 1.5,
                use_zscores=FALSE,
                make_zero_NA=FALSE) {
    
    flog.info(paste("::process_data:Start", sep=""))

    if(out_path != "." & !file.exists(out_path)){
        dir.create(out_path)
    }



    flog.info(paste("\n\n\tSTEP 01: incoming data\n"))

    # Split the reference data into groups if requested
    if (!is.null(num_ref_groups)) {
        ##TODO: update to use infercnv_obj
        groups_ref <- split_references(average_data=data, #data_smoothed,
                                       ref_obs=reference_obs,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)


        flog.info(paste("::process_data:split_reference. ",
                               "found ",length(groups_ref)," reference groups.",
                               sep=""))
        
    }
    
    
    # Plot incremental steps.
    if (plot_steps) {        

        infercnv_obj_01 <- infercnv_obj
        save('infercnv_obj_01', file=file.path(out_path, "01_incoming_data.infercnv_obj"))
        
        plot_cnv(infercnv_obj=infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="01_incoming_data",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.01_incoming_data",
                 write_expr_matrix=TRUE)
    }
    
    
    # Make sure data is log transformed + 1
    if (transform_data){

        flog.info(paste("\n\n\tSTEP 02: log transformation of data\n"))
        
        infercnv_obj <- log2xplus1(infercnv_obj)
                
        # Plot incremental steps.
        if (plot_steps){

            infercnv_obj_02 <- infercnv_obj
            save('infercnv_obj_02', file=file.path(out_path, "02_logtransformed.infercnv_obj"))

            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_path,
                     color_safe_pal=FALSE,
                     x.center=0,
                     title="02_log_transformed_data",
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename="infercnv.02_log_transformed",
                     write_expr_matrix=TRUE
                     )
        }
    }
    

    if (make_zero_NA) {
        infercnv_obj <- make_zero_NA(infercnv_obj)
    }
    
    ###################################################
    ## Step 03: removing insufficiently expressed genes

    flog.info(paste("\n\n\tSTEP 03: Removing lowly expressed genes\n"))
        
    # Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    # Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    #  with mean values below cutoff.

    infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, cutoff)
    
    ## require each gene to be present in a min number of cells for ref sets

    infercnv_obj <- require_above_min_cells_ref(infercnv_obj, min_cells_per_gene=min_cells_per_gene)


    if (plot_steps){
        
        infercnv_obj_03 <- infercnv_obj
        
        save('infercnv_obj_03', file=file.path(out_path, "03_reduced_by_cutoff.infercnv_obj"))
        
        plot_cnv(infercnv_obj=infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="03_reduced_by_cutoff",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.03_reduced_by_cutoff",
                 write_expr_matrix=TRUE)
        
    }
    

    ##########################################################################################
    ## Step 4: Centering data (w/ or w/o z-score transform) and max-centered threshold applied

    flog.info(paste("\n\n\tSTEP 04a: centering gene expression\n"))
        
    if (use_zscores) {

        infercnv_obj <- transform_to_reference_based_Zscores(infercnv_obj)
        
    }
    else {
        # just center
        infercnv_obj <- mean_center_gene_expr(infercnv_obj)
    }

    if (plot_steps){
        
        infercnv_obj_04a <- infercnv_obj
        
        save('infercnv_obj_04a', file=file.path(out_path, "04a_centering_gene_expr.infercnv_obj"))
        
        plot_cnv(infercnv_obj=infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="04a_centering_gene_expr",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.04a_centering_gene_expr",
                 write_expr_matrix=TRUE)
        
    }


    
        
    #######################################################
    ## Apply maximum centered expression thresholds to data
    # Cap values between threshold and -threshold, retaining earlier center

    flog.info(paste("\n\n\tSTEP 04b: apply max centered expression threshold\n"))
    
    threshold = max_centered_threshold
    if (is.na(max_centered_threshold)) {
        threshold = mean(abs(get_average_bounds(infercnv_obj)))
    }
    
    infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, threshold=threshold)
    
    # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_04b <- infercnv_obj
        
        save('infercnv_obj_04b', file=file.path(out_path, "04b_apply_max_centered_expr_threshold.infercnv_obj"))

        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="04b_apply_max_centered_expr_threshold",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.04b_apply_max_centred_expr_threshold",
                 write_expr_matrix=TRUE)
        
    }
    

    ###########################################################################
    # Step 5: For each cell, smooth the data along chromosome with gene windows
    
    flog.info(paste("\n\n\tSTEP 05: Smoothing data per cell by chromosome\n"))

    infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=TRUE)
    
    
    # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_05 <- infercnv_obj
        save('infercnv_obj_05', file=file.path(out_path, "05_smoothed_by_chr.infercnv_obj"))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="05_smoothed_by_chr",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.05_smoothed_by_chr")
    }
    


    ## 
    # Step 6: 
    # Center cells/observations after smoothing. This helps reduce the
                                        # effect of complexity.

    
    flog.info("\n\n\tSTEP 06: re-centering data across chromosome after smoothing\n")
    
    infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj)
    
    
    # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_06 <- infercnv_obj
        
        save('infercnv_obj_06', file=file.path(out_path, "06_recentered_cells_by_chr.infercnv_obj"))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="06_centering_of_smoothed",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.06_centering_of_smoothed")
        
    }


    ####################################
    ## Step 07: Remove average reference

    flog.info("\n\n\tSTEP 07: removing average of reference data\n")
        
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, method=ref_subtract_method)
    
    
    # Plot incremental steps.
    if (plot_steps){
                
        infercnv_obj_07 <- infercnv_obj
        
        save('infercnv_obj_07', file=file.path(out_path, "07_remove_ref_avg_from_obs.infercnv_obj"))

        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
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

    flog.info("\n\n\tSTEP 08: removing genes at chr ends\n")
        
    infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, window_length)
    
                                        # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_08 <- infercnv_obj
        
        save('infercnv_obj_08', file=file.path(out_path, "08_remove_gene_at_chr_ends.infercnv_obj"))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="08_remove_genes_at_chr_ends",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.08_remove_genes_at_chr_ends",
                 write_expr_matrix=TRUE)
        
    }
    
    
    ################################
    # Step 10: de-noising 

    flog.info("\n\n\tSTEP 10: Denoising\n")
        
    if (! is.na(noise_filter)) {

        if (noise_filter > 0) {
            flog.info(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
            infercnv_obj <- clear_noise(infercnv_obj,
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
        flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", sd_amplifier))
        infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj,
                                                    sd_amplifier = sd_amplifier,
                                                    adjust_towards_zero=FALSE)
    }
    
    if (plot_steps){
        
        infercnv_obj_10 <- infercnv_obj
        
        save('infercnv_obj_10', file=file.path(out_path, "10_denoise.infercnv_obj"))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="10_denoised",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.10_denoised")
        
    }

    ##################################
    # STEP 10: Remove outliers for viz

    flog.info("\n\n\tSTEP 11: Removing outliers\n")
        
    infercnv_obj = remove_outliers_norm(infercnv_obj,
                                        out_method=method_bound_vis,
                                        lower_bound=lower_bound_vis,
                                        upper_bound=upper_bound_vis)
    
        
    # Plot incremental steps.
    if (TRUE) { #plot_steps){


        infercnv_obj_11 <- infercnv_obj

        save('infercnv_obj_11', file=file.path(out_path, "11_remove_outlier.infercnv_obj"))

        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="11_removed_outliers",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.11_removed_outliers")
    }
    
    return(infercnv_obj)
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
    
    flog.info("Creating NGCHM as infercnv.ngchm")
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
subtract_ref_expr_from_obs <- function(infercnv_obj,
                                       method="by_mean"
                                       ) {
    
                                        # r = genes, c = cells
    flog.info(paste("::subtract_ref_expr_from_obs:Start", sep=""))




    # Max and min mean gene expression within reference groups.
    average_max <- NULL
    average_min <- NULL
    # average_reference_obs <- average_data[,ref_observations, drop=FALSE]
    # Reference gene within reference groups
    # now reference indexes of ref_groups are relative to the full average_data matrix and not the average_reference_obs references submatrix
    
    #infercnv_obj <- invert_log2xplus1(infercnv_obj) 

    
    ref_groups = infercnv_obj@reference_grouped_cell_indices

    for (ref_group in ref_groups) {
        
        if (method == "by_mean") {

            grp_average <- rowMeans(infercnv_obj@processed.data[ , ref_group, drop=FALSE], na.rm=TRUE)
            if(is.null(average_max)){
                average_max <- grp_average
            }
            if(is.null(average_min)){
                average_min <- grp_average
            }
            average_max <- pmax(average_max, grp_average)
            average_min <- pmin(average_min, grp_average)

        } else if (method == "by_quantiles") {

            grp_expression_data = infercnv_obj@processed.data[, ref_group, drop=FALSE, na.rm=TRUE]
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
    # TODO:  can we vectorize this?
    #     and if not, set up a progress bar?
    for(gene_i in 1:nrow(infercnv_obj@processed.data)){
        current_col <- infercnv_obj@processed.data[gene_i, ]

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
        infercnv_obj@processed.data[gene_i, ] <- row_init
    }

   # infercnv_obj <- log2xplus1(infercnv_obj)    
    
    return(infercnv_obj)

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
    flog.info(paste("::split_references:Start", sep=""))
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

remove_outliers_norm <- function(infercnv_obj,
                                 out_method="average_bound",
                                 lower_bound=NA,
                                 upper_bound=NA) {

    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))
    
    data <- infercnv_obj@processed.data
    
    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        logging::logerror("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }
    
    
    
    if ( (! is.na(lower_bound)) & (! is.na(upper_bound)) ) {
        
        ## using the values as specified.
        flog.info(paste("::remove_outlier_norm: using hard thresholds: ",
                        "lower_bound:" , lower_bound,
                        "upper_bound:", upper_bound) )
        
    } else if (! is.na(out_method)) {
        
        # using out_method instead of specified bounds.
        flog.info(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))
                
        if (out_method == "average_bound"){

            bounds = get_average_bounds(infercnv_obj)
            lower_bound = bounds[1]
            upper_bound = bounds[2]

        } else {
            flog.error(paste("::remove_outlier_norm:Error, please",
                                    "provide an approved method for outlier",
                                    "removal for visualization."))
            stop(991)
        }
    } else {
        flog.error("::remove_outlier_norm:Error, must specify outmethod or define exact bounds")
        stop(992)
    }
    
    # apply bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound
    
    infercnv_obj@processed.data <- data
    
    return(infercnv_obj)
    
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
center_cell_expr_across_chromosome <- function(infercnv_obj, method="mean") { # or median

    flog.info(paste("::center_smooth across chromosomes per cell"))

    # Center within columns (cells)
    if (method == "median") {
        row_median <- apply(infercnv_obj@processed.data, 2, function(x) { median(x, na.rm=T) } )
        
        infercnv_obj@processed.data <- t(apply(infercnv_obj@processed.data, 1, "-", row_median))
    }
    else {
        # by mean
        row_means <- apply(infercnv_obj@processed.data, 2, function(x) { mean(x, na.rm=T) } )
        
        infercnv_obj@processed.data <- t(apply(infercnv_obj@processed.data, 1, "-", row_means))
    }
    return(infercnv_obj)
}


# Return the indices of the rows that average above the cut off
#
# Args:
# data Data to measure the average row and evaluate
#                 against the cutoff. Row = Genes, Col = Cells.
# cutoff Threshold to be above to be kept.
                                        #
                                        # assumes infercnv_obj@processed.data are log2(x+1) transformed.
                                        # considers the log2(mean(inv_log(data))) < threshold to be removed.

# Returns:
# Returns a vector of row indices to keep (are above the cutoff).
require_above_min_mean_expr_cutoff <- function(infercnv_obj, min_mean_expr_cutoff) {

    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))

    # restrict to reference cells:
    ref_cells_data <- infercnv_obj@processed.data[ , get_reference_grouped_cell_indices(infercnv_obj) ]
    
    average_gene <- rowMeans(ref_cells_data)

    flog.info(paste("::process_data:Averages (counts).", sep=""))
    # Find averages above a certain threshold
    indices <- which(average_gene < min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(paste("Removing ", length(indices),
                        " genes from matrix as below mean expr threshold: ",
                        min_mean_expr_cutoff), sep="")

        infercnv_obj@processed.data = infercnv_obj@processed.data[ -1 * indices, ]

                                        # match w/ gene_order info
        infercnv_obj@gene_order = infercnv_obj@gene_order[ -1 * indices, ]

        validate_infercnv_obj(infercnv_obj)
    }
    
    return(infercnv_obj)
        
}


#' indicate which genes (rows) have at least specified min_cells_per_gene
#'
#' Args
#' @param data Data (expression) matrix
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#' @param obs_idx vector containing the column indices for the observed (tumor) cells
#' @param ref_idx vector containing the column indices for teh reference (normal) cells

require_above_min_cells_ref <- function(infercnv_obj, min_cells_per_gene) {

    ref_cell_indices = get_reference_grouped_cell_indices(infercnv_obj)
    
    ref_data = infercnv_obj@processed.data[,ref_cell_indices]
    
    ref_genes_passed = which(apply(ref_data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))

    num_genes_total = length(ref_cell_indices)
    num_removed = num_genes_total - length(ref_genes_passed)
    if (num_removed > 0) {
        flog.info(paste("Removed ", num_removed, " genes having fewer than ",
                        min_cells_per_gene, " min cells per gene. = ",
                        num_removed / num_genes_total * 100, " % genes removed here."), sep="")

        if (num_removed == num_genes_total) {

            flog.warn(paste("::All genes removed! Must revisit your data..., cannot continue here."))
            stop(998)
        }
        
        infercnv_obj@processed.data <- infercnv_obj@processed.data[ref_genes_passed, ]

                                        # match w/ gene_order
        infercnv_obj@gene_order <- infercnv_obj@gene_order[ref_genes_passed, ]

        validate_infercnv_obj(infercnv_obj)
        
    }
    else {

        flog.info("no genes removed due to min cells/gene filter")
        
    }
    
    return(infercnv_obj)
    
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
    flog.info(paste("::order_reduce:Start.", sep=""))
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
        flog.info(paste("::process_data:order_reduce:The position file ",
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
    flog.info(paste("::process_data:order_reduce:Reduction from positional ",
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
clear_noise <- function(infercnv_obj, threshold, adjust_towards_zero) {
    
    flog.info(paste("********* ::clear_noise:Start. threshold: ", threshold,  " adj_towards_zero: ", adjust_towards_zero, sep=""))

    if (threshold == 0) {
        return(infercnv_obj); # nothing to do
    }

    smooth_matrix = infercnv_obj@processed.data
    
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

    infercnv_obj@processed.data <- smooth_matrix
    
    return(infercnv_obj)
}


# clear_noise_via_ref_quantiles: define noise levels based on quantiles within the ref (normal cell) distribution.
# Any data points within this defined quantile are set to zero.

clear_noise_via_ref_mean_sd <- function(infercnv_obj, sd_amplifier=1.5, adjust_towards_zero=FALSE) {

    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    vals = infercnv_obj@processed.data[,ref_idx]
    
    #vals[vals==0] = NA  # use remaining ref vals that weren't already turned to zeros
    
    upper_bound <- mean(apply(vals, 2, function(x) sd(x, na.rm=T))) * sd_amplifier
    
    lower_bound <- -1 * upper_bound
        
    flog.info(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                           lower_bound, "-", upper_bound, sep=" "))


    smooth_matrix <- infercnv_obj@processed.data
    
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
    
    infercnv_obj@processed.data <- smooth_matrix
    
    return(infercnv_obj)
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
.remove_tails <- function(smooth_matrix, chr, tail_length){

    #flog.info(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
         tail_length <- floor(chr_length / 3)
    }
    remove_indices <- chr[1:tail_length]
    remove_indices <- c(remove_indices,
                        chr[ ( (chr_length + 1) - tail_length):
                             chr_length])
    
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
smooth_by_chromosome <- function(infercnv_obj, window_length, smooth_ends=TRUE){

    data = infercnv_obj@processed.data
    
    flog.info(paste("::smooth_window:Start.", sep=""))
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(infercnv_obj)
    }
    if (window_length > nrow(data)){
        flog.warn("window length exceeds number of rows in data, returning original unmodified data")
        return(infercnv_obj)
    }
    
    tail_length <- (window_length - 1) / 2
    num_genes <- nrow(data)
    data_sm <- apply(data,
                     2,
                     .smooth_window_helper,
                     window_length=window_length)
    logging::logdebug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))

    if (smooth_ends) {
        # Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
        data_sm <- apply(data_sm,
                         2,
                         .smooth_ends_helper,
                         tail_length=tail_length)

    }
    
    # Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)

    infercnv_obj@processed.data <- data_sm

    return(infercnv_obj)
}

# Helper function for smoothing the ends of a moving average.
#
# Args:
# obs_data: Data to smooth
# tail_length:  Length of the tail to smooth.
#
# Returns:
# Data smoothed.
.smooth_ends_helper <- function(obs_data, tail_length) {

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
.smooth_window_helper <- function(obs_data, window_length){

    nas = is.na(obs_data)
    vals = obs_data[! nas]

    smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)

    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]

    obs_data[! nas] = vals

    return(obs_data)
}




get_average_bounds <- function (infercnv_obj) {

    lower_bound <- mean(apply(infercnv_obj@processed.data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(infercnv_obj@processed.data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))
    
    return(c(lower_bound, upper_bound))

}


log2xplus1 <- function(invercnv_obj) {

    infercnv_obj@processed.data <- log2(infercnv_obj@raw.data + 1)

    return(infercnv_obj)
        
}

invert_log2xplus1 <- function(infercnv_obj) {

    infercnv_obj@processed.data <- 2^infercnv_obj@processed.data - 1

    return(infercnv_obj)
}


make_zero_NA <- function(infercnv_obj) {

    infercnv_obj@processed.data <- infercnv_obj@processed.data[infercnv_obj@processed.data == 0] <- NA
    
    return(infercnv_obj)

}

transform_to_reference_based_Zscores <- function(infercnv_obj) {

    # center and convert to z-scores
    flog.info(paste("::center_and_Zscore_conversion", sep=""))
    
    # remember, genes are rows, cells are cols
    
    # centering and z-scores based on the reference (normal) cells:
    
    # ref data represent the null distribution
    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    
    ref_data = infercnv_obj@processed.data[,ref_idx]
    
    gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=T)})
    gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=T)})

    # assume at least Poisson level variation
    gene_ref_sd = pmax(gene_ref_sd, gene_ref_mean)
    
    # center all genes at the ref (normal) center:
    infercnv_obj@processed.data = sweep(infercnv_obj@processed.data, 1, gene_ref_mean, FUN="-")
    
    # convert to z-scores based on the ref (normal) distribution
    infercnv_obj@processed.data = sweep(infercnv_obj@processed.data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.
    
    
    return(infercnv_obj)
    
}


mean_center_gene_expr <- function(infercnv_obj) {

    flog.info(paste("::centering", sep=""))

    infercnv_obj@processed.data <- sweep(infercnv_obj@processed.data, 1, rowMeans(infercnv_obj@processed.data, na.rm=T), FUN="-")
        
    return(infercnv_obj)
}


get_reference_grouped_cell_indices <- function(infercnv_obj) {

    return( unlist(infercnv_obj@reference_grouped_cell_indices) )

}

apply_max_threshold_bounds <- function(infercnv_obj, threshold) {

    flog.info(paste("::process_data:setting max centered expr, threshold set to: +/-: ", threshold))

    infercnv_obj@processed.data[infercnv_obj@processed.data > threshold] <- threshold
    infercnv_obj@processed.data[infercnv_obj@processed.data < (-1 * threshold)] <- -1 * threshold

    return(infercnv_obj)
}


remove_genes_at_ends_of_chromosomes <- function(infercnv_obj, window_length) {

    contig_tail= (window_length - 1) / 2
    
    remove_indices <- c()
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    for (chr in chrs){
        #flog.info(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- .remove_tails(infercnv_obj@processed.data,
                                    which(gene_chr_listing == chr),
                                    contig_tail)
        
        #logging::logdebug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))

        remove_indices <- c(remove_indices, remove_chr)

    }
    if (length(remove_indices) > 0){

        infercnv_obj@processed.data <- infercnv_obj@processed.data[ -1 * remove_indices, ]
        infercnv_obj@gene_order <- infercnv_obj@gene_order[ -1 * remove_indices, ]

        validate_infercnv_obj(infercnv_obj)

        flog.info(paste("::process_data:Remove genes at chr ends, ",
                        "new dimensions (r,c) = ",
                        paste(dim(infercnv_obj@processed.data), collapse=","),
                        " Total=", sum(infercnv_obj@processed.data, na.rm=TRUE),
                        " Min=", min(infercnv_obj@processed.data, na.rm=TRUE),
                        " Max=", max(infercnv_obj@processed.data, na.rm=TRUE),
                        ".", sep=""))
        
    }
    else {
        flog.error("No genes removed at chr ends.... something wrong here")
        stop(1234)
    }

    validate_infercnv_obj(infercnv_obj)
    
    return(infercnv_obj)

}

    

validate_infercnv_obj <- function(infercnv_obj) {

    flog.info("validating infercnv_obj")
    
    if (all.equal(rownames(infercnv_obj@processed.data), rownames(infercnv_obj@gene_order))) {
        # all good.
        return();
        
    }
    else {

        flog.error("hmm.... rownames(infercnv_obj@processed.data != rownames(infercnv_obj@gene_order))")
        broken.infercnv_obj = infercnv_obj
        save('broken.infercnv_obj', file="broken.infercnv_obj")
        
    }
    
    
    genes = setdiff(rownames(infercnv_obj@processed.data), rownames(infercnv_obj@gene_order))
    if (length(genes) != 0) {
        flog.error(paste("The following genes are in infercnv_obj@processed.data and not @gene_order:", paste(genes, collapse=","),
                         sep=" "))
        
    }

    genes = setdiff(rownames(infercnv_obj@gene_order), rownames(infercnv_obj@processed.data))
    if (length(genes) != 0) {
        flog.error(paste("The following genes are in @gene_order and not infercnv_obj@processed.data:", paste(genes, collapse=","),
                         sep=" "))
        
    }

    stop("Problem detected w/ infercnv_obj")
    
}


reinit_infercnv <- function(infercnv_obj) {
    
    # restore the processed.data with the original raw.data

    infercnv_obj@processed.data <- infercnv_obj@raw.data

    return(infercnv_obj)

}

symmetrical_logxplus1 <- function(infercnv_obj) {

    data = infercnv_obj@processed.data

    zero_val_pos = (data == 0)

    data[data>0] = log2(data[data>0] + 1)

    data[data<0] = -1 * log2(-1 * data[data<0] + 1)

    data[zero_val_pos] = NA

    infercnv_obj@processed.data = data

    return(infercnv_obj)
}

