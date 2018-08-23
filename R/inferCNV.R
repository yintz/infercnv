#!/usr/bin/env Rscript


options(error = function() traceback(2))

CHR <- "chr"
START <- "start"
STOP <- "stop"

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


# Returns the color palette for contigs.
#
# Returns:
# Color Palette
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}


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
#' @param name_ref_groups Names of groups from the "annotations" table whose cells
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
process_data <- function(data,
                      gene_order,
                      cutoff,
                      reference_obs,
                      transform_data,
                      window_length,
                      max_centered_threshold,
                      noise_filter=NA,
                      noise_quantiles=c(0.025, 0.975),
                      name_ref_groups,
                      num_ref_groups,
                      out_path,
                      obs_annotations_groups,
                      obs_annotations_names,
                      grouping_key_coln,
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
        groups_ref <- split_references(average_data=data, #data_smoothed,
                                       ref_obs=reference_obs,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)
    }
    else {
        groups_ref <- name_ref_groups
    }
    logging::loginfo(paste("::process_data:split_reference. ",
                           "found ",length(groups_ref)," reference groups.",
                           sep=""))

    # Return list of matrices
    ret_list <- list()

    ret_list[["REF_GROUPS"]] = groups_ref
    ret_list[["REF_OBS_IDX"]] = reference_obs
    chr_order_for_plotting <- paste(as.vector(as.matrix(gene_order[1])))

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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
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
                     reference_idx=ret_list[["REF_OBS_IDX"]],
                     ref_contig=NULL,
                     contig_cex=1,
                     ref_groups=ret_list[["REF_GROUPS"]],
                     name_ref_groups=name_ref_groups,
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
                     reference_idx=ret_list[["REF_OBS_IDX"]],
                     ref_contig=NULL,
                     contig_cex=1,
                     ref_groups=ret_list[["REF_GROUPS"]],
                     name_ref_groups=name_ref_groups,
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


    if (use_zscores) {

        # center and convert to z-scores
        logging::loginfo(paste("::center_and_Zscore_conversion", sep=""))
        data = t(scale(t(data), center=T, scale=T))
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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="04_center_with_threshold",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.04_center_with_threshold")

    }


    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_length)

    data <- NULL
    logging::loginfo(paste("::process_data:Smoothed data.", sep=""))
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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="05_smoothed",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.05_smoothed")
    }

    # Center cells/observations after smoothing. This helps reduce the
    # effect of complexity.
    data_smoothed <- center_smoothed(data_smoothed)

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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
                 out_dir=out_path,
                 color_safe_pal=FALSE,
                 x.center=0,
                 title="07_remove_average",
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename="infercnv.07_remove_average")

    }


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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
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

    # Clear noise: set values to zero that are in the defined noise range

    if ( (! is.na(noise_filter)) && noise_filter > 0) {
        logging::loginfo(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
        data_smoothed <- clear_noise(smooth_matrix=data_smoothed,
                                      threshold=noise_filter)
    }
    else {
        logging::loginfo(paste("::process_data:Remove noise, noise threshold defined via quantiles: ", noise_quantiles))
        data_smoothed <- clear_noise_via_ref_quantiles(smooth_matrix=data_smoothed,
                                                        ref_idx=unlist(ret_list[["REF_GROUPS"]]),
                                                        quantiles=noise_quantiles)
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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
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
                 reference_idx=ret_list[["REF_OBS_IDX"]],
                 ref_contig=NULL,
                 contig_cex=1,
                 ref_groups=ret_list[["REF_GROUPS"]],
                 name_ref_groups=name_ref_groups,
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

# Not testing, params ok
# Log intermediate step with a plot and text file of the steps.
#
# Args:
# data: The data frame to plot.
# plot_name: The absolute path to the pdf to be plotted.
#
# Returns:
# No return
plot_step <- function(data, plot_name){
    text_file <- unlist(strsplit(plot_name, "\\."))
    text_file <- paste(c(text_file[1:length(text_file)], "txt"),
                       collapse=".", sep=".")
    pdf(plot_name, useDingbats=FALSE)
    image(as.matrix(data),
          col=colorRampPalette(c("cyan","white","white","magenta"))(n=100))
    dev.off()
    write.table(data, file=text_file, quote=F, sep="\t")
}


#' Formats the data and sends it for plotting.
#'
#' @title Plot the matrix as a heatmap. Clustering is on observation only, gene position is preserved.
#'
#' @param plot_data Data matrix to plot (columns are observations).
#' @param contigs The contigs the data is group in in order of rows.
#' @param reference_idx Vector of reference indices.
#' @param ref_contig If given, will focus cluster on only genes in this contig.
#' @param ref_groups Groups of vector indices (as indices in reference_idx).
#' @param out_dir Directory in which to save pdf and other output.
#' @param title Plot title.
#' @param obs_title Title for the observations matrix.
#' @param ref_title Title for the reference matrix.
#' @param obs_annotations_groups Vector of observations annotations group assignation (as indices).
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores k_obs_groups.
#' @param contig_cex Contig text size.
#' @param k_obs_groups Number of groups to break observation into.
#' @param x.center Value on which to center expression.
#' @param hclust_method Clustering method to use for hclust.
#' @param color_safe_pal Logical indication of using a color blindness safe
#'                          palette.
#' @param pdf_filename Filename to save the figure to.
#'
#' @return
#' No return, void.
#'
plot_cnv <- function(plot_data,
                     contigs,
                     reference_idx,
                     ref_contig,
                     ref_groups,
                     name_ref_groups,
                     out_dir,
                     title,
                     obs_title,
                     ref_title,
                     obs_annotations_groups,
                     obs_annotations_names,
                     grouping_key_coln,
                     cluster_by_groups=FALSE,
                     contig_cex=1,
                     k_obs_groups=1,
                     x.center=0,
                     hclust_method='average',
                     color_safe_pal=TRUE,
                     output_filename="infercnv",
                     output_format="png"){

    logging::loginfo(paste("::plot_cnv:Start", sep=""))
    logging::loginfo(paste("::plot_cnv:Current data dimensions (r,c)=",
                            paste(dim(plot_data), collapse=","),
                            " Total=", sum(plot_data, na.rm=T),
                            " Min=", min(plot_data, na.rm=T),
                            " Max=", max(plot_data, na.rm=T),
                            ".", sep=""))
    logging::loginfo(paste("::plot_cnv:Depending on the size of the matrix",
                           " this may take a moment.",
                           sep=""))

    # Contigs
    unique_contigs <- unique(contigs)
    n_contig <- length(unique_contigs)
    ct.colors <- get_group_color_palette()(n_contig)
    names(ct.colors) <- unique_contigs

    # Select color palette
    custom_pal <- color.palette(c("purple3", "white", "darkorange2"),
                                c(2, 2))
    if (color_safe_pal == FALSE){
        custom_pal <- color.palette(c("darkblue", "white", "darkred"),
                                    c(2, 2))
    }

    # Row seperation based on reference
    ref_idx <- NULL
    if (!is.null(reference_idx) && length(reference_idx) < ncol(plot_data)){
        reference_idx <- which(colnames(plot_data) %in% reference_idx)
        if(length(reference_idx) > 0){
            ref_idx <- reference_idx
        }
    }

    # Column seperation by contig and label axes with only one instance of name
    contig_tbl <- table(contigs)[unique_contigs]
    col_sep <- cumsum(contig_tbl)
    col_sep <- col_sep[-1 * length(col_sep)]   ## FIXME:  removing last entry?
    # These labels are axes labels, indicating contigs at the first column only
    # and leaving the rest blank.
    contig_labels <- c()
    contig_names <-c()
    for (contig_name in names(contig_tbl)){
        contig_labels <- c(contig_labels,
                           contig_name,
                           rep("", contig_tbl[contig_name] - 1))
        contig_names <- c(contig_names,rep(contig_name,contig_tbl[contig_name]))
    }

    # Calculate how many rows will be made for the number of columns in the grouping key
    grouping_key_rown <- c()
    grouping_key_rown[1] <- ceiling(length(obs_annotations_names)/grouping_key_coln[1])
    grouping_key_rown[2] <- ceiling(length(name_ref_groups)/grouping_key_coln[2])
    # Calculate how much bigger the output needs to be to accodomate for the grouping key
    grouping_key_height <- c((grouping_key_rown[2] + 2) * 0.175, (grouping_key_rown[1] + 3) * 0.175)

    # Rows observations, Columns CHR
    if (output_format == "pdf") {
        pdf(paste(out_dir, paste(output_filename, ".pdf", sep=""), sep="/"),
            useDingbats=FALSE,
            width=10,
            height=(8.13 + sum(grouping_key_height)),
            paper="USr")
    }
    else if (output_format == "png") {
        png(paste(out_dir, paste(output_filename, ".png", sep=""), sep="/"),
            width=10,
            height=(8.13 + sum(grouping_key_height)),
            units="in",
            res=600)
    }

    # Plot observations
    ## Make Observation Samples
    ## Remove observation col names, too many to plot
    ## Will try and keep the reference names
    ## They are more informative anyway
    obs_data <- plot_data
    if (!is.null(ref_idx)){
        obs_data <- plot_data[, -1 * ref_idx, drop=FALSE]
        if (ncol(obs_data) == 1){
                plot_data <- cbind(obs_data, obs_data)
                names(obs_data) <- c("", names(obs_data)[1])
        }
    }

    obs_data_t <- t(obs_data)

    # Subsample the data to only the references and update the ref_group indexes to match their new indexes
    # ref_data_t <- plot_data[, ref_idx, drop=FALSE]
    ref_data_t <- NULL
    updated_ref_groups <- list()
    current_ref_count <- 1
    current_grp_idx <- 1
    plot_data <- as.matrix(plot_data)
    for (ref_grp in ref_groups) {
        ref_data_t <- cbind(ref_data_t, plot_data[, ref_grp, drop=FALSE])
        updated_ref_groups[[current_grp_idx]] = seq(current_ref_count, current_ref_count + length(ref_grp) - 1)
        current_ref_count <- current_ref_count + length(ref_grp)
        current_grp_idx <- current_grp_idx + 1
    }
    ref_groups <- updated_ref_groups

    nb_breaks <- 16
    breaksList_t <-
        seq(min(min(obs_data_t, na.rm=TRUE), min(ref_data_t, na.rm=TRUE)),
        max(max(obs_data_t,na.rm=TRUE), max(ref_data_t, na.rm=TRUE)),
        length.out=nb_breaks)


    # Create file base for plotting output
    force_layout <- plot_observations_layout(grouping_key_height=grouping_key_height)
    plot_cnv_observations(obs_data=obs_data_t,
                          file_base_name=out_dir,
                          cluster_contig=ref_contig,
                          contig_colors=ct.colors[contigs],
                          contig_labels=contig_labels,
                          contig_names=contig_names,
                          col_pal=custom_pal,
                          contig_seps=col_sep,
                          num_obs_groups=k_obs_groups,
                          obs_annotations_groups=obs_annotations_groups,
                          obs_annotations_names=obs_annotations_names,
                          grouping_key_coln=grouping_key_coln[1],
                          cluster_by_groups=cluster_by_groups,
                          cnv_title=title,
                          cnv_obs_title=obs_title,
                          contig_lab_size=contig_cex,
                          breaksList=breaksList_t,
                          x.center=x.center,
                          hclust_method=hclust_method,
                          layout_lmat=force_layout[["lmat"]],
                          layout_lhei=force_layout[["lhei"]],
                          layout_lwid=force_layout[["lwid"]])
    obs_data <- NULL

    if(!is.null(ref_idx)){
        plot_cnv_references(ref_data=ref_data_t,
                            ref_groups=ref_groups,
                            name_ref_groups=name_ref_groups,
                            grouping_key_coln=grouping_key_coln[2],
                            col_pal=custom_pal,
                            contig_seps=col_sep,
                            file_base_name=out_dir,
                            cnv_ref_title=ref_title,
                            breaksList=breaksList_t,
                            x.center=x.center,
                            layout_add=TRUE)
    }
    dev.off()
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the observational samples
#
# Args:
#' @param obs_data Data to plot as observations. Rows = Cells, Col = Genes.
#' @param col_pal The color palette to use.
#' @param contig_colors The colors for the contig bar.
#' @param contig_labels The labels for the contigs.
#' @param contig_names Names of the contigs.
#' @param contig_seps Indices for line seperators of contigs.
#' @param num_obs_groups Number of groups of observations to create.
#' @param file_base_name Base of the file to used to make output file names.
#' @param cnv_title Title of the plot.
#' @param cnv_obs_title Title for the observation matrix.
#' @param contig_lab_size Text size for contigs.
#' @param cluster_contig A value directs cluster to only genes on this contig.
#' @param obs_annotations_groups Vector of observations annotations group assignation (as indices).
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores num_obs_groups.
#' @param breaksList List of values used as splitters on coloring range.
#' @param hclust_method Method to use for hclust.
#' @param testing If TRUE, does not plot anything.
#' @param layout_lmat lmat values to use in layout.
#' @param layout_lhei lhei values to use in layout.
#' @param layout_lwid lwid values to use in layout.
#
#' @return Void.
# Returns:
# Void
plot_cnv_observations <- function(obs_data,
                                  col_pal,
                                  contig_colors,
                                  contig_labels,
                                  contig_names,
                                  contig_seps,
                                  num_obs_groups,
                                  file_base_name,
                                  cnv_title,
                                  cnv_obs_title,
                                  contig_lab_size=1,
                                  cluster_contig=NULL,
                                  obs_annotations_groups,
                                  obs_annotations_names,
                                  grouping_key_coln,
                                  cluster_by_groups,
                                  breaksList,
                                  x.center,
                                  hclust_method="average",
                                  testing=FALSE,
                                  layout_lmat=NULL,
                                  layout_lhei=NULL,
                                  layout_lwid=NULL){

    logging::loginfo("plot_cnv_observation:Start")
    logging::loginfo(paste("Observation data size: Cells=",
                           nrow(obs_data),
                           "Genes=",
                           ncol(obs_data),
                           sep=" "))
    observation_file_base <- paste(file_base_name, "observations.txt", sep=.Platform$file.sep)

    # Output dendrogram representation as Newick
    # Need to precompute the dendrogram so we can manipulate
    # it before the heatmap plot
    ## Optionally cluster by a specific contig
    hcl_desc <- "General"
    hcl_group_indices <- 1:ncol(obs_data)
    if(!is.null(cluster_contig)){
        hcl_contig_indices <- which(contig_names == cluster_contig)
        if(length(hcl_group_indices)>0){
            hcl_group_indices <- hcl_contig_indices
            hcl_desc <- cluster_contig
            logging::loginfo(paste("plot_cnv_observation:Clustering only by contig ", cluster_contig))
        } else {
           logging::logwarn(paste("plot_cnv_observations: Not able to cluster by",
                                     cluster_contig,
                                     "Clustering by all genomic locations.",
                                     "To cluster by local genomic location next time",
                                     "select from:",
                                     unique(contig_names),
                                     collapse=",",
                                     sep=" "))
        }
    }
                                        # HCL with a inversely weighted euclidean distance.
    logging::loginfo(paste("clustering observations via method: ", hclust_method, sep=""))
    # obs_hcl <- NULL
    obs_dendrogram <- list()
    ordered_names <- NULL
    isfirst <- TRUE
    hcl_obs_annotations_groups <- vector()
    obs_seps <- c()
    if (cluster_by_groups) {
        for (i in seq(1, max(obs_annotations_groups))) {
            group_obs_hcl <- hclust(dist(obs_data[which(obs_annotations_groups == i), hcl_group_indices]), method=hclust_method)
            ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups == i), hcl_group_indices])[group_obs_hcl$order])
            if (isfirst) {
                write.tree(as.phylo(group_obs_hcl),
                   file=paste(file_base_name, "observations_dendrogram.txt", sep=.Platform$file.sep))
                isfirst <- FALSE
            }
            else {
                write.tree(as.phylo(group_obs_hcl),
                   file=paste(file_base_name, "observations_dendrogram.txt", sep=.Platform$file.sep), append=TRUE)
            }
            group_obs_dend <- as.dendrogram(group_obs_hcl)
            obs_dendrogram[[length(obs_dendrogram) + 1]] <- group_obs_dend
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, rep(i, length(which(obs_annotations_groups == i))))
            obs_seps <- c(obs_seps, length(ordered_names))
        }
        if (length(obs_dendrogram) > 1) {
            obs_dendrogram <- do.call(merge, obs_dendrogram)
        } else {
            obs_dendrogram <- obs_dendrogram[[1]]
        }
        split_groups <- rep(1, dim(obs_data)[1])
        names(split_groups) <- ordered_names
    }
    else {
        obs_hcl <- hclust(dist(obs_data[,hcl_group_indices]), method=hclust_method)
        write.tree(as.phylo(obs_hcl),
                   file=paste(file_base_name, "observations_dendrogram.txt", sep=.Platform$file.sep))
        obs_dendrogram <- as.dendrogram(obs_hcl)
        ordered_names <- row.names(obs_data)[obs_hcl$order]
        split_groups <- cutree(obs_hcl, k=num_obs_groups)
        split_groups <- split_groups[ordered_names]
        hcl_obs_annotations_groups <- obs_annotations_groups[obs_hcl$order]

        # Make a file of members of each group
        logging::loginfo("plot_cnv_observation:Writing observations by grouping.")
        for (cut_group in unique(split_groups)){
          group_memb <- names(split_groups)[which(split_groups == cut_group)]
          # Write group to file
          memb_file <- file(paste(file_base_name,
                                  paste(hcl_desc,"HCL",cut_group,"members.txt",sep="_"),
                                  sep=.Platform$file.sep))
          write.table(obs_data[group_memb,], memb_file)
          # Record seperation
          ordered_memb <- which(ordered_names %in% group_memb)
          if (is.null(obs_seps)) {
            obs_seps <- c(length(ordered_memb))
          }
          else {
            obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + length(ordered_memb)))
          }
        }
        obs_seps <- c(obs_seps, length(ordered_names))
    }

    if (length(obs_seps) > 1) {
        obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) - 1):1]
    }

    # Output HCL group membership.
    # Record locations of seperations

    # Make colors based on groupings
    row_groupings <- get_group_color_palette()(length(table(split_groups)))[split_groups]
    row_groupings <- cbind(row_groupings, get_group_color_palette()(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])
    annotations_legend <- cbind(obs_annotations_names, get_group_color_palette()(length(table(hcl_obs_annotations_groups))))

    # Make a file of coloring and groupings
    logging::loginfo("plot_cnv_observation:Writing observation groupings/color.")
    groups_file_name <- file.path(file_base_name, "observation_groupings.txt")
    # file_groups <- rbind(split_groups,row_groupings)
    file_groups <- cbind(split_groups, row_groupings[,1], hcl_obs_annotations_groups, row_groupings[,2])
    #row.names(file_groups) <- c("Group","Color")
    colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color", "Annotation Group", "Annotation Color")
    # write.table(t(file_groups), groups_file_name)
    write.table(file_groups, groups_file_name)



    # obs_seps <- unique(obs_seps)
    # obs_seps <- sort(obs_seps)

    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(obs_data),
                                     col_count=ncol(obs_data),
                                     row_seps=obs_seps,
                                     col_seps=contig_seps)

    obs_data <- obs_data[ordered_names, ]

    # Remove row/col labels, too cluttered
    # and print.
    orig_row_names <- row.names(obs_data)
    row.names(obs_data) <- rep("", nrow(obs_data))

    data_observations <- heatmap.cnv(obs_data,
                                        Rowv=obs_dendrogram,
                                        Colv=FALSE,
                                        cluster.by.row=TRUE,
                                        cluster.by.col=FALSE,
                                        main=cnv_title,
                                        ylab=cnv_obs_title,
                                        margin.for.labCol=2,
                                        xlab="Genomic Region",
                                        key=TRUE,
                                        labCol=contig_labels,
                                        cexCol=contig_lab_size,
                                        notecol="black",
                                        density.info="histogram",
                                        denscol="blue",
                                        trace="none",
                                        dendrogram="row",
                                        cexRow=0.8,
                                        breaks=breaksList,
                                        scale="none",
                                        x.center=x.center,
                                        color.FUN=col_pal,
                                        if.plot=!testing,
                                        # Seperate by contigs
                                        sepList=contigSepList,
                                        sep.color=c("black","black"),
                                        sep.lty=1,
                                        sep.lwd=1,
                                        # Color rows by user defined cut
                                        RowIndividualColors=row_groupings,
                                        annotations_legend=annotations_legend,
                                        grouping_key_coln=grouping_key_coln,
                                        # Color by contigs
                                        ColIndividualColors=contig_colors,
                                        # Legend
                                        key.title="Distribution of Expression",
                                        key.xlab="Modified Expression",
                                        key.ylab="Count",
                                        # Layout
                                        force_lmat=layout_lmat,
                                        force_lwid=layout_lwid,
                                        force_lhei=layout_lhei)

    # Write data to file.
    logging::loginfo(paste("plot_cnv_references:Writing observation data to",
                           observation_file_base,
                           sep=" "))
    row.names(obs_data) <- orig_row_names
    write.table(obs_data[data_observations$rowInd,data_observations$colInd],
                file=observation_file_base)
}

# Not Testing, params ok.
# Create the layout for the plot
# This is a modification of the original
# layout from the GMD heatmap.3 function
#
# Returns:
# list with slots "lmat" (layout matrix),
#                             "lhei" (height, numerix vector),
#                             and "lwid" (widths, numeric vector)
plot_observations_layout_original <- function()
{
    ## Plot observational samples
    obs_lmat <- c(0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                  7, 9, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                  0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  5, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    obs_lmat <- matrix(obs_lmat,ncol=14,byrow=TRUE)

    obs_lhei <- c(1.125, 2.215, .15,
                   .5, .5, .5,
                   .5, .5, .5,
                   .5, .5, .5,
                  0.0075, 0.0075, 0.0075)

    return(list(lmat=obs_lmat,
           lhei=obs_lhei,
           lwid=NULL))
}

plot_observations_layout <- function(grouping_key_height)
{
    ## Plot observational samples
    obs_lmat <- c(0,  0,  0,  0,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
                  # 8,  0, 10,  0,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                  8, 11, 10,  0,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                  0,  0,  0,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
    obs_lmat <- matrix(obs_lmat,ncol=14,byrow=TRUE)

    obs_lhei <- c(1.125, 2.215, .15,
                   .5, .5, .5,
                   .5, .5, .5,
                   .5, .5, .5,
                  0.03, grouping_key_height[1]+0.04, grouping_key_height[2]+0.04, 0.03)

    return(list(lmat=obs_lmat,
           lhei=obs_lhei,
           lwid=NULL))
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the reference samples
#
# Args:
# ref_data Data to plot as references. Rows = Cells, Col = Genes
# ref_groups Groups of references to plot together.
# col_pal The color palette to use.
# contig_seps Indices for line seperators of contigs.
# file_base_name Base of the file to used to make output file names.
# cnv_ref_title Title for reference matrix.
# layout_lmat lmat values to use in the layout.
# layout_lwid lwid values to use in the layout.
# layout_lhei lhei values to use in the layout.
# layout_add Indicates the ref image shoudl be added to the previous plot.
# testing: Turns off plotting when true.
#
# Returns:
# Void
plot_cnv_references <- function(ref_data,
                                ref_groups,
                                name_ref_groups,
                                grouping_key_coln,
                                col_pal,
                                contig_seps,
                                file_base_name,
                                cnv_ref_title,
                                breaksList,
                                x.center=x.center,
                                layout_lmat=NULL,
                                layout_lwid=NULL,
                                layout_lhei=NULL,
                                layout_add=FALSE,
                                testing=FALSE){

    logging::loginfo("plot_cnv_references:Start")
    logging::loginfo(paste("Reference data size: Cells=",
                           ncol(ref_data),
                           "Genes=",
                           nrow(ref_data),
                           sep=" "))
    number_references <- ncol(ref_data)
    reference_ylab <- NA
    reference_data_file <- paste(file_base_name, "references.txt", sep=.Platform$file.sep)

    ref_seps <- c()
    # Handle only one reference
    # heatmap3 requires a 2 x 2 matrix, so with one reference
    # I just duplicate the row and hid the second name so it
    # visually looks like it is just taking up the full realestate.
    if(number_references == 1){
        ref_data <- cbind(ref_data, ref_data)
        names(ref_data) <- c("",names(ref_data)[1])
        colnames(ref_data) <- c("", colnames(ref_data)[1])
    }

    # Handle reference groups
    # If there is more than one reference group, visually break
    # up the groups with a row seperator. Also plot the rows in
    # order so the current groups are shown and seperated.
    if(length(ref_groups) > 1){
        i_cur_idx <- 0
        order_idx <- c()
        grouping <- 0
        for (ref_grp in ref_groups){
            i_cur_idx <- i_cur_idx + length(ref_grp)
            ref_seps <- c(ref_seps, i_cur_idx)
            order_idx <- c(order_idx, ref_grp)
        }
        ref_seps <- ref_seps[1:(length(ref_seps) - 1)]
        ref_data <- ref_data[, order_idx, drop=FALSE]
    }

    split_groups <- c()
    i_cur_idx <- 1
    for (ref_grp in ref_groups){
        split_groups <- c(split_groups, rep(i_cur_idx, length(ref_grp)))
        i_cur_idx <- i_cur_idx + 1
    }

    # Make row color column colors from groupings
    logging::loginfo(paste("plot_cnv_references:Number reference groups=",
                           length(ref_groups)),
                           sep=" ")

    # Transpose data.
    ref_data <- t(ref_data)
    # Remove labels if too many.
    ref_orig_names <- row.names(ref_data)
    if (number_references > 20){
        # The reference labs can become clustered
        # Dynamically change labels given a certain number of labels.
        reference_ylab <- cnv_ref_title
        row.names(ref_data) <- rep("", number_references)
    }


    row_groupings <- as.matrix(get_group_color_palette()(length(table(split_groups)))[split_groups])
    annotations_legend <- cbind(name_ref_groups, get_group_color_palette()(length(name_ref_groups)))

    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(ref_data),
                                     col_count=ncol(ref_data),
                                     row_seps=ref_seps,
                                     col_seps=contig_seps)

    # Print controls
    logging::loginfo("plot_cnv_references:Plotting heatmap.")
    data_references <- heatmap.cnv(ref_data,
                                   main=NULL, #NA,
                                   ylab=reference_ylab,
                                   xlab=NULL, #NA,
                                   key=FALSE,
                                   labCol=rep("", nrow(ref_data)),
                                   notecol="black",
                                   trace="none",
                                   dendrogram="none",
                                   Colv=FALSE,
                                   Rowv=FALSE,
                                   cexRow=0.4,
                                   breaks=breaksList,
                                   scale="none",
                                   x.center=x.center,
                                   color.FUN=col_pal,
                                   sepList=contigSepList,
                                   # Row colors and legend
                                   RowIndividualColors=row_groupings,
                                   annotations_legend=annotations_legend,
                                   grouping_key_coln=grouping_key_coln,
                                   # Seperate by contigs
                                   sep.color=c("black","black"),
                                   sep.lty=1,
                                   sep.lwd=1,
                                   if.plot=!testing,
                                   # Layout
                                   force_lmat=layout_lmat,
                                   force_lwid=layout_lwid,
                                   force_lhei=layout_lhei,
                                   force_add=layout_add)

    # Write data to file
    row.names(ref_data) <- ref_orig_names
    logging::loginfo(paste("plot_cnv_references:Writing reference data to",
                           reference_data_file,
                           sep=" "))
    write.table(t(ref_data[data_references$rowInd,data_references$colInd]),
                file=reference_data_file)
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
clear_noise <- function(smooth_matrix, threshold){

    logging::loginfo(paste("********* ::clear_noise:Start. threshold: ", threshold,  sep=""))
    if (threshold > 0){
        smooth_matrix[abs(smooth_matrix) < threshold] <- 0
    }
    return(smooth_matrix)
}

# clear_noise_via_ref_quantiles: define noise levels based on quantiles within the ref (normal cell) distribution.
# Any data points within this defined quantile are set to zero.

clear_noise_via_ref_quantiles <- function(smooth_matrix, ref_idx, quantiles=c(0.025, 0.975) ) {

    vals = smooth_matrix[,ref_idx]

    vals[vals==0] = NA  # use remaining ref vals that weren't already turned to zeros

    lower_quantile = quantiles[1]
    upper_quantile = quantiles[2]

    logging::loginfo(paste("::clear_noise_via_ref_quantiles: using noise quantiles set at: ",
                           lower_quantile, "-", upper_quantile, sep=""))

    lower_bound <- mean(apply(vals, 2,
                              function(x) quantile(x, probs=lower_quantile, na.rm=TRUE)))

    upper_bound <- mean(apply(vals, 2,
                              function(x) quantile(x, probs=upper_quantile, na.rm=TRUE)))

    logging::loginfo(paste("::clear_noise_via_ref_quantiles: removing noise between bounds: ",
                           lower_bound, "-", upper_bound, sep=" "))

    smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = 0

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
smooth_window <- function(data, window_length, smooth_ends=TRUE, re_center=TRUE){

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
    if (re_center) {

        # re-center genes now after the smoothing:
        data_sm = sweep(data_sm, 1, rowMeans(data_sm, na.rm=TRUE))
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


#########################################
# Plotting from GMD
###################
###
# This package relies on heatmap.3 from the library GMD
# This code is a modified version of heatmap.3, some changes
# were required for our plotting.
# The heatmap.cnv function should be considered a modification
# of th GMD library function heatmap.3, all credit goes to
# their authors.
## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## A copy of gtools::invalid
##
## see \code{invalid} in package:gtools for details
## Test if a value is missing, empty, or contains only NA or NULL values
## param: x value to be tested
.invalid <-
  function(x)
{
  if (missing(x) || is.null(x) || length(x) == 0)
    return(TRUE)
  if (is.list(x))
    return(all(sapply(x, .invalid)))
  else if (is.vector(x))
    return(all(is.na(x)))
  else return(FALSE)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
.is.grouped <-
  function(x)
{
  x <- as.character(x)
  x.f <- factor(x,levels=unique(x),ordered=TRUE)
  identical(as.character(sort(x.f)),x)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Call a function with arguments
##
## Call a function with arguments
## FUN function or function name
## ... unnameed function arguments
## MoreArgs named (or unnameed) function arguments
.call.FUN <-
  function(FUN,...,MoreArgs)
{
  FUN <- match.fun(FUN)
  tmp.MoreArgs <- list(...)
  if (!.invalid(MoreArgs)){
    if (length(MoreArgs)>=1){
      for (i in 1:length(MoreArgs)) tmp.MoreArgs[[names(MoreArgs)[i]]] <- MoreArgs[[i]]
    }
  }
  ret <- do.call(FUN, tmp.MoreArgs)
  if ("call" %in% names(ret)){
    ret$call <- match.call()
  }
  if ("call" %in% names(attributes(ret))){
    attr(ret,"call") <- match.call()
  }
  return(ret)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Scale values to make them follow Standard Normal Distribution
##
## Scale values to make them follow Standard Normal Distribution
## param x numeric
## param scale character, indicating the type to scale.
## param na.rm logical
## return an object with the same dimention of `x'.
.scale.data <-
  function(x,scale,na.rm=TRUE)
{
  if(scale=="row"){
    x <- sweep(x,1,rowMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,1,sd,na.rm=na.rm)
    x <- sweep(x,1,sx,FUN="/")
  } else if(scale=="column"){
    x <- sweep(x,2,colMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,2,sd,na.rm=na.rm)
    x <- sweep(x,2,sx,,FUN="/")
  }
  x
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Scale values to a new range: c(low, high)
## x numeric
## low numeric, lower bound of target values
## high numeric, higher bound of target values
## return an object with the same dimention of `x'.
.scale.x <-
  function(x,low=0,high=1,na.rm=TRUE)
{
  if(identical(max(x,na.rm=na.rm),min(x,na.rm=na.rm))) NA
  a <- 1/(max(x)-min(x))
  b <- -min(x)/(max(x)-min(x))
  a*x+b
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Plot text
##
## Plot text
## x character, text to plot
## cex
## forecolor color of foreground
## bg color of background
## bordercolor color of border
## axes as in \code{graphics:::plot}
## ... additional arguments for \code{graphics:::text}
.plot.text <- function(x,xlim=c(0,1),ylim=c(0,1),cex=1,forecolor=par("fg"),bg=par("bg"),bordercolor=NA,axes=FALSE,...){
  if (.invalid(x)){
    x <- NULL
  }
  if (is.null(x)){
    x <- ""
  } else if (is.na(x)){
    x <- 'NA'
  }

  plot(xlim,ylim,type="n",ylab="",xlab="",xaxt="n",yaxt="n",axes=axes)
  rect(xleft=0, ybottom=0, xright=1, ytop=1, col=bg, border=bordercolor)
  text(0.5,0.5,x,cex=cex,...)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## This was originally heatmap.3.
heatmap.cnv <-
  function(x,

           ## whether a dissimilarity matrix
           diss=inherits(x,"dist"),

           ## dendrogram control
           Rowv=TRUE,
           Colv=TRUE,
           dendrogram=c("both","row","column","none"),

           ## dist object
           dist.row,
           dist.col,
           dist.FUN=gdist,
           dist.FUN.MoreArgs=list(method="euclidean"),

           ## hclust object
           hclust.row,
           hclust.col,
           hclust.FUN=hclust,
           hclust.FUN.MoreArgs=list(method="ward.D"),

           ## data scaling
           scale=c("none","row","column"),
           na.rm=TRUE,

           ## clustering control
           cluster.by.row=FALSE,
           cluster.by.col=FALSE,
           kr=NA,
           kc=NA,
           row.clusters=NA,
           col.clusters=NA,

           ## image plot
           revR=FALSE,
           revC=FALSE,
           add.expr,

           ## mapping data to colors
           breaks,
           ## centering colors to a value
           x.center,
           ## colors
           color.FUN=gplots::bluered,
           ##
           ## block sepration
           sepList=list(NULL,NULL),
           sep.color=c("gray45","gray45"),
           sep.lty=1,
           sep.lwd=2,

           ## cell labeling
           cellnote,
           cex.note=1.0,
           notecol="cyan",
           na.color=par("bg"),

           ## level trace
           trace=c("none","column","row","both"),
           tracecol="cyan",
           hline,
           vline,
           linecol=tracecol,

           ## Row/Column Labeling
           labRow=TRUE, ## shown by default
           labCol=TRUE, ## shown by default
           srtRow=NULL,
           srtCol=NULL,
           sideRow=4,
           sideCol=1,
           ##
           margin.for.labRow,
           margin.for.labCol,
           ColIndividualColors,
           RowIndividualColors,
           annotations_legend,
           grouping_key_coln,
           cexRow,
           cexCol,
           labRow.by.group=FALSE,
           labCol.by.group=FALSE,


           ## plot color key + density info
           key=TRUE,
           key.title="Color Key",
           key.xlab="Value",
           key.ylab="Count",

           keysize=1.5,
           mapsize=9,
           mapratio=4/3,
           sidesize=3,
           cex.key.main=0.75,
           cex.key.xlab=0.75,
           cex.key.ylab=0.75,
           density.info=c("histogram","density","none"),
           denscol=tracecol,
           densadj=0.25,

           ## plot titles/labels
           main="Heatmap",
           sub="",
           xlab="",
           ylab="",
           cex.main=2,
           cex.sub=1.5,
           font.main=2,
           font.sub=3,
           adj.main=0.5,
           mgp.main=c(1.5,0.5,0),
           mar.main=3,
           mar.sub=3,
           ## plot ##
           if.plot=TRUE,

           ## plot of partition (left/top of heatmap)
           plot.row.partition=FALSE,
           plot.col.partition=FALSE,
           cex.partition=1.25,
           color.partition.box="gray45",
           color.partition.border="#FFFFFF",

           ## plot of summary (right/bottom of heatmap)
           plot.row.individuals=FALSE,
           plot.col.individuals=FALSE,
           plot.row.clusters=FALSE,
           plot.col.clusters=FALSE,
           plot.row.clustering=FALSE,
           plot.col.clustering=FALSE,

           ##
           plot.row.individuals.list=FALSE,
           plot.col.individuals.list=FALSE,
           plot.row.clusters.list=FALSE,
           plot.col.clusters.list=FALSE,
           plot.row.clustering.list=FALSE,
           plot.col.clustering.list=FALSE,

           ## for plot of clusters - row
           row.data=FALSE,
           ## for plot of clusters - col
           col.data=FALSE,

           ##
           if.plot.info=FALSE,
           text.box,
           cex.text=1.0,

           ## Force in layout info
           force_lmat=NULL,
           force_lwid=NULL,
           force_lhei=NULL,
           force_add=FALSE,

           ## extras
           ...
           )
{
  ## check input - take1 ##
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  x.ori <- x

  if(!inherits(x,"dist") & !is.matrix(x)){
    stop("`x' should either be a matrix, a data.frame or a `dist' object.")
  }

  if (! sideRow %in% c(2,4)){
    stop('sideRow must be either 2 or 4.')
  }

  if (! sideCol %in% c(1,3)){
    stop('sideCol must be either 1 or 3.')
  }

  ## store input
  Rowv.ori <- Rowv
  Colv.ori <- Colv

  ## check
  dendrogram <- match.arg(dendrogram)
  if ( (dendrogram %in% c("both","row")) & !inherits(Rowv,"dendrogram") ){
    warning("Discrepancy: row dendrogram is asked;  Rowv is set to `TRUE'.")
    Rowv <- TRUE
  }

  if ( (dendrogram %in% c("both","col")) & !inherits(Colv,"dendrogram") ){
    warning("Discrepancy: col dendrogram is asked;  Colv is set to `TRUE'.")
    Colv <- TRUE
  }

  if (identical(Rowv, FALSE) | missing(Rowv)){
    if(!identical(cluster.by.row,FALSE)){
      warning("Discrepancy: No row dendrogram is asked; cluster.by.row is set to `FALSE'.")
      cluster.by.row <- FALSE
    }
  } else {
    if(!identical(cluster.by.row,TRUE)){
      warning("Discrepancy: row dendrogram is asked; cluster.by.row is set to `TRUE'.")
      cluster.by.row <- TRUE
    }
  }

  if (identical(Colv, FALSE) | .invalid(Colv)){
    if(!identical(cluster.by.col,FALSE)){
      warning("Discrepancy: No col dendrogram is asked; cluster.by.col is set to `FALSE'.")
      cluster.by.col <- FALSE
    }
  } else {
    if(!identical(cluster.by.col,TRUE)){
      warning("Discrepancy: col dendrogram is asked; cluster.by.col is set to `TRUE'.")
      cluster.by.col <- TRUE
    }
  }

  if (!.invalid(kr)){
    if (is.numeric(kr)){
      if(!plot.row.partition){
        warning("Discrepancy: kr is set, therefore plot.row.partition is set to `TRUE'.")
        plot.row.partition <- TRUE
      }
    }
  }

  if (!.invalid(kc)){
    if (is.numeric(kc)){
      if(!plot.col.partition){
        warning("Discrepancy: kc is set, therefore plot.col.partition is set to `TRUE'.")
        plot.col.partition <- TRUE
      }
    }
  }

  ## generate dist.obj - row/col ##
  if (inherits(x,"dist")){
    dist.row <- dist.col <- x ## dist.obj
    x <- as.matrix(x)
    mat.row <- mat.col <- x ## mat.obj
    symm <- TRUE
  } else if (is.matrix(x)){
    symm <- isSymmetric(x)
    if (diss){
      if (!symm){
        stop("Dissimilary matrix should be symmetric. Please set `diss' to FALSE if `x' is not dissimilary matrix.")
      } else {
        flush.console()
      }
      mat.row <- mat.col <- x
      dist.row <- dist.col <- as.dist(x)
    } else{
      if (cluster.by.row) {
        if (.invalid(dist.row)){
          dist.row <- .call.FUN(dist.FUN,x,MoreArgs=dist.FUN.MoreArgs)
        }
        mat.row <- as.matrix(dist.row)
      } else {
        dist.row <- NULL
        mat.row <- NULL
      }
      if (cluster.by.col) {
        if (.invalid(dist.col)){
          dist.col <- .call.FUN(dist.FUN,t(x),MoreArgs=dist.FUN.MoreArgs)
        }
        mat.col <- as.matrix(dist.col)
      } else {
        dist.col <- NULL
        mat.col <- NULL
      }
    }
  }
  ## check input - take2: di ##
  di <- dim(x)

  if(length(di)!=2 || !is.numeric(x)){
    stop("`x' should only contain `numeric' values and can be converted to a 2-D matrix.")
  }

  ## parse param ##
  scale <- if(symm && .invalid(scale)) "none" else match.arg(scale) ## no scale on symmetric matrix
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  dist.FUN <- match.fun(dist.FUN)
  hclust.FUN <- match.fun(hclust.FUN)
  color.FUN <- match.fun(color.FUN)

  ## NG if both breaks and scale are specified ##
  if(!.invalid(breaks) & (scale!="none")){
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  }

  ## nr and nc ##
  nr <- di[1]
  nc <- di[2]

  ## check input - take3: nr,nc ##
  if(nr <=1 || nc <=1)
    stop("`x' must have at least 2 rows and 2 columns")

  ## font size of row/col labels ##
  cexRow0 <- 0.2+1/log10(nr)
  cexCol0 <- 0.2+1/log10(nc)

  if (.invalid(cexRow)) {
    cexRow <- cexRow0
  } else {
    # message('heatmap.3 | From GMD 0.3.3, please use relative values for cexRow.')
    cexRow <- cexRow0*cexRow
  }
  if (.invalid(cexCol)) {
    cexCol <- cexCol0
  } else {
    # message('heatmap.3 | From GMD 0.3.3, please use relative values for cexCol.')
    cexCol <- cexCol0*cexCol
  }

  ## cellnote ##
  ## ##if(.invalid(cellnote)) cellnote <- matrix("",ncol=ncol(x),nrow=nrow(x))

  ## ------------------------------------------------------------------------
  ## parse dendrogram ##
  ## ------------------------------------------------------------------------

  if (missing(Rowv)) Rowv <- FALSE
  if (.invalid(Colv)) Colv <- if(symm) Rowv else FALSE
  if (Colv=="Rowv") {
    if ((!isTRUE(Rowv) | !symm) ){
      Colv <- FALSE
      warning("`Colv' is specified to use \"Rowv\", but either `Rowv' is invalid or `x' is not symmetric; Colv is suppressed.")
    } else{
      Colv <- Rowv
    }
  }
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  flush.console()

  if ( (!inherits(Rowv,"dendrogram") & !identical(Rowv,FALSE)) | (cluster.by.row & .invalid(row.clusters))){
    if (.invalid(hclust.row)){
      hclust.row <- .call.FUN(hclust.FUN,dist.row,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.row$order) != nr){
        stop("`hclust.row' should have equal size as the rows.")
      }
    }
  } else{
    hclust.row <- NULL
  }

  if(symm){
    hclust.col <- hclust.row
  }

  if ( !inherits(Colv,"dendrogram") & !identical(Colv,FALSE) | (cluster.by.col & .invalid(col.clusters))){
    if (.invalid(hclust.col)){
      hclust.col <- .call.FUN(hclust.FUN,dist.col,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.col$order) != nc){
        stop("`hclust.col' should have equal size as the cols.")
      }
    }
  } else {
    hclust.col <- NULL
  }
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  ddr <- ddc <- NULL
  ## get the dendrograms and reordering row/column indices - row ##
  if(inherits(Rowv,"dendrogram")){
    if (attr(Rowv,"members") != nr){
      stop("`Rowv' should contain equal size of members as the rows.")
    }
    ddr <- Rowv ## use Rowv 'as-is',when it is dendrogram
    # rowInd <- order.dendrogram(ddr)
    rowInd <- 1:nr  ### data has already been sorted in the proper order, and coloring vectors too
  } else if (is.integer(Rowv)){ ## compute dendrogram and do reordering based on given vector
    ddr <- as.dendrogram(hclust.row)
    ddr <-  reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else if (isTRUE(Rowv)){ ## if TRUE,compute dendrogram and do reordering based on rowMeans

    Rowv <- rowMeans(x,na.rm=TRUE)
    ddr <- as.dendrogram(hclust.row)
    ddr <- reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr !=length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else{
    rowInd <- nr:1 ## from bottom.
  }
  ## get the dendrograms and reordering row/column indices - col ##
  if(inherits(Colv,"dendrogram")){
    if (attr(Colv,"members") != nc){
      stop("`Colv' should contain equal size of members as the cols.")
    }
    ddc <- Colv ## use Colv 'as-is',when it is dendrogram
    colInd <- order.dendrogram(ddc)
  } else if(identical(Colv,"Rowv")) {
    if(exists("ddr")){
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else{
      colInd <- rowInd
    }
  } else if(is.integer(Colv)){## compute dendrogram and do reordering based on given vector
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("`colInd' is of wrong length.")
  } else if (isTRUE(Colv)){## if TRUE,compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x,na.rm=TRUE)
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc !=length(colInd))
      stop("`colInd' is of wrong length.")
  } else{
    colInd <- 1:nc ## from left
  }
  ## ------------------------------------------------------------------------
  ## check consistency
  ## ------------------------------------------------------------------------

  ## Xmisc::logme(dendrogram)
  ## Xmisc::logme(Colv)
  ## Xmisc::logme(Rowv)

  ## dendrogram - check consistency: Rowv ##
  if ( is.null(ddr) & (dendrogram %in% c("both","row"))){
    warning("Discrepancy: Rowv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting row dendogram.")
    if (is.logical(Colv) & (Colv.ori) & dendrogram=="both")
      dendrogram <- "column"
    else
      dendrogram <- "none"
  }

  ## dendrogram - check consistency: Colv ##
  if ( is.null(ddc) & (dendrogram %in% c("both","column"))){
    warning("Discrepancy: Colv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting column dendogram.")
    if (is.logical(Rowv) & (identical(Rowv.ori,TRUE) | is.numeric(Rowv.ori) | inherits(Rowv.ori,"dendrogram")) & dendrogram=="both")
      dendrogram <- "row"
    else
      dendrogram <- "none"
  }

  ## check consistency
  if (is.null(ddr)){
    if(isTRUE(cluster.by.row) | isTRUE(plot.row.partition) | isTRUE(plot.row.clusters) | isTRUE(plot.row.clustering) ){
      warning("Using invalid `Rowv' while allowing",
              "`cluster.by.row' or `plot.row.partition' or `plot.row.clusters' or `plot.row.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }
  if (is.null(ddc)){
    if(isTRUE(cluster.by.col) | isTRUE(plot.col.partition) | isTRUE(plot.col.clusters) | isTRUE(plot.col.clustering) ){
      warning("Using invalid `Colv' while allowing",
              "`cluster.by.col' or `plot.col.partition' or `plot.col.clusters' or `plot.col.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }

  if (is.null(ddr)) cluster.by.row <- plot.row.partition <- plot.row.clusters <- plot.row.clustering <- FALSE
  if (is.null(ddc)) cluster.by.col <- plot.col.partition <- plot.col.clusters <- plot.col.clustering <- FALSE
  ## ------------------------------------------------------------------------
  ## Reordering
  ## ------------------------------------------------------------------------
  flush.console()
  ## reorder x and cellnote ##
  x <- x[rowInd,colInd]

  if (!.invalid(cellnote)) cellnote <- cellnote[rowInd,colInd]

  ## reorder labels - row ##
  if(identical(labRow,TRUE)){ ## Note: x is already reorderred
    labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  } else if(identical(labRow,FALSE) | .invalid(labRow)){
    labRow <- rep("",nrow(x))
  } else if(is.character(labRow)){
    labRow <- labRow[rowInd]
  }

  ## reorder cellnote/labels - col ##
  if (identical(labCol,TRUE)){
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  } else if(identical(labCol,FALSE) | .invalid(labCol)){
    labCol <- rep("",ncol(x))
  } else if(is.character(labCol)){
    labCol <- labCol[colInd]
  }

  ## ------------------------------------------------------------------------
  ## scale
  ## center to 0 and scale to 1 in row or col but not both! ##
  ## ------------------------------------------------------------------------
  flush.console()
  x <- .scale.data(x,scale,na.rm)
  ## ------------------------------------------------------------------------
  ## labels for observations/clusters/
  ## ------------------------------------------------------------------------
  ## margin for labels

  margin.for.labRow0 <- max(nchar(labRow, keepNA = FALSE))*0.75+0.2
  margin.for.labCol0 <- max(nchar(labCol, keepNA = FALSE))*0.75+0.2  ## Change compared to GMD to work with R 3.3+ that changed the default value of keepNA to "NA"

  if (.invalid(margin.for.labRow)){
    margin.for.labRow <- margin.for.labRow0
  }
  #} else {
  #  message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labRow.')
  #  margin.for.labRow <- margin.for.labRow0*margin.for.labRow
  #}
  if (margin.for.labRow < 2) {
      margin.for.labRow <- 2
  }

  # if (.invalid(margin.for.labCol)){
  margin.for.labCol <- margin.for.labCol0
  # }
  #  else {
  #   # message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labCol.')
  #   margin.for.labCol <- margin.for.labCol0*margin.for.labCol
  # }

  if (margin.for.labCol < 2) {
    margin.for.labCol <- 2
  }

  ## group unique labels - row ## ##??check
  if (!.invalid(labRow.by.group) & !identical(labRow.by.group,FALSE)){
    group.value <- unique(labRow)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labRow)
    labRow <- rep("",length(labRow))
    labRow[group.index] <- group.value
  }

  ## group unique labels - col ## ##??check
  if (!.invalid(labCol.by.group) & !identical(labCol.by.group,FALSE)){
    group.value <- unique(labCol)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labCol)
    labCol <- rep("",length(labCol))
    labCol[group.index] <- group.value
  }
  ## ------------------------------------------------------------------------
  ## color breaks
  ## ------------------------------------------------------------------------
  flush.console()

  ## set breaks for binning x into colors ##
  if(.invalid(breaks)){
    breaks <- 16
  } else {
      logging::logdebug(paste("inferCNV::heatmap.cnv, breaks parameter set to: [", paste(breaks, collapse=","), "]", sep=""))
  }

  ## get x.range according to the value of x.center ##
  if (!.invalid(x.center)){ ## enhanced
    if (is.numeric(x.center)){
      x.range.old <- range(x,na.rm=TRUE)
      if (length(breaks) > 1) {
          # important, use specified breakpoint info here if set by user
          x.range.old = range(breaks)
      }
      dist.to.x.center <- max(abs(x.range.old-x.center))
      x.range <- c(x.center-dist.to.x.center,x.center+dist.to.x.center)
      if (length(breaks) > 1) {
          # re-set the breaks according to the new x.range
          breaks=seq(x.range[1], x.range[2], length=16)
          logging::logdebug(paste("inferCNV::heatmap.cnv, resetting breaks to adjusted x.range: [",
                                  paste(breaks, collapse=","), "]", sep=""))
      }

    } else {
      stop("`x.center' should be numeric.")
    }
  } else{
      x.range <- range(x,na.rm=TRUE)
      if (length(breaks) > 1) {
          # important, use specified breakpoint info here if set by user
          x.range = range(breaks)
      }

  }
  logging::logdebug( paste("inferCNV::heatmap.cnv x range set to: ",
                          paste(x.range, collapse=",")), sep="" )

  ## set breaks for centering colors to the value of x.center ##
  if(length(breaks)==1){
    breaks <-
      seq(min(min(x,na.rm=TRUE),x.range[1]),
          max(max(x,na.rm=TRUE),x.range[2]),
          length.out=breaks)
  }

  ## count of breaks and colors ##
  nbr <- length(breaks)
  ncolor <- length(breaks)-1

  ## generate colors ##
  colors <- color.FUN(ncolor)
  ## set up breaks and force values outside the range into the endmost bins ##
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x<min.breaks,min.breaks,x)
  x[] <- ifelse(x>max.breaks,max.breaks,x)
  ## ------------------------------------------------------------------------
  ## check if it is sufficient to draw side plots ##
  ## ------------------------------------------------------------------------
  if (cluster.by.row){
    if (!.invalid(row.clusters)) {## suppress kr
      if(!is.numeric(row.clusters) | length(row.clusters)!=nr | !(.is.grouped(row.clusters))){
        warning("`row.clusters' is not a grouped numeric vector of length nrow(x); cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      } else{
        row.clusters <- row.clusters[rowInd]
        kr <- length(unique(row.clusters))
      }
    } else {
      if (.invalid(kr)) kr <- 2
      if (is.numeric(kr) & length(kr)==1){
        row.clusters <- cutree(hclust.row,k=kr)
        row.clusters <- row.clusters[rowInd]
      } else {
        warning("`kr' should be numeric of length one; cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      }
    }
  }

  if (cluster.by.col){
    if (!.invalid(col.clusters)) {## suppress kc
      if(!is.numeric(col.clusters) | length(col.clusters)!=nc | !(.is.grouped(col.clusters))){
        warning("`col.clusters' is not a grouped numeric vector of length ncol(x); cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      } else{
        col.clusters <- col.clusters[colInd]
        kc <- length(unique(col.clusters))
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      }
    } else {
      if (.invalid(kc)) kc <- 2
      if (is.numeric(kc) & length(kc)==1){
        col.clusters <- cutree(hclust.col,k=kc)
        col.clusters <- col.clusters[colInd]
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      } else {
        warning("`kc' should be numeric of length one; cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      }
    }
  }
  ## ------------------------------------------------------------------------
  ## Plotting
  ## ------------------------------------------------------------------------
  if (if.plot){

    ir <- length(plot.row.individuals.list)
    ic <- length(plot.col.individuals.list)
    cr <- length(plot.row.clustering.list)
    cc <- length(plot.col.clustering.list)

    flush.console()
    if(mapratio<=1){
      sr <- 12
      sc <- sr*mapratio
    } else {
      sc <- 12
      sr <- sc/mapratio
    }

    ## calculate the plot layout ##
    ## 1) for heatmap
    lmat <- matrix(1,nrow=sr,ncol=sc)
    lwid <- c(rep(mapsize/sc,sc))
    lhei <- c(rep(mapsize/mapratio/sr,sr))

    ## 2) row.clusters
    if (plot.row.partition | plot.row.clusters){
      lmat <- cbind(max(lmat,na.rm=TRUE)+1,lmat)
      lwid <- c(0.3,lwid)
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid)

    }
    ## 3) col.clusters
    if (plot.col.partition | plot.col.clusters){
      lmat <- rbind(c(NA,rep(max(lmat,na.rm=TRUE)+1,sc)),lmat)
      lhei <- c(0.3/mapratio,lhei)
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }

    if(!.invalid(RowIndividualColors)) { ## 4) add middle column to layout for vertical sidebar ##??check
      # if(!is.character(RowIndividualColors) || length(RowIndividualColors) !=nr)
      if(!is.character(RowIndividualColors) || dim(RowIndividualColors)[1] !=nr)
        stop("'RowIndividualColors' must be a character vector of length nrow(x)")
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)
      lwid <- c(0.2,lwid)
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)  # not really useful as this if is only true in plot_cnv_observations which forces layout_lmat and layout_lhei
      lwid <- c(0.2,lwid) # layout_lwid is however not forced, so need to actually update it here
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid)
    }

    if(!.invalid(ColIndividualColors)) { ## 5) add middle row to layout for horizontal sidebar ##??check
      if(!is.character(ColIndividualColors) || length(ColIndividualColors) !=nc){
        stop("'ColIndividualColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(c(rep(NA,ncol(lmat)-sc),rep(max(lmat,na.rm=TRUE)+1,sc)),lmat)
      lhei <- c(0.2/mapratio,lhei)
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }
    ## 6) for row-dend
    lmat <- cbind(c(rep(NA,nrow(lmat)-sr),
                    rep(max(lmat,na.rm=TRUE)+1,sr)),
                  lmat
                  )
    lwid <- c(keysize,lwid)

    ## 7) for col-dend, 8) for kay
    lmat <- rbind(c(
                    max(lmat,na.rm=TRUE)+2,
                    rep(NA,ncol(lmat)-sc-1),
                    rep(max(lmat,na.rm=TRUE)+1,sc)
                    ),
                  lmat
                  )
    lhei <- c(keysize/mapratio,lhei)
    ## text.box##
    ## numbered 999 ##
    ## 9) for RowPlot (from bottom)
    if(.invalid(text.box)){
      text.box <- "made by\nFunction: heatmap.3\nPackage: GMD\nin R"
    }
    if(plot.row.individuals) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),nrow(lmat)-sr),# text
                      rep((ir:1)+max(lmat,na.rm=TRUE)+(1),each=floor(sr/ir)),rep(NA,sr%%ir)
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    ## 10) for ColPlot from right
    if(plot.col.individuals) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),ncol(lmat)-sc-1),# text
                      rep((1:ic)+max(lmat,na.rm=TRUE)+(1),each=floor(sc/ic)),rep(NA,sc%%ic),
                      ##NA # change to numeric if text.box
                      999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }
    ## 11) for RowPlot (from bottom)
    if(plot.row.clusters) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-1), # text
                      rep((kr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sr/kr)),rep(NA,sr%%kr),
                      ##NA
                      999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    ## 12) for ColPlot from right
    if(plot.col.clusters) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-2),# text
                      rep((1:kc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/kc)),rep(NA,sc%%kc),
                      ##NA,NA # change to numeric if text.box
                      999,999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    ## 13) for RowPlot (from bottom)
    if(plot.row.clustering) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-2), # text
                      rep(c((cr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1)),each=floor(sr/cr)),rep(NA,sr%%cr),
                      ##NA,NA
                      999,999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }


    ## 14) for ColPlot from right
    if(plot.col.clustering) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-3),# text
                      rep((1:cc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/cc)),rep(NA,sc%%cc),
                      ##NA,NA,NA # change to numeric if text.box
                      999,999,999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    lmat[is.na(lmat)] <- 0
    if (any(lmat==999)) flag.text <- TRUE else flag.text <- FALSE
    lmat[lmat==999] <- max(lmat[lmat!=999])+1

    ## Graphics `output' ##
    ## layout
    if(!is.null(force_lmat)){
        lmat <- force_lmat
    }
    if(!is.null(force_lwid)){
        lwid <- force_lwid
    }
    if(!is.null(force_lhei)){
        lhei <- force_lhei
    }
    if(!force_add){
        layout(lmat,widths=lwid,heights=lhei,respect=FALSE)
    }

    ## reverse columns
    if(revC){ ## x columns reversed
      iy <- nr:1
      ddc <- rev(ddc)
      x <- x[iy,]
      if (!.invalid(cellnote)) cellnote <- cellnote[iy,]
    } else {
      iy <- 1:nr
    }

    ## reverse rows
    if(revR){ ## x columns reversed
      ix <- nc:1
      ddr <- rev(ddr)
      x <- x[,ix]
      if (!.invalid(cellnote)) cellnote <- cellnote[,ix]
    } else {
      ix <- 1:nc
    }

    ## 1) draw the main carpet/heatmap
    margins <- c(margin.for.labCol,0,0,margin.for.labRow)
    mgp <- c(3,1,0)
    par(mar=margins,mgp=mgp);outer=FALSE

    x.save <- x
    if(!symm || scale !="none"){ ##??
      x <- t(x)
      if (!.invalid(cellnote)) cellnote <- t(cellnote)
    }
    image(1:nc,1:nr,
          x,
          xlim=0.5+c(0,nc),ylim=0.5+c(0,nr),
          axes=FALSE,xlab="",ylab="",col=colors,breaks=breaks,
          ...)


    ## plot/color NAs
    if(!.invalid(na.color) & any(is.na(x))){
      mmat <- ifelse(is.na(x),1,NA)
      image(1:nc,1:nr,mmat,axes=FALSE,xlab="",ylab="",
            col=na.color,add=TRUE)
    }

    ##
    ## labCol (?)
    if ((dendrogram %in% c("both","col")) & sideCol==3) {
      warning("Discrepancy: col dendrogram is asked; srtCol is set to 1.")
      sideCol <- 1
    }
    if (!length(srtCol)) {
      axis(sideCol,1:nc,labels=labCol,las=2,line=-0.5,tick=0,cex.axis=cexCol,outer=outer)
    } else {
      if (sideCol==1){
        if (sideCol==1) .sideCol <- par("usr")[3]-0.5*srtCol/90 else .sideCol <- par("usr")[4]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }

    if(!.invalid(xlab)) mtext(xlab,side=1,line=margins[1]-1.25)

    ## labRow (?)
    if ((dendrogram %in% c("both","row")) & sideRow==2) {
      warning("Discrepancy: row dendrogram is asked; sideRow is set to 4.")
      sideRow <- 4
    }
    if (!length(srtRow)) {
      axis(sideRow,iy,labels=labRow,las=2,line=-0.5,tick=0,cex.axis=cexRow,outer=outer)
    } else {
      if (sideRow==4){
        if (sideRow==4) .sideRow <- par("usr")[2]+0.5*srtRow/90 else .sideRow <- par("usr")[1]-0.5*srtRow/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }

    if(!.invalid(ylab)) mtext(ylab,side=4,line=margins[4]-1.25)

    if (!.invalid(add.expr))
      eval(substitute(add.expr))
    ## Enhanced: add 'sep.color' colored spaces to visually separate sections
    if (plot.row.partition | plot.row.clusters){ ##??
      plot.row.partitionList <- get.sep(clusters=row.clusters,type="row")
    } else {
      plot.row.partitionList <- NULL
    }
    if (plot.col.partition | plot.col.clusters){ ##??
      plot.col.partitionList <- get.sep(clusters=col.clusters,type="column")
    } else {
      plot.col.partitionList <- NULL
    }

    row.sepList <- sepList[[1]]
    if (!.invalid(row.sepList)){
      for (i in 1:length(row.sepList)){
        i.sep <- row.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[1]
             )
      }
    }
    col.sepList <- sepList[[2]]
    if (!.invalid(col.sepList)){
      for (i in 1:length(col.sepList)){
        i.sep <- col.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[2]
             )
      }
    }
    ## show traces
    min.scale <- min(breaks)
    max.scale <- max(breaks)

    x.scaled  <- .scale.x(t(x),min.scale,max.scale)

    if(.invalid(hline)) hline=median(breaks)
    if(.invalid(vline)) vline=median(breaks)

    if(trace %in% c("both","column")){
      for( i in colInd ){
        if(!.invalid(vline)){
          vline.vals <- .scale.x(vline,min.scale,max.scale)
          abline(v=i-0.5+vline.vals,col=linecol,lty=2)
        }
        xv <- rep(i,nrow(x.scaled))+x.scaled[,i]-0.5
        xv <- c(xv[1],xv)
        yv <- 1:length(xv)-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }

    if(trace %in% c("both","row")){
      for( i in rowInd ){
        if(!.invalid(hline)){
          hline.vals <- .scale.x(hline,min.scale,max.scale)
          abline(h=i+hline,col=linecol,lty=2)
        }
        yv <- rep(i,ncol(x.scaled))+x.scaled[i,]-0.5
        yv <- rev(c(yv[1],yv))
        xv <- length(yv):1-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }
    ## cellnote
    if(!.invalid(cellnote)){
      text(x=c(row(cellnote)),
           y=c(col(cellnote)),
           labels=c(cellnote),
           col=notecol,
           cex=cex.note)
    }

    ## 2) plot.row.partition
    if(plot.row.partition |plot.row.clusters) { ##row.clusters
      par(mar=c(margins[1],0.5,0,0.1))

      row.clusters.unique <- unique(row.clusters)
      row.clusters.unique <- row.clusters.unique[!is.na(row.clusters.unique)]

      image(rbind(1:nr),
            xlim=0.5+c(0,1),ylim=0.5+c(0,nr),
            col=par("bg"),
            axes=FALSE,add=force_add)
      if (!.invalid(plot.row.partitionList)){
        for (i in 1:length(plot.row.partitionList)){
          i.sep <- plot.row.partitionList[[i]]
          rect(
               xleft=0+0.5,
               ybottom=i.sep[2]+0.5,
               xright=1+0.5,
               ytop=i.sep[4]+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- row.clusters.unique[i]
          s <- g
          text(x=1,y=(i.sep[2]+0.5+i.sep[4]+0.5)/2,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=90
               )
        }
      }

    }
    ## 3) plot.col.partition
    if(plot.col.partition | plot.col.clusters) {
      par(mar=c(0.1,0,0,margins[4]))
      col.clusters.unique <- unique(col.clusters)
      col.clusters.unique <- col.clusters.unique[!is.na(col.clusters.unique)]

      image(cbind(1:nc),
            xlim=0.5+c(0,nc),ylim=0.5+c(0,1),
            col=par("bg"),
            axes=FALSE,add=force_add)

      if (!.invalid(plot.col.partitionList)){
        for (i in 1:length(plot.col.partitionList)){
          i.sep <- plot.col.partitionList[[i]]
          rect(
               xleft=i.sep[1]+0.5,
               ybottom=0+0.5,
               xright=i.sep[3]+0.5,
               ytop=1+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- col.clusters.unique[i]
          s <- g
          text(x=(i.sep[1]+0.5+i.sep[3]+0.5)/2,y=1,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=0
               )
        }
      }

    }
    ## 4) draw the side color bars - for row
    if(!.invalid(RowIndividualColors)) {
      par(mar=c(margins[1],0,0,0.5))
      # image(rbind(1:nr),col=RowIndividualColors[rowInd, 1],axes=FALSE,add=force_add)
      image(rbind(1:nr),col=RowIndividualColors[rowInd, 1],axes=FALSE,add=FALSE)

        if (dim(RowIndividualColors)[2] > 1) {
        par(mar=c(margins[1],0,0,0.5))
            image(rbind(1:nr),col=RowIndividualColors[rowInd, 2], axes=FALSE, add=force_add)
        }
    }

    ## 5) draw the side color bars - for col
    if(!.invalid(ColIndividualColors)) {
      par(mar=c(0.5,0,0,margins[4]))
      image(cbind(1:nc),col=ColIndividualColors[colInd],axes=FALSE,add=force_add)
    }

    ## 6) row-dend
    par(mar=c(margins[1],0,0,0))
    if(dendrogram %in% c("both","row")){
      plot(ddr,horiz=TRUE,axes=FALSE,yaxs="i",leaflab="none")
    }else{
      .plot.text(ylim=range(iy))
      if (sideRow==2){
        .sideRow <- par("usr")[2]-0.5*srtCol/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }


    par(mar=c(0,0,0,0))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend(x=0.6,y=1.2 ,legend=annotations_legend[,1],
        cex=1.2,
        fill=annotations_legend[,2],
        ncol=grouping_key_coln)


    ## 7) col-dend and title
    mar3 <- (if(!is.null(main)) mar.main else 0) +
      (if(!is.null(sub)) mar.sub else 0)
    par(mar=c(0,0,mar3,margins[4]))

    if(dendrogram %in% c("both","column")){
      plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    } else{
      if(key) {
        .plot.text(xlim=range(1:nc))
      }
      if (sideCol==3){
        .sideCol <- par("usr")[3]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }


    ## title
    if (is.null(sub)) main.line <- 1 else main.line <- 3
    if(!is.null(main)) title(main,cex.main=cex.main,adj=adj.main,mgp=mgp.main,font.main=font.main,line=main.line)
    if(!is.null(sub)) title(sub,cex.main=cex.sub,adj=adj.main,mgp=mgp.main,font.main=font.sub,line=0)
    ##if(!is.null(main)) title(main,cex.main=1.5*op[["cex.main"]])
    ## 8) plot the color-key
    if(key){
      cex.key <- 0.75
      op.ori <- par()

      par(mar=c(2,1.5,0.75,1)*keysize,cex=cex.key,mgp=c(0.75,0,0),tcl=-0.05)
      z <- seq(x.range[1],x.range[2],length=length(colors))
      logging::logdebug(paste("::inferCNV::heatmap.cnv colorkey z range: ", paste(z, collapse=","), sep=""))
      logging::logdebug(paste("::inferCNV::heatmap.cnv colorkey breaks range: ", paste(breaks, collapse=","), sep=""))
      logging::logdebug(paste("::inferCNV::heatmap.cnv colorkey colors range: ", paste(colors, collapse=","), sep=""))

      image(z=matrix(z,ncol=1),
            col=colors,
            breaks=breaks,
            xaxt="n",
            yaxt="n",
            xlab=key.xlab,
            ylab="",
            main="",add=force_add
            )

      par(usr=c(0,1,0,1))
      lv <- pretty(breaks)
      xv <- .scale.x(as.numeric(lv),x.range[1],x.range[2])
      axis(1,at=xv,labels=lv,cex.axis=cex.key*1)
      if(density.info=="density"){
        ## Experimental : also plot density of data
        dens <- density(x,adjust=densadj,na.rm=TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- .scale.x(dens$x,x.range[1],x.range[2])
        lines(dens$x,dens$y / max(dens$y) * 0.95,col=denscol,lwd=1)
        axis(2,at=pretty(dens$y)/max(dens$y) * 0.95,pretty(dens$y),cex.axis=cex.key*1)
        ##title("Color Key and Density",cex.lab=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,"Density",line=0.75,cex=cex.key)
      } else if(density.info=="histogram"){
        h <- hist(x,plot=FALSE,breaks=breaks)
        hx <- .scale.x(breaks,x.range[1],x.range[2])
        hy <- c(h$counts,h$counts[length(h$counts)])
        lines(hx,hy/max(hy)*0.95,lwd=1,type="s",col=denscol)
        axis(2,at=pretty(hy)/max(hy)*0.95,pretty(hy),cex.axis=cex.key*1)
        ##title("Color Key and Histogram",cex.main=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,key.ylab,line=0.75,cex=cex.key)
      } else{
        title(key.title,cex.main=cex.key,font.main=1)
      }
    } else{
      if(!force_add){
      .plot.text()
      }
    }


    ## 9)
    if(plot.row.individuals) {
      .plot.text("Row\nIndividuals",cex=cex.text,bg="white")
      for (i in 1:ir) {
        ##.plot.text()
        tmp <- plot.row.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 10)
    if(plot.col.individuals) {
      .plot.text("Column\nIndividuals",cex=cex.text,bg="white",srt=90)
      for (i in 1:ic) {
        ##.plot.text()
        tmp <- plot.col.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 11) for RowPlot from bottom
    if (plot.row.clusters){
      .plot.text("Row\nClusters",cex=cex.text,bg="white")

      tmp <- plot.row.clusters.list[[1]]
      row.data <- row.data[rowInd]
      for (i in unique(row.clusters)){
        i.x <- row.data[row.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Row group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }
    ## 12) for ColPlot from left
    if (plot.col.clusters){
      .plot.text("Col\nClusters",cex=cex.text,bg="white",srt=90)

      tmp <- plot.col.clusters.list[[1]]
      col.data <- if(revC) col.data[rev(colInd)] else col.data[colInd]
      for (i in unique(col.clusters)){
        i.x <- col.data[col.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Col group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }

    ## 13)
    if(plot.row.clustering) {
      .plot.text("Row\nClustering",cex=cex.text,bg="white")

      for (i in 1:cr) {
        ##.plot.text()
        tmp <- plot.row.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    ## 14)
    if(plot.col.clustering) {
      .plot.text("Column\nClustering",cex=cex.text,bg="white",srt=90)

      for (i in 1:cc) {
        ##.plot.text()
        tmp <- plot.col.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 15) text
    if (!.invalid(text.box) & if.plot.info){
      .plot.text(text.box,cex=cex.text,bg="gray75")
    } else{
      if (flag.text){
        .plot.text()
      }
    }

  } else {
    x.save <- x
  }

  ret <-
    list(x.ori=x.ori,
         x=x.save,
         rowInd=rowInd,colInd=colInd,
         row.clusters=row.clusters,col.clusters=col.clusters,
         dist.row=dist.row,dist.col=dist.col,
         hclust.row=hclust.row,hclust.col=hclust.col,
         kr=kr,kc=kc
         )
  class(ret) <- c("hclustering",class(ret))
  invisible(ret)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## The official version from the package GMD is now archived
## https://cran.r-project.org/web/packages/GMD/index.html
gdist <-
  function(x,
           method="euclidean",
           MoreArgs=NULL,
           diag=FALSE,
           upper=FALSE
           )
{
  if(method %in% c("correlation","correlation.of.observations")){
    FUN <- function(x,...){
      as.dist(1-cor(t(x),y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  } else if(method %in% c("correlation.of.variables")){
    FUN <- function(x,...){
      as.dist(1-cor(x,y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  }

  COMMON_METHODS <-
    c("euclidean","maximum",
      "manhattan","canberra",
      "binary","minkowski"
      )

  if(method %in% COMMON_METHODS){
    d <- dist(x=x,method=method,diag=diag,upper=upper,p=MoreArgs$p)
  } else if (method %in% c("correlation","correlation.of.observations","correlation.of.variables")){
    ##d <- .call.FUN(FUN,x,MoreArgs)
    d <- FUN(x,method=MoreArgs$method,use=MoreArgs$use)
    attr(d,"method") <- method
  } else {
    FUN <- match.fun(method)
    MoreArgs[["diag"]] <- diag
    MoreArgs[["upper"]] <- upper
    d <- .call.FUN(FUN,x,MoreArgs)

    ## check attributes of the dist object ##
    if(is.null(attr(d,"method"))){
      attr(d,"method") <- method
    }
    if(is.null(attr(d,"call"))){
      attr(d,"call") <- match.call()
    }
    if(is.null(attr(d,"Size"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Size'.",attr(d,"method")))
    }
    if(is.null(attr(d,"Labels"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Labels'.",attr(d,"method")))
    }
  }
  attr(d,"Diag") <- diag
  attr(d,"Upper") <- upper
  class(d) <- "dist"
  return(d)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Get row or column lines of separation for \code{heatmap.3} according to clusters
## param clusters a numerical vector, indicating the cluster labels of observations.
## param type string, one of the following: \code{c("row","column","both")}
get.sep <-
  function(clusters,type=c("row","column","both"))
{
  ##   if(!is.numeric(clusters) | !(.is.grouped(clusters))){
  ## stop("`clusters' should be a grouped numeric vector.")
  ##   }
  tmp.whichna <- which(is.na(clusters))
  tmp.which <- which(!duplicated(clusters))

  tmp.sep <- data.frame(start=tmp.which,end=c(tmp.which[-1],length(clusters)+1)-1)
  tmp.sep2 <- tmp.sep[tmp.sep$start<=tmp.sep$end,]

  ## lines of separation
  sepList <- list()
  if (type=="row"){
    xmax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(0,tmp.sep2[i.s,1]-1,xmax,tmp.sep2[i.s,2])
    }
  } else if (type=="column"){
    ymax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,0,tmp.sep2[i.s,2],ymax)
    }
  } else if (type=="both"){
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,tmp.sep2[i.s,1]-1,tmp.sep2[i.s,2],tmp.sep2[i.s,2])
    }
  }
  sepList
}




#' Wrapper method that verifies settings and prepares the input for the other methods.
#' Use this function unless you know what you are doing.
#'
#' @title Run the infercnv analysis.
#'
#' @param x Data to plot (columns are cells and rows are genes). Can be either a file path or the actual data frame.
#' @param gene_order Position information for each gene used to smooth over chromosomes.
#' Can be either a file path or an actual data frame whose rows are genes and columns are the chromosome name,
#' the start position, and the end position.
#' @param annotations A data frame with an annotation for each cell. It can be used with the name_ref_groups option
#' to define which cells are to be used as references (and how they are grouped),
#' and is used for a color bar on the side of the observations heatmap.
#' @param use_color_safe To support the needs of those who see colors differently,
#' use this logical option to change the colors to a palette visibly distinct to all color blindness.
#' @param contig_label_size Used to increase or decrease the text labels for the X axis (contig names).
#' @param cutoff A number >= 0 is expected. A cut off for the average expression of genes to be used
#' for CNV inference (use the value before log2 transformation).
#' @param log_transform Matrix is assumed to be Log2(TPM+1) transformed.
#' If instead it is raw TPMs use this flag so that the data will be transformed.
#' @param delim Delimiter for reading expression matrix and writing matrices output.
#' @param noise_filter value for which +/- from zero will be set to zero to clear on heatmap
#' @param noise_quantiles quantile range for average bounds of residual reference values to
#' be set to zero to clear on heatmap (alternative to noise_filter)
#' @param max_centered_expression This value and -1 * this value are used as the maximum value
#' expression that can exist after centering data. If a value is outside of this range,
#' it is truncated to be within this range
#' @param num_obs_groups Number of groups in which to break the observations.
#' @param output_dir Directory in which to save plot and other output.
#' @param output_format Format in which to save plot. One of either "pdf" or "png"
#' @param num_ref_groups Define a number of groups to make automatically
#' by unsupervised clustering. This ignores annotations within references,
#' but does not mix them with observations.
#' @param name_ref_groups Names of groups from the "annotations" table whose cells
#' are to be used as reference groups. If num_ref_groups is not provided, groups will be kept.
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores num_obs_groups.
#' @param ref_subtract_method Method used to subtract the reference values from the observations.
#' Valid choices are: "by_mean", "by_quantiles".
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' @param clustering_contig When clustering observation samples, all genomic locations
#' are used unless this option is given. The expected value is one of the contigs (Chr)
#' in the genomic positions file (case senstive). All genomic positions will be plotted
#' but only the given contig will be used in clustering / group creation.
#' @param plot_steps Using this argument turns on plotting intemediate steps.
#' The plots will occur in the same directory as the output pdf.
#' Please note this option increases the time needed to run
#' @param bound_method_vis Method to automatically detect and bound outliers.
#' Used for visualizing. If both this argument and --vis_bound_threshold are given,
#' this will not be used. Valid choice is : "average_bound".
#' @param bound_threshold_vis Used as upper and lower bounds for values
#' in the visualization. If a value is outside this bound it will be replaced by the
#' closest bound. Should be given in the form of 1,1 (upper bound, lower bound).
#' @param window_length Window length for the smoothing.
#' @param contig_tail Contig tail to be removed.
#' @param fig_title Plot title.
#' @param obs_title Title for the observations matrix.
#' @param ref_title Title for the reference matrix.
#' @param load_workspace Try to load workspace saved from a previous run as far in the run as possible (checks for that important arguments match). Possible values are 0 "don't load any previous data", 1 "load preprocessed data", >2 "load processed data".
#' @param log_level Logging level.
#' @param ngchm Logical to decide whether to create a Next Generation Clustered Heat Map.
#' @param path_to_shaidyMapGen Path to the java application ShaidyMapGen.jar.
#' @param gene_symbol Specify the label type that is given to the gene needed to create linkouts, default is NULL,
#' @param min_cells_per_gene Minimum number of cells to be expressed for a gene in the reference set to be retained.
#' @param use_zscores Whether to calculate z scores from expression values for analysis or not.
#' @param make_zero_NA Replaces 0 values in the matrix by NA.
#'
#' @return
#' No return, void.
#' @export
infercnv <-
  function(x,  ## input_matrix, accepted both as a data.frame or a text file to read
        gene_order,  ## gene_order, accepted both as a data.frame or a text file to read
        annotations=NULL,
        use_color_safe=FALSE,
        contig_label_size=1,
        cutoff=0,
        log_transform=FALSE,
        delim="\t",
        noise_filter=NA,
        noise_quantiles=c(0.025,0.975),
        max_centered_expression=NA,
        num_obs_groups=1,
        output_dir=NULL,
        output_format="png",
        ## reference_observations=NULL,
        num_ref_groups=NULL,
        name_ref_groups=NULL,
        cluster_by_groups=FALSE,
        ref_subtract_method="by_mean",
        hclust_method="complete",
        clustering_contig=NULL,
        plot_steps=FALSE,
        bound_method_vis="average_bound",
        bound_threshold_vis=NULL, # example:  "-1,1"
        window_length=101,
        contig_tail=NULL,
        fig_title="Copy Number Variation Inference",
        obs_title="Observations (Cells)",
        ref_title="References (Cells)",
        load_workspace=2,
        log_level="INFO",
        ngchm=FALSE,
        path_to_shaidyMapGen=NULL,
        gene_symbol=NULL,
        min_cells_per_gene=3,
        use_zscores=FALSE,
        make_zero_NA=FALSE
        )
{

    if (!requireNamespace("fastcluster", quietly=TRUE)) {
        warning("fastcluster library not available, using the default hclust method instead.")
    }

    C_VIS_OUTLIER_CHOICES <- c("average_bound")
    C_REF_SUBTRACT_METHODS <- c("by_mean", "by_quantiles")
    C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
    C_OUTPUT_FORMAT <- c("pdf", "png")

    if (!(log_level %in% names(loglevels))) {
        stop("Error, log_level given not accepable, please use one from the following: ", paste(names(loglevels), collapse=", "))
    }
    logging::basicConfig()
    logging::setLevel(log_level, logging::getHandler('basic.stdout'))
    logging::basicConfig(level=log_level)
    if ( (output_dir == "") || (is.null(output_dir) || (!is.character(output_dir))) ) {
        stop("Error, no output_dir given. Please enter a file path to save the heatmap.")
    }
    if ( (!is.numeric(cutoff)) || (cutoff < 0) ) {
        stop("Error, cutoff < 0. Please enter a value greater or equal to zero.")
    }
    if ( !(bound_method_vis %in% C_VIS_OUTLIER_CHOICES)) {
        stop("Error, must specify acceptable vis_bound_method from one of the following: ", paste(C_VIS_OUTLIER_CHOICES, collapse=", "))
    }
    if ( !(ref_subtract_method %in% C_REF_SUBTRACT_METHODS)) {
        stop("Error, must specify accepable ref_subtract_method from one of the following: ", paste(C_REF_SUBTRACT_METHODS, collapse=", "))
    }
    if ( !(hclust_method %in% C_HCLUST_METHODS)) {
        stop("Error, must specify acceptable hclust_method from one of the following: ", paste(C_HCLUST_METHODS, collapse=", "))
    }
    if ( !(output_format %in% C_OUTPUT_FORMAT)) {
        stop("Error, must specify acceptable output_format from one of the following: ", paste(C_OUTPUT_FORMAT, collapse=", "))
    }

    if(!file.exists(output_dir)){
        dir.create(output_dir)
    }

    max_centered_expression <- abs(max_centered_expression)
    noise_filter <- abs(noise_filter)

    if (window_length < 0) {
        stop("Error, please provide a value greater or equal to zero for the window_length.")
    }
    if (is.null(contig_tail)) {
        contig_tail <- (window_length - 1)/2
    }
    if(contig_tail < 0) {
        stop("Error, please provide a value greater or equal to zero for the window_length.")
    }

    if (is.null(annotations)) {
        warning("Warning, no annotations provided.")
        if (!is.null(name_ref_groups)) {
            stop("Error, provided name_ref_groups but no annotations to read the grouping information from.")
        }
    }
    else {
        if (!is.null(name_ref_groups)) {
            name_ref_groups <- unlist(strsplit(name_ref_groups,","))  ## TOCHECK Could maybe require a list directly instead of a string with comas in between names?
        }
    }

    if (ngchm == TRUE){ ## check if required java application ShaidyMapGen.jar exists.
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
    }

    bounds_viz <- c(NA,NA)
    if (! (is.null(bound_threshold_vis) || is.na(bound_threshold_vis) ) ) {
        bounds_viz <- as.numeric(unlist(strsplit(bound_threshold_vis,",")))  ## TOCHECK maybe require a list directly instead of a string to split?
    }

    # do_work = 0

    passed_args <- list(x=x, gene_order=gene_order, annotations=annotations, cutoff=cutoff,
        log_transform=log_transform, delim=delim, noise_filter=noise_filter,
        max_centered_expression=max_centered_expression, num_obs_groups=num_obs_groups, num_ref_groups=num_ref_groups,
        name_ref_groups=name_ref_groups, ref_subtract_method=ref_subtract_method,
        hclust_method=hclust_method, clustering_contig=clustering_contig, window_length=window_length,
        contig_tail=contig_tail, bound_method_vis=bound_method_vis, bounds_viz=bounds_viz)

    to_save_processed_args <- c("x", "gene_order", "annotations", "cutoff", "log_transform", "delim", "noise_filter",
        "max_centered_expression", "num_obs_groups", "num_ref_groups", "name_ref_groups", "ref_subtract_method",
        "hclust_method", "clustering_contig", "window_length", "contig_tail", "bound_method_vis", "bounds_viz")

    processed_save_path <- paste(output_dir, "/infercnv.processed.", hclust_method, ".Rdata", sep="")
    processed_args_save_path <- paste(output_dir, "/infercnv.processed.args.", hclust_method, ".Rdata", sep="")
    if (file.exists(processed_args_save_path) && file.exists(processed_save_path) && load_workspace >= 2) {
        load(processed_args_save_path)
        if (identical(passed_args$x, x) &&
            identical(passed_args$gene_order, gene_order) &&
            identical(passed_args$annotations, annotations) &&
            identical(passed_args$cutoff, cutoff) &&
            identical(passed_args$log_transform, log_transform) &&
            identical(passed_args$delim, delim) &&
            identical(passed_args$noise_filter, noise_filter) &&
            identical(passed_args$max_centered_expression, max_centered_expression) &&
            identical(passed_args$num_obs_groups, num_obs_groups) &&
            identical(passed_args$num_ref_groups, num_ref_groups) &&
            identical(passed_args$name_ref_groups, name_ref_groups) &&
            identical(passed_args$ref_subtract_method, ref_subtract_method) &&
            identical(passed_args$hclust_method, hclust_method) &&  # not really needed since the filename implies it is the same hclust method
            identical(passed_args$clustering_contig, clustering_contig) &&
            identical(passed_args$window_length, window_length) &&
            identical(passed_args$contig_tail, contig_tail) &&
            identical(passed_args$bound_method_vis, bound_method_vis) &&
            identical(passed_args$bounds_viz, bounds_viz)
            ) {

            load(processed_save_path)
            logging::loginfo("Successfully loaded previously processed workspace.")
            # load_workspace <- 2 # only plot
        } else {  # restore arguments
            x <- passed_args$x
            gene_order <- passed_args$gene_order
            annotations <- passed_args$annotations
            cutoff <- passed_args$cutoff
            log_transform <- passed_args$log_transform
            delim <- passed_args$delim
            noise_filter <- passed_args$noise_filter
            max_centered_expression <- passed_args$max_centered_expression
            num_obs_groups <- passed_args$num_obs_groups
            num_ref_groups <- passed_args$num_ref_groups
            name_ref_groups <- passed_args$name_ref_groups
            ref_subtract_method <- passed_args$ref_subtract_method
            hclust_method <- passed_args$hclust_method
            clustering_contig <- passed_args$clustering_contig
            window_length <- passed_args$window_length
            contig_tail <- passed_args$contig_tail
            bound_method_vis <- passed_args$bound_method_vis
            bounds_viz <- passed_args$bounds_viz

            load_workspace <- 1 # try to load preprocessed data
        }
    } else if (load_workspace >= 2) {
        load_workspace <- 1
    }

    to_save_preprocess_args <- c("x", "gene_order", "annotations", "delim", "num_obs_groups", "num_ref_groups", "name_ref_groups")

    preprocess_save_path <- paste(output_dir, "/infercnv.preprocess.Rdata", sep="")
    preprocess_args_save_path <- paste(output_dir, "/infercnv.preprocess.args.Rdata", sep="")
    if (load_workspace == 1) {
        preprocess_save_path <- paste(output_dir, "/infercnv.preprocess.Rdata", sep="")
        if (file.exists(preprocess_args_save_path) && file.exists(preprocess_save_path)) {
            load(preprocess_args_save_path)
            if (identical(passed_args$x, x) &&
                identical(passed_args$gene_order, gene_order) &&
                identical(passed_args$annotations, annotations) &&
                identical(passed_args$delim, delim) &&
                identical(passed_args$num_obs_groups, num_obs_groups) &&
                identical(passed_args$num_ref_groups, num_ref_groups) &&
                identical(passed_args$name_ref_groups, name_ref_groups)
                ) {

                load(preprocess_save_path)
                logging::loginfo("Successfully loaded previously prepared workspace for processing.")
                # load_workspace <- 1  # process and plot
            } else {  # restore arguments
                x <- passed_args$x
                gene_order <- passed_args$gene_order
                annotations <- passed_args$annotations
                delim <- passed_args$delim
                num_obs_groups <- passed_args$num_obs_groups
                num_ref_groups <- passed_args$num_ref_groups
                name_ref_groups <- passed_args$name_ref_groups

                load_workspace <- 0 # can't usefully load any data
            }
        } else {
            load_workspace <- 0
        }
    }

    if (load_workspace <= 0) {
        if (!(is.data.frame(x))) {
            expression_data <- read.table(x, sep=delim, header=TRUE, row.names=1, check.names=FALSE)
        }
        else {
            expression_data <- x
        }

        # remove rows with NAs
        expression_data <- expression_data[complete.cases(expression_data), , drop=FALSE]

        input_gene_order <- seq(1, nrow(expression_data), 1)

        if (!.invalid(gene_order)) {
            if(!is.data.frame(gene_order)) {
                input_gene_order <- read.table(gene_order, header=FALSE, row.names=1, sep="\t")
                names(input_gene_order) <- c(CHR, START, STOP)
            }
            else {
                if (names(input_gene_order) != c(CHR, START, STOP)) {
                    stop("Error, gene_order given but names do not match expected format.")
                }
            }
        }

        input_reference_samples <- colnames(expression_data)
        observations_annotations_groups <- NULL

        if (!is.null(annotations)) {
            if ( is.data.frame(annotations)) {
                input_classifications <- annotations
            }
            else {
                input_classifications <- read.table(annotations, header=FALSE, row.names=1, sep=delim, stringsAsFactors=FALSE)
            }
            input_classifications <- input_classifications[order(match(row.names(input_classifications), colnames(expression_data))), , drop=FALSE]
            name_ref_groups_indices <- list()
            refs <- c()
            for (name_group in name_ref_groups) {
                name_ref_groups_indices[length(name_ref_groups_indices) + 1] <- list(which(input_classifications[,1] == name_group))
                refs <- c(refs, row.names(input_classifications[which(input_classifications[,1] == name_group), , drop=FALSE]))
            }
            input_reference_samples <- unique(refs)

            all_annotations <- unique(input_classifications[,1])
            observations_annotations_names <- setdiff(all_annotations, name_ref_groups)
        }

        if (!is.null(num_ref_groups)) {
            name_ref_groups_indices <- c(num_ref_groups)
        }

        if (length(input_reference_samples) !=
            length(intersect(input_reference_samples, colnames(expression_data)))){
            missing_reference_sample <- setdiff(input_reference_samples,
                                                colnames(expression_data))
            error_message <- paste("Please make sure that all the reference sample",
                                   "names match a sample in your data matrix.",
                                   "Attention to: ",
                                   paste(missing_reference_sample, collapse=","))
            stop(error_message)
        }

        order_ret <- order_reduce(data=expression_data, genomic_position=input_gene_order)
        expression_data <- order_ret$expr
        input_gene_order <- order_ret$order

        if(is.null(expression_data)){
            error_message <- paste("None of the genes in the expression data",
                                   "matched the genes in the reference genomic",
                                   "position file. Analysis Stopped.")
            stop(error_message)
        }

        obs_annotations_groups <- input_classifications[,1]
        counter <- 1
        obs_annotations_names <- c()
        for (classification in observations_annotations_names) {
          obs_annotations_groups[which(obs_annotations_groups == classification)] <- counter
          obs_annotations_names[counter] = classification
          counter <- counter + 1
        }
        names(obs_annotations_groups) <- rownames(input_classifications)
        obs_annotations_groups <- obs_annotations_groups[input_classifications[,1] %in% observations_annotations_names]  # filter based on initial input in case some input annotations were numbers overlaping with new format and to remove references indexes
        obs_annotations_groups <- as.integer(obs_annotations_groups)  ## they should already all be integers because they are based on "counter"

        grouping_key_coln <- c()
        grouping_key_coln[1] <- floor(123/(max(nchar(obs_annotations_names)) + 4))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
        if (grouping_key_coln[1] < 1) {
            grouping_key_coln[1] <- 1
        }
        grouping_key_coln[2] <- floor(123/(max(nchar(name_ref_groups)) + 4))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
        if (grouping_key_coln[2] < 1) {
            grouping_key_coln[2] <- 1
        }

        logging::loginfo("Saving workspace")

        to_save_preprocess <- c("expression_data", "input_gene_order", "input_reference_samples", "name_ref_groups", "name_ref_groups_indices", "obs_annotations_groups", "obs_annotations_names", "grouping_key_coln", "num_obs_groups")
        save(list=to_save_preprocess_args, file=preprocess_args_save_path)
        save(list=to_save_preprocess, file=preprocess_save_path)
    }



    #if (all.equal(bounds_viz, c(NA,NA))) {

        # determine viz bounds in a data-driven way.
    #    d = as.numeric(as.matrix(expression_data))
    #    q = quantile(d[d>0], c(0.9))
    #    bounds_viz=c(-1*q, q)


        #if (length(bounds_viz) != 2){
        #stop(paste("Please use the correct format for the argument",
        #                   "bound_threshold_vis . Two numbers seperated",
        #                   "by a comma is expected (lowerbound,upperbound)",
        #                   ". As an example, to indicate that outliers are",
        #                   "outside of -1 and 1 give the following.",
        #                   "vis_bound_threshold=\"-1,1\""))
    #}




    if (load_workspace < 2) {  # 0 or 1
        ret_list = process_data(data=expression_data,
                               gene_order=input_gene_order,
                               cutoff=cutoff,
                               reference_obs=input_reference_samples,
                               transform_data=log_transform,
                               window_length=window_length,
                               max_centered_threshold=max_centered_expression,
                               noise_filter=noise_filter,
                               name_ref_groups=name_ref_groups,
                               num_ref_groups=name_ref_groups_indices,
                               obs_annotations_groups=obs_annotations_groups,
                               obs_annotations_names=obs_annotations_names,
                               grouping_key_coln=grouping_key_coln,
                               out_path=output_dir,
                               k_obs_groups=num_obs_groups,
                               plot_steps=plot_steps,
                               contig_tail=contig_tail,
                               method_bound_vis=bound_method_vis,
                               lower_bound_vis=bounds_viz[1],
                               upper_bound_vis=bounds_viz[2],
                               ref_subtract_method=ref_subtract_method,
                               hclust_method=hclust_method,
                               min_cells_per_gene=min_cells_per_gene,
                               use_zscores=use_zscores,
                               make_zero_NA=make_zero_NA)


        to_save_processed <- c("ret_list", "num_obs_groups", "obs_annotations_groups", "obs_annotations_names", "grouping_key_coln", "hclust_method", "input_gene_order", "name_ref_groups_indices")  # only save what will be needed for plotting and can not be changed without the processing having to be changed too

        save(list=to_save_processed_args, file=processed_args_save_path)
        save(list=to_save_processed, file=processed_save_path)
    }

    # Output data before viz outlier
    write.table(ret_list["PREVIZ"], sep=delim,
                file=file.path(output_dir,
                           "expression_pre_vis_transform.txt"))

    # Output data after viz outlier
    write.table(ret_list["VIZ"], sep=delim,
                file=file.path(output_dir,
                           "expression_post_viz_transform.txt"))

    plot_cnv(plot_data=ret_list[["VIZ"]],
             contigs=ret_list[["CONTIGS"]],
             k_obs_groups=num_obs_groups,
             obs_annotations_groups=obs_annotations_groups,
             obs_annotations_names=obs_annotations_names,
             grouping_key_coln=grouping_key_coln,
             cluster_by_groups=cluster_by_groups,
             reference_idx=ret_list[["REF_OBS_IDX"]],
             ref_contig=clustering_contig,
             contig_cex=contig_label_size,
             ref_groups=ret_list[["REF_GROUPS"]],
             name_ref_groups=name_ref_groups,
             out_dir=output_dir,
             color_safe_pal=use_color_safe,
             hclust_method=hclust_method,
             title=fig_title,
             obs_title=obs_title,
             ref_title=ref_title,
             output_format=output_format)

    if (ngchm) {

        if (!requireNamespace("NGCHM", quietly=TRUE)) {
            stop("The \"NGCHM\" library is required to use \"-ngchm=TRUE\" but it is not available.", .call=FALSE)
        }

        logging::loginfo("Creating NGCHM as infercnv.ngchm")
        Create_NGCHM(plot_data = ret_list[["VIZ"]],
                       path_to_shaidyMapGen = path_to_shaidyMapGen,
                       reference_idx = ret_list[["REF_OBS_IDX"]],
                       ref_index = name_ref_groups_indices,
                       location_data = input_gene_order,
                       out_dir = output_dir,
                       contigs = ret_list[["CONTIGS"]],
                       ref_groups = ret_list[["REF_GROUPS"]],
                       title = fig_title,
                       gene_symbol = gene_symbol)
    }
}



get_average_bounds = function (data) {

    lower_bound <- mean(apply(data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))

    return(c(lower_bound, upper_bound))

}
