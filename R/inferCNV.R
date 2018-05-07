#!/usr/bin/env Rscript

CHR = "chr"
START = "start"
STOP = "stop"

# Remove the average of the genes of the reference observations from all 
# observations' expression. Normalization by column.
#
# Args:
# average_data: Matrix containing the data to remove average
#               from (this includes the reference observations).
#               Row = Genes, Col = Cells.
# ref_observations: Indices of reference observations.
#                   Only these are used in the average.
# ref_groups: A list of vectors of indices refering to the
#             different groups of the reference indices.
#
# Returns:
# Expression with the average gene expression in the reference 
#          observations removed.
average_over_ref <- function(average_data,
                             ref_observations,
                             ref_groups){
    # r = genes, c = cells
    logging::loginfo(paste("::average_over_ref:Start", sep=""))
    # Max and min mean gene expression within reference groups.
    average_max <- NULL
    average_min <- NULL
    average_reference_obs <- average_data[,ref_observations, drop=FALSE]
    # Reference gene within reference groups
    for (ref_group in ref_groups){
        grp_average <- rowMeans(average_reference_obs[,ref_group,
                                                      drop=FALSE],
                                na.rm=TRUE)
        if(is.null(average_max)){
            average_max <- grp_average
        }
        if(is.null(average_min)){
            average_min <- grp_average
        }
        average_max <- pmax(average_max, grp_average)
        average_min <- pmin(average_min, grp_average)
    }
    # Remove the Max and min averages of the reference groups from the
    # For each gene.
    for(gene_i in 1:nrow(average_data)){
        current_col <- average_data[gene_i, ]
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
# steps: Vector of colors to change use in the palette
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
# row_count: Total number of rows
# col_count: Total number of columns
# row_seps: Vector of integers indices for row breaks
# col_seps: Vector of integer indices for column breaks
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
# average_data: Matrix containing data. Row = Genes, Col = Cells.
# ref_obs: Indices of reference obervations.
# num_groups: The number of groups to partition nodes in or a list
#                       of already partitioned indices.
#
# Returns:
# Returns a list of grouped reference observations given as
#             vectors of groups. These are indices relative to the reference
#             observations only, so a return 1 indicates the first reference
#             row, not the first row.
split_references <- function(average_data,
                             ref_obs,
                             num_groups){
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
    if (length(num_groups) == 1){
        num_groups <- unlist(num_groups)
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
        hc <- hclust(dist(average_reference_obs))

        split_groups <- cutree(hc, k=num_groups)
        split_groups <- split_groups[hc$order]
        # Keep the sort of the hclust
        for(cut_group in unique(split_groups)){
            group_idx <- which(split_groups == cut_group)
            ret_groups[[cut_group]] <- hc$order[group_idx]
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
# out_method: Method to remove outliers [(average_bound, NA (hard threshold))]
# lower_bound: Lower bound which identifies a measurement
#                        as an outlier. 
# upper_bound: Upper bound which identifies a measurement
#                        as an outlier.
# plot_step: True will plot this analysis step.
#
# Returns:
# Return data matrix with outliers removed
remove_outliers_norm <- function(data,
                                 out_method=NA,
                                 lower_bound=NA,
                                 upper_bound=NA,
                                 plot_step=NA){
    logging::loginfo(paste("::remove_outlier_norm:Start"))
    if(is.na(data) || is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        return(data)
    }
    if (is.na(lower_bound) || is.na(upper_bound)){
        if(is.na(out_method)){
            logging::loginfo("::remove_outlier_norm:WARNING outlier removal was not performed.")
            return(data)
        }
        if (out_method == "average_bound"){
            lower_bound <- mean(apply(data, 2,
                                      function(x) quantile(x)[[1]]))
            upper_bound <- mean(apply(data, 2,
                                      function(x) quantile(x)[[5]]))
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
            data[data < lower_bound] <- lower_bound
            data[data > upper_bound] <- upper_bound
            return(data)
        } else {
            logging::logerror(paste("::remove_outlier_norm:Error, please",
                                    "provide an approved method for outlier",
                                    "removal for visualization."))
            stop(991)
        }
    }

    # Hard threshold given bounds
    if(!is.na(plot_step)){
        pdf(plot_step, useDingbats=FALSE)
        boxplot(data)
        points(1:ncol(data), rep(lower_bound, ncol(data)),
               pch=19, col="orange")
        points(1:ncol(data), rep(upper_bound, ncol(data)),
               pch=19, col="orange")
        dev.off()
    }
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound
    return(data)
}

# Center data after smoothing. Center with in cells using median.
#
# Args:
# data_smoothed: Matrix to center.
#                          Row = Genes, Col = cells.
#
# Returns:
# Matrix that is median centered.
#             Row = Genes, Col = cells.
center_smoothed <- function(data_smoothed){

    logging::loginfo(paste("::center_smoothed:Start"))
    # Center within columns
    row_median <- apply(data_smoothed, 2, median)
    return(t(apply(data_smoothed, 1, "-", row_median)))
}

# Center data and threshold (both negative and postive values)
#
# Args:
# center_data: Matrix to center. Row = Genes, Col = Cells.
# threshold: Values will be required to be with -/+1 * 
#                      threshold after centering.
# Returns:
# Centered and thresholded matrix
center_with_threshold <- function(center_data, threshold){

    logging::loginfo(paste("::center_with_threshold:Start", sep=""))
    # Center data (automatically ignores zeros)
    center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
    # Cap values between threshold and -threshold and recenter
    center_data[center_data > threshold] <- threshold
    center_data[center_data < (-1 * threshold)] <- -1 * threshold
    center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
    return(center_data)
}

# Returns the color palette for contigs.
#
# Returns:
# Color Palette
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}


#' @title Infer CNV changes given a matrix of RNASeq counts. Output a pdf and matrix of final values.
#'
#' @param data: Expression matrix (genes X samples),
#'                 assumed to be log2(TPM+1) .
#' @param gene_order: Ordering of the genes (data's rows) 
#'                       according to their genomic location
#'                       To include all genes use 0.
#' @param cutoff: Cut-off for the average expression of genes to be 
#'                   used for CNV inference.
#' @param reference_obs: Column names of the subset of samples (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all samples will 
#'                          be the reference.
#' @param transform_data: Indicator to log2 + 1 transform
#' @param window_length: Length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#' @param max_centered_threshold: The maximum value a a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value.
#' @param noise_threshold: The minimum difference a value can be from the 
#'                            average reference in order for it not to be
#'                            removed as noise.
#' @param num_ref_groups: The number of reference groups of a list of
#'                           indicies for each group of reference indices in
#'                           relation to reference_obs.
#' @param out_path: The path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
#' @param plot_steps: If true turns on plotting intermediate steps.
#' @param contig_tail: Length of the tail removed from the ends of contigs.
#' @param method_bound: Method to use for bounding values in the visualization.
#' @param lower_bound_vis: Lower bound to normalize data to for visualization.
#' @param upper_bound_vis: Upper bound to normalize data to for visualization.
#'
#' @return
#' Returns a list including:
#'     CNV matrix before visualization.
#'     CNV matrix after outlier removal for visualization.
#'     Contig order
#'     Column names of the subset of samples that should be used as references.
#'     Names of samples in reference groups. 
#' @export
infer_cnv <- function(data,
                      gene_order,
                      cutoff,
                      reference_obs,
                      transform_data,
                      window_length,
                      max_centered_threshold,
                      noise_threshold,
                      num_ref_groups,
                      out_path,
                      plot_steps=FALSE,
                      contig_tail= (window_length - 1) / 2,
                      method_bound_vis=NA,
                      lower_bound_vis=NA,
                      upper_bound_vis=NA){

    logging::loginfo(paste("::infer_cnv:Start", sep=""))
    # Return list of matrices
    ret_list <- list()

    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                  plot_name=file.path(out_path,
                                      "00_reduced_data.pdf"))
    }

    # Make sure data is log transformed + 1
    if (transform_data){
        data <- log2(data + 1)
    }
    
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                  plot_name=file.path(out_path,
                                      "02_transformed.pdf"))
    }

    # Reduce by cutoff
    keep_gene_indices <- above_cutoff(data, cutoff)
    if (!is.null(keep_gene_indices)){
        data <- data[keep_gene_indices, , drop=FALSE]
        gene_order <- gene_order[keep_gene_indices, , drop=FALSE]
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff, ",
                      "new dimensions (r,c) = ",
                      paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
        logging::logdebug(paste("::infer_cnv:Keeping indices.", sep=""))
    } else {
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff.", sep=""))
        logging::logwarn(paste("::No indicies left to keep.",
                               " Stoping."))
        stop(998)
    }
    
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(out_path,
                                                "03_reduced_by_cutoff.pdf"))
    }

    # Reduce contig info
    chr_order <- gene_order[1]
    gene_order <- NULL

    # Center data (automatically ignores zeros)
    data <- center_with_threshold(data, max_centered_threshold)
    logging::loginfo(paste("::infer_cnv:Outlier removal, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(out_path,
                                                "05_center_with_threshold.pdf"))
    }

    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_length)
    data <- NULL
    logging::loginfo(paste("::infer_cnv:Smoothed data.", sep=""))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "06_smoothed.pdf"))
    }

    # Center cells/observations after smoothing. This helps reduce the
    # effect of complexity.
    data_smoothed <- center_smoothed(data_smoothed)
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "07_recentered.pdf"))
    }

    # Split the reference data into groups if requested
    groups_ref <- split_references(average_data=data_smoothed,
                                             ref_obs=reference_obs,
                                             num_groups=num_ref_groups)
    logging::loginfo(paste("::infer_cnv:split_reference. ",
                           "found ",length(groups_ref)," reference groups.",
                           sep=""))

    # Remove average reference
    i_ref_obs <- which(colnames(data_smoothed) %in% reference_obs)
    data_smoothed <- average_over_ref(average_data=data_smoothed,
                                                ref_observations=i_ref_obs,
                                                ref_groups=groups_ref)
    logging::loginfo(paste("::infer_cnv:Remove average, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "08_remove_average.pdf"))
    }

    # Remove Ends
    logging::logdebug(chr_order)
    remove_indices <- c()
    for (chr in unlist(unique(chr_order))){
        logging::loginfo(paste("::infer_cnv:Remove tail contig ",
                               chr, ".", sep=""))
        remove_chr <- remove_tails(data_smoothed,
                                             which(chr_order == chr),
                                             contig_tail)
        remove_indices <- c(remove_indices, remove_chr)
    }
    if (length(remove_indices) > 0){
        chr_order <- chr_order[remove_indices,]
        data_smoothed <- data_smoothed[remove_indices, , drop=FALSE]
    }
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "09_remove_ends.pdf"))
    }
    logging::loginfo(paste("::infer_cnv:Remove ends, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))

    # Remove noise
    data_smoothed <- remove_noise(smooth_matrix=data_smoothed,
                                            threshold=noise_threshold)
    logging::loginfo(paste("::infer_cnv:Remove moise, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(out_path,
                                                "11_denoise.pdf"))
    }

    # Output before viz outlier
    ret_list[["PREVIZ"]] = data_smoothed

    # Remove outliers for viz
    remove_outlier_viz_pdf <- NA
    if (plot_steps){
        remove_outlier_viz_pdf <- file.path(out_path,
                                            "10A_remove_outlier.pdf")
    }
    ret_list[["VIZ"]] <- remove_outliers_norm(data=data_smoothed,
                                          out_method=method_bound_vis,
                                          lower_bound=lower_bound_vis,
                                          upper_bound=upper_bound_vis,
                                          plot_step=remove_outlier_viz_pdf)
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=ret_list[["VIZ"]],
                  plot_name=file.path(out_path,
                                      "10B_remove_outlier.pdf"))
    }
    logging::loginfo(paste("::infer_cnv:remove outliers, ",
                           "new dimensions (r,c) = ",
                           paste(dim(ret_list[["VIZ"]]), collapse=","),
                           " Total=", sum(ret_list[["VIZ"]]),
                           " Min=", min(ret_list[["VIZ"]]),
                           " Max=", max(ret_list[["VIZ"]]),
                           ".", sep=""))

    ret_list[["CONTIGS"]] = paste(as.vector(as.matrix(chr_order)))
    ret_list[["REF_OBS_IDX"]] = reference_obs
    ret_list[["REF_GROUPS"]] = groups_ref
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
    write.table(data, file=text_file)
}

#' @title Plot the matrix as a heatmap. Clustering is on observation only, gene position is preserved.
#'
#' @param plot_data: Data matrix to plot (columns are observations).
#' @param contigs: The contigs the data is group in in order of rows.
#' @param reference_idx: Vector of reference indices.
#' @param ref_contig: If given, will focus cluster on only genes in this contig
#' @param reg_groups: Groups of vector indices (as indices in reference_idx)
#' @param out_dir: Directory in which to save pdf and other output.
#' @param title: Plot title.
#' @param obs_title: Title for the observations matrix.
#' @param ref_title: Title for the reference matrix.
#' @param contig_cex: Contig text size.
#' @param k_obs_groups: Number of groups to break observation into
#' @param color_safe_pal: Logical indication of using a color blindness safe
#'                          palette.
#'
#' @return
#' No return, void.
#' @export
plot_cnv <- function(plot_data,
                     contigs,
                     reference_idx,
                     ref_contig,
                     ref_groups,
                     out_dir,
                     title,
                     obs_title,
                     ref_title,
                     contig_cex=1,
                     k_obs_groups=1,
                     color_safe_pal=TRUE){

    logging::loginfo(paste("::plot_cnv:Start", sep=""))
    logging::loginfo(paste("::plot_cnv:Current data dimensions (r,c)=",
                            paste(dim(plot_data), collapse=","),
                            " Total=", sum(plot_data),
                            " Min=", min(plot_data),
                            " Max=", max(plot_data),
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
    col_sep <- col_sep[-1 * length(col_sep)]
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

    # Rows observations, Columns CHR
    pdf(paste(out_dir,"infercnv.pdf",sep="/"),
        useDingbats=FALSE,
        width=10,
        height=7.5,
        paper="USr")

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
    ref_data_t <- plot_data[, ref_idx, drop=FALSE]
    nb_breaks <- 16
    breaksList_t <-
        seq(min(min(obs_data_t, na.rm=TRUE), min(ref_data_t, na.rm=TRUE)),
        max(max(obs_data_t,na.rm=TRUE), max(ref_data_t, na.rm=TRUE)),
        length.out=nb_breaks)

    # Create file base for plotting output
    force_layout <- plot_observations_layout()
    plot_cnv_observations(obs_data=obs_data_t,
                          file_base_name=out_dir,
                          cluster_contig=ref_contig,
                          contig_colors=ct.colors[contigs],
                          contig_label=contig_labels,
                          contig_names=contig_names,
                          col_pal=custom_pal,
                          contig_seps=col_sep,
                          num_obs_groups=k_obs_groups,
                          cnv_title=title,
                          cnv_obs_title=obs_title,
                          contig_lab_size=contig_cex,
                          breaksList=breaksList_t,
                          layout_lmat=force_layout[["lmat"]],
                          layout_lhei=force_layout[["lhei"]],
                          layout_lwid=force_layout[["lwid"]])
    obs_data <- NULL

    if(!is.null(ref_idx)){
        plot_cnv_references(ref_data=ref_data_t,
                            ref_groups=ref_groups,
                            col_pal=custom_pal,
                            contig_seps=col_sep,
                            file_base_name=out_dir,
                            cnv_ref_title=ref_title,
                            breaksList=breaksList_t,
                            layout_add=TRUE)
    }
    dev.off()
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the observational samples
#
# Args:
# obs_data: Data to plot as observations. Rows = Cells, Col = Genes
# col_pal: The color palette to use.
# contig_colors: The colors for the contig bar.
# contig_labels: The labels for the contigs.
# contig_names: Names of the contigs
# contig_seps: Indices for line seperators of contigs.
# num_obs_groups: Number of groups of observations to create
# file_base_name: Base of the file to used to make output file names.
# cnv_title: Title of the plot.
# cnv_obs_title: Title for the observation matrix.
# contig_lab_size: Text size for contigs.
# cluster_contig: A value directs cluster to only genes on this contig
# layout_lmat: lmat values to use in layout
# layout_lhei: lhei values to use in layout
# layout_lwid: lwid values to use in layout
#
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
                                  breaksList,
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
           logging::logwarning(paste("plot_cnv_observations: Not able to cluster by",
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
    obs_hcl <- hclust(dist(obs_data[,hcl_group_indices]),"average")
    obs_dendrogram <- as.dendrogram(obs_hcl)
    write.tree(as.phylo(obs_hcl),
               file=paste(file_base_name, "observations_dendrogram.txt", sep=.Platform$file.sep))

    # Output HCL group membership.
    # Record locations of seperations
    obs_seps <- c(0)
    ordered_names <- rev(row.names(obs_data)[obs_hcl$order])
    split_groups <- cutree(obs_hcl, k=num_obs_groups)

    # Make colors based on groupings
    row_groupings <- get_group_color_palette()(length(table(split_groups)))[split_groups]
    
    # Make a file of coloring and groupings
    logging::loginfo("plot_cnv_observation:Writing observation groupings/color.")
    groups_file_name <- file.path(file_base_name, "observation_groupings.txt")
    file_groups <- rbind(split_groups,row_groupings)
    row.names(file_groups) <- c("Group","Color")
    write.table(t(file_groups), groups_file_name)

    # Make a file of members of each group
    logging::loginfo("plot_cnv_observation:Writing observations by grouping.")
    for (cut_group in unique(split_groups)){
        group_memb = names(split_groups)[which(split_groups == cut_group)]
        # Write group to file
        memb_file <- file(paste(file_base_name,
                                paste(hcl_desc,"HCL",cut_group,"members.txt",sep="_"),
                                sep=.Platform$file.sep))
        write.table(obs_data[group_memb,], memb_file)
        # Record seperation
        ordered_memb <- which(ordered_names %in% group_memb)
        obs_seps <- c(obs_seps, max(ordered_memb),max(ordered_memb))
    }
    obs_hcl <- NULL
    obs_seps <- unique(obs_seps)
    obs_seps <- sort(obs_seps)
    
    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(obs_data),
                                     col_count=ncol(obs_data),
                                     col_seps=contig_seps)
    

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
                                        x.center=0,
                                        color.FUN=col_pal,
                                        if.plot=!testing,
                                        # Seperate by contigs
                                        sepList=contigSepList,
                                        sep.color="black",
                                        sep.lty=1,
                                        sep.lwd=1,
                                        # Color rows by user defined cut
                                        RowIndividualColors=row_groupings,
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
plot_observations_layout <- function()
{
    ## Plot observational samples
    obs_lmat <- c(0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                  6, 8, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                  0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  4, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    obs_lmat <- matrix(obs_lmat,ncol=13,byrow=TRUE)

    obs_lhei <- c(1.125, 2.215, .15,
                   .5, .5, .5,
                   .5, .5, .5,
                   .5, .5, .5,
                  0.0075, 0.0075, 0.0075)

    return(list(lmat=obs_lmat,
           lhei=obs_lhei,
           lwid=NULL))
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the reference samples
#
# Args:
# ref_data: Data to plot as references. Rows = Cells, Col = Genes
# ref_groups: Groups of references to plot together.
# col_pal: The color palette to use.
# contig_seps: Indices for line seperators of contigs.
# file_base_name: Base of the file to used to make output file names.
# cnv_ref_title: Title for reference matrix.
# layout_lmat: lmat values to use in the layout.
# layout_lwid: lwid values to use in the layout.
# layout_lhei: lhei values to use in the layout.
# layout_add: Indicates the ref image shoudl be added to the previous plot.
# testing: Turns off plotting when true.
#
# Returns:
# Void
plot_cnv_references <- function(ref_data,
                                ref_groups,
                                col_pal,
                                contig_seps,
                                file_base_name,
                                cnv_ref_title,
                                breaksList,
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

    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(ref_data),
                                     col_count=ncol(ref_data),
                                     row_seps=ref_seps,
                                     col_seps=contig_seps)

    # Print controls
    logging::loginfo("plot_cnv_references:Plotting heatmap.")
    data_references <- heatmap.cnv(ref_data,
                                   main=NA,
                                   ylab=reference_ylab,
                                   xlab=NA,
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
                                   x.center=0,
                                   color.FUN=col_pal,
                                   sepList=contigSepList,
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
# data: Data to measure the average row and evaluate
#                 against the cutoff. Row = Genes, Col = Cells.
# cutoff: Threshold to be above to be kept.
#
# Returns:
# Returns a vector of row indicies to keep (are above the cutoff).
above_cutoff <- function(data, cutoff){

    logging::loginfo(paste("::above_cutoff:Start", sep=""))
    average_gene <- log2(rowMeans( ( (2 ^ data) - 1), na.rm=TRUE) + 1 )
    logging::loginfo(paste("::infer_cnv:Averages (counts).", sep=""))
    # Find averages above a certain threshold
    indicies <- which(average_gene > cutoff)
    if (length(indicies) > 0){
        return(indicies)
    } else {
        return(NULL)
    }
}

#' Order the data and subset the data to data in the genomic position file.
#'
#' Args:
#' @param data: Data (expression) matrix where the row names should be in
#'                 the row names of the genomic_position file.
#' @param genomic_position: Data frame read in from the genomic position file
#'
#' @return Returns a matrix of expression in the order of the
#'            genomic_position file. NULL is returned if the genes in both
#'            data parameters do not match.
#' @export
order_reduce <- function(data, genomic_position){
    logging::loginfo(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL, chr_order=NULL)
    if (is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }

    # Drop pos_gen entries that are position 0
    remove_by_position <- -1 * which(genomic_position[2] + genomic_position[3] == 0)
    if (length(remove_by_position)){
        logging::logdebug(paste("::infer_cnv:order_reduce: removing genes specified by pos == 0, count: ",
                                length(remove_by_position), sep=""))
        
        genomic_position <- genomic_position[remove_by_position, , drop=FALSE]
    }

    # Reduce to genes in pos file
    
    logging::logdebug(paste("::infer_cnv:order_reduce: gene identifers in expression matrix: ",
                            row.names(data), collapse="\n", sep=""))
    logging::logdebug(paste("::infer_cnv:order_reduce: gene identifers in genomic position table: ",
                            row.names(data), collapse="\n", sep=""))


    
    keep_genes <- row.names(data)[which(row.names(data)
                                  %in% row.names(genomic_position))]
    logging::logdebug(paste("::infer_cnv:order_reduce: keep_genes size: ", length(keep_genes),
                            sep=""))
    
    # Keep genes found in position file
    if(length(keep_genes)){
        ret_results$expr <- data[keep_genes, , drop=FALSE]
        ret_results$order <- genomic_position[keep_genes, , drop=FALSE]
    } else {
        logging::loginfo(paste("::infer_cnv:order_reduce:The position file ",
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
    logging::loginfo(paste("::infer_cnv:order_reduce:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    logging::logdebug(paste("::infer_cnv:order_reduce end."))
    return(ret_results)
}

# Remove values that are too close to the average and are considered noise.
#
# Args:
# smooth_matrix: A matrix of values, smoothed, and with average 
#                          reference removed. Row = Genes, Col = Cells.
# threshold: The amount of difference a value must be from the
#                      reference before the value can be kept and not 
#                      removed as noise.
# Returns:
# Denoised matrix
remove_noise <- function(smooth_matrix, threshold){

    logging::loginfo(paste("::remove_noise:Start.", sep=""))
    if (threshold > 0){
        smooth_matrix[abs(smooth_matrix) < threshold] <- 0
    }
    return(smooth_matrix)
}

# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix: Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr: Indices of the chr in which the tails are to be removed.
# tail_length: Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.
remove_tails <- function(smooth_matrix, chr, tail_length){

    logging::loginfo(paste("::remove_tails:Start.", sep=""))
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
# data: Data matrix to smooth. Row = Genes, Col = Cells.
# window_length: Length of window to use for the moving average.
#        Should be a positive, odd integer.
#
# Returns:
# Matrix with columns smoothed with a simple moving average.
smooth_window <- function(data, window_length){

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
    # Fix ends
    data_end <- apply(data,
                      2,
                      smooth_ends_helper,
                      obs_tails=tail_length)

    for (row_end in 1:tail_length){
        end_bound <- (num_genes - row_end) + 1
        data_sm[row_end, ] <- data_end[row_end, ]
        data_sm[end_bound, ] <- data_end[end_bound, ]
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
# obs_tails: Length of the tail to smooth.
#
# Returns:
# Data smoothed.
smooth_ends_helper <- function(obs_data, obs_tails){
    end_data <- rep(NA,length(obs_data))
    obs_count <- length(obs_data)
    for (tail_end in 1:obs_tails){
        bounds <- tail_end - 1
        end_tail <- obs_count - bounds
        end_data[tail_end] <- mean(obs_data[(tail_end - bounds):
                                           (tail_end + bounds)],
                                           na.rm=TRUE)
        end_data[end_tail] <- mean(obs_data[(end_tail - bounds):
                                           (end_tail + bounds)],
                                           na.rm=TRUE)
    }
    return(end_data)
}

# Smooth vector of values over the given window length.
#
# Args:
# obs_data: Vector of data to smooth with a moving average.
# window_length: Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.
smooth_window_helper <- function(obs_data, window_length){

    return(filter(obs_data, rep(1 / window_length, window_length), sides=2))
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
           hclust.FUN.MoreArgs=list(method="ward"),

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
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexRow.')
    cexRow <- cexRow0*cexRow
  }
  if (.invalid(cexCol)) {
    cexCol <- cexCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for cexCol.')
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
    rowInd <- order.dendrogram(ddr)
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
      margin.for.labRow = 2
  }

  if (.invalid(margin.for.labCol)){
    margin.for.labCol <- margin.for.labCol0
  } else {
    message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labCol.')
    margin.for.labCol <- margin.for.labCol0*margin.for.labCol
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
  }

  ## get x.range according to the value of x.center ##
  if (!.invalid(x.center)){ ## enhanced
    if (is.numeric(x.center)){
      x.range.old <- range(x,na.rm=TRUE)
      dist.to.x.center <- max(abs(x.range.old-x.center))
      x.range <- c(x.center-dist.to.x.center,x.center+dist.to.x.center)
    } else {
      stop("`x.center' should be numeric.")
    }
  } else{
    x.range <- range(x,na.rm=TRUE)
  }


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
      if(!is.character(RowIndividualColors) || length(RowIndividualColors) !=nr)
        stop("'RowIndividualColors' must be a character vector of length nrow(x)")
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)
      lwid <- c(0.2,lwid)
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
      image(rbind(1:nr),col=RowIndividualColors[rowInd],axes=FALSE,add=force_add)
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
    ## 7) col-dend and title
    mar3 <- (if(!is.null(main)) mar.main else 0) +
      (if(!is.null(sub)) mar.sub else 0)
    par(mar=c(0,0,mar3,margins[4]))

    if(dendrogram %in% c("both","column")){
      plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    } else{
      .plot.text(xlim=range(1:nc))
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
