#!/usr/bin/env Rscript

CHR = "chr"
START = "start"
STOP = "stop"

#' Remove the average of the genes of the reference observations from all 
#' observations' expression. Normalization by column.
#'
#' Args:
#'    @param average_data: Matrix containing the data to remove average from
#'                         (this includes the reference observations).
#'                         Row = Genes, Col = Cells.
#'    @param reference_observations: Indices of reference observations.
#'                                   Only these are used in the average.
#'    @param ref_groups: A list of vectors of indices refering to the
#'                       different groups of the reference indices.
#'
#' Returns:
#'    @return Expression with the average gene expression in the reference 
#'            observations removed.
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

#' Helper function allowing greater control over the steps in a color palette.
#' Source:http://menugget.blogspot.com/2011/11/define-color-steps-for-
#'               colorramppalette.html#more
#'
#' Args:
#'    @param steps: Vector of colors to change use in the palette
#'    @param between: Steps where gradients change
#' Returns:
#'    @returns: Color palette
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

#'
#'
#'
#create_column_plot <- function(data, col_pal, col_height=5){
#    print("dims")
#    print(dim(data))
#    # Get means and adjust them to a positive range of integers
#    # For the color palette.
#    col_means <- colMeans(data)
#    max_cnv <- max(max(col_means), abs(min(col_means)))
#    relative_values <- round(col_means / max_cnv * 100)
#    min_rel <- min(relative_values)
#    if(min_rel < 0){
#        relative_values <- relative_values + abs(min_rel)
#    }
#
#    # Make a color palette
#    col_color_pal <- col_pal(max(relative_values)+1)
#    cnv_colors <- col_color_pal[relative_values]
#    cnv_colors[col_means == 0] <- "white"
#
#    # Make a number of columns to simulate a bar graph given the
#    # magnitude of the number.
#    col_columns <- list()
#    col_increment <- max(abs(col_means))/col_height
#    print("max(abs(col_means))")
#    print(max(abs(col_means)))
#    for(height in 1:col_height){
#        cur_col_increment <- col_increment * ( height - 1 )
#        print("cur_col_increment")
#        print(cur_col_increment)
#        cnv_colors[col_means > 0 & col_means < cur_col_increment] <- "white"
#        cnv_colors[col_means < 0 & col_means > (-1 * cur_col_increment)] <- "white"
#        col_columns[[height]] <- cnv_colors
#    }
#    ret_col_info <- as.matrix(as.data.frame(col_columns, stringsAsFactors=FALSE))
#    colnames(ret_col_info) <- rep(".", ncol(ret_col_info))
#    print(ret_col_info)
#    return(ret_col_info)
#}

#' Create a sepList forthe heatmap.3 plotting function given integer vectors
#' of rows and columns where speration should take place.
#' The expected input to thebheatmap function is a list of 2 lists.
#' The first list are column based rectangles, and the second row.
#' To define a rectagle the index of the row or column where the line of the rectagle
#' should be placed is done with a vector of integers, left, bottom, right and top line.
#' Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
create_sep_list <- function(row_count, col_count,
                            col_seps=NULL,
                            row_seps=NULL){
    sepList <- list()
    # Heatmap.3 wants a list of boxes for seperating columns
    # Column data
    if(!is.null(col_seps) && !is.na(col_seps) && (length(col_seps)>0)){
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
    if(!is.null(row_seps) && !is.na(row_seps) && (length(row_seps)>0)){
        rowList <- list()
        for(sep in 1:length(row_seps)){
            rowList[[sep]] <- c(row_seps[sep],0,row_seps[sep],col_count)
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}

#' Split up reference observations in to k groups and return indices
#' for the different groups.
#'
#' Args:
#'    @param average_data: Matrix containing data. Row = Genes, Col = Cells.
#'    @param ref_obs: Indices of reference obervations.
#'    @param num_groups: The number of groups to partition nodes in or a list
#'                       of already partitioned indices.
#'
#' Returns:
#'    @return Returns a list of grouped reference observations given as
#'            vectors of groups. These are indices relative to the reference
#'            observations only, so a return 1 indicates the first reference
#'            row, not the first row.
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

#' Set outliers to some upper or lower bound. Then normalize values to
#' approximately [-1, 1]. This is to prep the data for visualization.
#'
#' Args:
#'    @param data: data to remove outliers. Outliers removed within columns.
#'    @param lower_bound: Lower bound which identifies a measurement
#'                        as an outlier. 
#'    @param upper_bound: Upper bound which identifies a measurement
#'                        as an outlier.
remove_outliers_norm <- function(data,
                                 out_method=NA,
                                 lower_bound=NA,
                                 upper_bound=NA,
                                 plot_step=NA){
    logging::loginfo(paste("::remove_outlier_norm:Start"))
    if (is.na(lower_bound) || is.na(upper_bound)){
        if(is.na(out_method)){
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

#' Center data after smoothing. Center with in cells using median.
#'
#' Args:
#'    @param data_smoothed: Matrix to center.
#'                          Row = Genes, Col = cells.
#'
#' Returns:
#'    @return: Matrix that is median centered.
#'             Row = Genes, Col = cells.
center_smoothed <- function(data_smoothed){

    logging::loginfo(paste("::center_smoothed:Start"))
    # Center within columns
    row_median <- apply(data_smoothed, 2, median)
    return(t(apply(data_smoothed, 1, "-", row_median)))
}

#' Center data and threshold (both negative and postive values)
#'
#' Args:
#'    @param center_data: Matrix to center. Row = Genes, Col = Cells.
#'    @param threshold: Values will be required to be with -/+1 * 
#'                      threshold after centering.
#' Returns:
#'    @return Centered and thresholded matrix
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

#' Check arguments and make sure the user input meet certain 
#' additional requirements.
#'
#' Args:
#'    @param arguments: Parsed arguments from user.
#' Returns:
#'    @return No return.
check_arguments <- function(arguments){

    logging::loginfo(paste("::check_arguments:Start", sep=""))
    # Require the name of a output pdf file
    if (!( "pdf_file" %in% names(arguments)) || (arguments$pdf_file == "")){
        logging::logerror(paste(":: --pdf: Please enter a file path to ",
                                "save the heatmap.",
                                 sep=""))
    }

    # Require the cut off to be above 0
    if (arguments$cutoff < 0){
        logging::logerror(paste(":: --cutoff: Please enter a value",
                                "greater or equal to zero for the cut off.",
                                sep=""))
    }

    # Require the logging level to be one handled by logging
    if (!(arguments$log_level %in% C_LEVEL_CHOICES)){
        logging::logerror(paste(":: --log_level: Please use a logging level ",
                                "given here: ", C_LEVEL_CHOICES,
                                collapse=",", sep=""))
    }

    # Require the visualization outlier detection to be a correct choice.
    if (!(arguments$bound_method_vis %in% C_VIS_OUTLIER_CHOICES)){
        logging::logerror(paste(":: --vis_bound_method: Please use a method ",
                                "given here: ", C_VIS_OUTLIER_CHOICES,
                                collapse=",", sep=""))
    }

    # Warn that an average of the samples is used in the absence of
    # normal / reference samples
    if (is.null(arguments$reference_observations)){
        logging::logwarn(paste(":: --reference_observations: No reference ",
                      "samples were given, the average of the samples ",
                      "will be used.",
                      sep=""))
    }

    # Make sure the threshold is centered.
    arguments$max_centered_expression <- abs(arguments$max_centered_expression)
    arguments$magnitude_filter <- abs(arguments$magnitude_filter)

    # Require the contig tail to be above 0
    if (is.na(arguments$contig_tail)){
        arguments$contig_tail <- (arguments$window_length - 1) / 2
    }

    if (arguments$contig_tail < 0){
        logging::logerror(paste(":: --tail: Please enter a value",
                                "greater or equal to zero for the tail.",
                                sep=""))
    }

    if (! is.na(suppressWarnings(as.integer(arguments$num_groups)))){
        arguments$num_groups <- list(as.integer(arguments$num_groups))
    } else {
        # Warn references must be given.
        if (is.null(arguments$reference_observations)){
            logging::logerror(paste(":: --ref_groups to use this function ",
                                    "references must be given. "))
        }

        # TODO need to check and make sure all reference indices are given.
        num_str <- unlist(strsplit(arguments$num_groups,","))
        if (length(num_str) == 1){
            logging::logerror(paste(":: --ref_groups. If explicitly giving ",
                                    "indices, make sure to give atleast ",
                                    "two groups", sep =""))
            stop(990)
        }

        num_groups <- list()
        for (num_token in num_str){
            token_numbers <- unlist(strsplit(num_token, ":"))
            number_count <- length(token_numbers)
            if (number_count == 1){
                singleton <- as.integer(number_count)
                num_groups[[length(num_groups) + 1]] <- singleton
            } else if (number_count == 2){
                from <- as.integer(token_numbers[1])
                to <- as.integer(token_numbers[2])
                num_groups[[length(num_groups) + 1]] <- seq(from, to)
            } else {
                logging::logerror(paste(":: --ref_groups is expecting either ",
                                        "one number or a comma delimited list ",
                                        "of numbers or spans using ':'. ",
                                        "Examples include: --ref_groups 3 or ",
                                        " --ref_groups 1,3,5,6,3 or ",
                                        " --ref_groups 1:5,6:20 or ",
                                        " --ref_groups 1,2:5,6,7:10 .", sep=""))
                stop(999)
            }
        }
        arguments$num_groups <- num_groups
    }
    return(arguments)
}

#' Returns the color palette for contigs.
#'
#' Returns:
#'    @return Color Palette
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}

#' Infer CNV changes given a matrix of RNASeq counts.
#' Output a pdf and matrix of final values.
#'
#' Args:
#'    @param data: Expression matrix (genes X samples),
#'                 assumed to be log2(TPM+1) .
#'    @param gene_order: Ordering of the genes (data's rows) 
#'                       according to their genomic location
#'                       To include all genes use 0.
#'    @param cutoff: Cut-off for the average expression of genes to be 
#'                   used for CNV inference.
#'    @param reference_obs: Column names of the subset of samples (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all samples will 
#'                          be the reference.
#'    @param transform_data: Inidcator to log2 + 1 transform
#'    @param window_length: Length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#'    @param max_centered_threshold: The maximum value a a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value.
#'    @param noise_threshold: The minimum difference a value can be from the 
#'                            average reference in order for it not to be
#'                            removed as noise.
#'    @param num_ref_groups: The number of reference groups of a list of
#'                           indicies for each group of reference indices in
#'                           relation to reference_obs.
#'    @param pdf_path: The path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
#'    @param plot_steps: If true turns on plotting intermediate steps.
#'    @contig_tail: Length of the tail removed from the ends of contigs.
#'    @color_safe: Logical indicator to use a color blind safe palette (TRUE) or
#'                 the original publication color scheme.
#'    @lower_bound_vis: Lower bound to normalize data to for visualization.
#'    @upper_bound_vis: Upper bound to normalize data to for visualization.
#'
#' Returns:
#'    @return No return.
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
                      num_obs_groups,
                      pdf_path,
                      plot_steps=FALSE,
                      contig_tail= (window_length - 1) / 2,
                      color_safe=TRUE,
                      method_bound_vis=NA,
                      lower_bound_vis=NA,
                      upper_bound_vis=NA){

    logging::loginfo(paste("::infer_cnv:Start", sep=""))
    plot_steps_path <- dirname(pdf_path)
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(plot_steps_path,
                                                "00_reduced_data.pdf"))
    }

    # Remove any gene without position information
    # Genes may be sorted correctly by not have position information
    # Here they are removed.
    remove_by_position <- -1 * which(gene_order[2] + gene_order[3] == 0)
    if (length(remove_by_position)){
        gene_order <- gene_order[remove_by_position, , drop=FALSE]
        data <- data[remove_by_position, , drop=FALSE]
    }
    logging::loginfo(paste("::infer_cnv:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    logging::logdebug(paste("::infer_cnv:Removed indices:"))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(plot_steps_path,
                                                "01_remove_by_pos_file.pdf"))
    }

    # Make sure data is log transformed + 1
    if (transform_data){
        data <- log2(data / 10 + 1)
    }
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(plot_steps_path,
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
                            plot_name=file.path(plot_steps_path,
                                                "03_reduced_by_cutoff.pdf"))
    }

    # Order genes by genomic region
    data <- data[with(gene_order, order(chr,start,stop)), , drop=FALSE]
    # This is the contig order, will be used in visualization.
    # Get the contig order in the same order as the genes.
    chr_order <- gene_order[with(gene_order,
                                 order(chr,start,stop)), , drop=FALSE][1]
    gene_order <- NULL
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data,
                            plot_name=file.path(plot_steps_path,
                                                "04_order_by_chr.pdf"))
    }

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
                            plot_name=file.path(plot_steps_path,
                                                "05_center_with_threshold.pdf"))
    }

    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_length)
    data <- NULL
    logging::loginfo(paste("::infer_cnv:Smoothed data.", sep=""))
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(plot_steps_path,
                                                "06_smoothed.pdf"))
    }

    # Center cells/observations after smoothing. This helps reduce the
    # effect of complexity.
    data_smoothed <- center_smoothed(data_smoothed)
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(plot_steps_path,
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
                            plot_name=file.path(plot_steps_path,
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
                            plot_name=file.path(plot_steps_path,
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
                            plot_name=file.path(plot_steps_path,
                                                "11_denoise.pdf"))
    }

    # Output before viz outlier
    write.table(data_smoothed, file=paste(pdf_path, ".txt", sep=""))
    logging::loginfo(paste("::infer_cnv:Writing final data to ",
                           paste(pdf_path, ".txt", sep=""), sep=""))

    # Remove outliers for viz
    remove_outlier_viz_pdf <- NA
    if (plot_steps){
        remove_outlier_viz_pdf <- file.path(plot_steps_path,
                                            "10A_remove_outlier.pdf")
    }
    data_smoothed <- remove_outliers_norm(data=data_smoothed,
                                          out_method=method_bound_vis,
                                          lower_bound=lower_bound_vis,
                                          upper_bound=upper_bound_vis,
                                          plot_step=remove_outlier_viz_pdf)
    # Plot incremental steps.
    if (plot_steps){
        plot_step(data=data_smoothed,
                            plot_name=file.path(plot_steps_path,
                                                "10B_remove_outlier.pdf"))
    }
    logging::loginfo(paste("::infer_cnv:remove outliers, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))

    # Plot and write data
    logging::loginfo(paste("::infer_cnv:Drawing plots to file:",
                           pdf_path, sep=""))
    logging::loginfo(paste("::infer_cnv:Current data dimensions (r,c)=",
                           paste(dim(data_smoothed), collapse=","), sep=""))
    plot_cnv(plot_data=data_smoothed,
                       contigs=paste(as.vector(as.matrix(chr_order))),
                       k_obs_groups=num_obs_groups,
                       reference_idx=reference_obs,
                       ref_groups=groups_ref,
                       pdf_path=pdf_path,
                       color_safe_pal=color_safe)
}

#' Log intermediate step with a plot and text file of the steps.
#'
#' Args:
#'    @params data: The data frame to plot.
#'    @params plot_name: The absolute path to the pdf to be plotted.
#'
#' Returns:
#'    @return No return
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

#' Plot the matrix as a heatmap.
#' Clustering is on observation only, gene position is preserved.
#'
#' Args:
#'    @param plot_data: Data matrix to plot (columns are observations).
#'    @param contigs: The contigs the data is group in in order of rows.
#'    @param reference_idx: Vector of reference indices.
#'    @param reg_groups: Groups of vector indices (as indices in reference_idx)
#'    @param pdf_path: Path to save pdf file.
#'    @param color_safe_pal: Logical indication of using a color blindness safe
#'                          palette.
#'
#' Returns:
#'    @return No return
plot_cnv <- function(plot_data,
                     contigs,
                     reference_idx,
                     ref_groups,
                     pdf_path,
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
    for (contig_name in names(contig_tbl)){
        contig_labels <- c(contig_labels,
                           contig_name,
                           rep("", contig_tbl[contig_name] - 1))
    }

    # Heatmap elements widths
    heatmap_widths <- c(1,6)

    # Rows observations, Columns CHR
    # Plot either high or low resolution
    pdf(pdf_path,
        useDingbats=FALSE,
        width=10,
        height=7.5,
        paper="USr")

    # Plot heatmap
    par(mar=c(.5,.5,.5,.5))

    # Plot observations
    ## Make Observation Samples
    ## Remove observation col names, too many to plot
    ## Will try and keep the reference names
    ## They are more informative anyway
    obs_data <- plot_data
    if (!is.null(ref_idx)){
        obs_data <- plot_data[, -1 * ref_idx, drop=FALSE]
        #TODO there will always be more obs, flip this with ref.
        if (ncol(obs_data) == 1){
                plot_data <- cbind(obs_data, obs_data)
                names(obs_data) <- c("", names(obs_data)[1])
        }
    }

    ## Plot observational samples
    plot_cnv_observations(obs_data=t(obs_data),
                          file_base_name=pdf_path,
                          heatmap_widths=heatmap_widths,
                          contig_colors=ct.colors[contigs],
                          contig_label=contig_labels,
                          col_pal=custom_pal,
                          contig_seps=col_sep,
                          num_obs_groups=k_obs_groups)
    obs_data <- NULL


    # Plot Reference Samples
    #if(!is.null(ref_idx)){
    #    plot_cnv_references(ref_data=plot_data[, ref_idx, drop=FALSE],
    #                        ref_groups=ref_groups,
    #                        col_pal=custom_pal,
    #                        contig_seps=col_sep,
    #                        file_base_name=pdf_path,
    #                        heatmap_widths=heatmap_widths)
    #}
    dev.off()
}

#' Plot the observational samples
#'
#' Args:
#'    @param obs_data: Data to plot as observations. Rows = Cells, Col = Genes
#'    @param col_pal: The color palette to use.
#'    @param contig_colors: The colors for the contig bar.
#'    @param contig_labels: The labels for the contigs.
#'    @param contig_seps: Indices for line seperators of contigs.
#'    @param file_base_name: Base of the file to used to make output file names.
#'    @param heatmap_widths: Width of heatmap column
#'
#' Returns:
#'    @return Void
plot_cnv_observations <- function(obs_data,
                                  col_pal,
                                  contig_colors,
                                  contig_labels,
                                  contig_seps,
                                  num_obs_groups,
                                  file_base_name,
                                  heatmap_widths){

    logging::loginfo("plot_cnv_observation:Start")
    logging::loginfo(paste("Observation data size: Cells=",
                           nrow(obs_data),
                           "Genes=",
                           ncol(obs_data),
                           sep=" "))
    observation_file_base <- paste(file_base_name, "_observations.txt", sep="")

    # Output dendrogram representation as Newick
    # Need to precompute the dendrogram so we can manipulate
    # it before the heatmap plot
    obs_hcl <- hclust(sparseEuclid(obs_data),"average")
    obs_dendrogram <- as.dendrogram(obs_hcl)
    write.tree(as.phylo(obs_hcl),
               file=paste(file_base_name, "_observations_dendrogram.txt", sep=""))

    # Output HCL group membership.
    # Record locations of seperations
    obs_seps <- c(0)
    ordered_names <- rev(row.names(obs_data)[obs_hcl$order])
    split_groups <- cutree(obs_hcl, k=num_obs_groups)
    row_groupings <- get_group_color_palette()(length(table(split_groups)))[split_groups]
    for (cut_group in unique(split_groups)){
        group_memb = names(split_groups)[which(split_groups == cut_group)]
        # Write group to file
        memb_file <- file(paste(file_base_name,cut_group,"members.txt",sep="_"))
        writeLines(group_memb, memb_file)
        close(memb_file)
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
    #par(fig=c(0,1,0,1), new=FALSE)
    data_observations <- GMD::heatmap.3(obs_data,
                                        Rowv=obs_dendrogram,
                                        Colv=FALSE,
                                        cluster.by.col=FALSE,
                                        main="Copy Number Variation Inference",
                                        ylab="Observations (Cells)",
                                        margin.for.labRow=10,
                                        margin.for.labCol=2,
                                        xlab="Genomic Region",
                                        key=TRUE,
                                        labCol=contig_labels,
                                        notecol="black",
                                        density.info="histogram",
                                        denscol="blue",
                                        trace="none",
                                        dendrogram="row",
                                        cexRow=0.8,
                                        scale="none",
                                        color.FUN=col_pal,
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
                                        key.ylab="Count")

    # Write data to file.
    logging::loginfo(paste("plot_cnv_references:Writing observation data to",
                           observation_file_base,
                           sep=" "))
    row.names(obs_data) <- orig_row_names
    write.table(obs_data[data_observations$rowInd,data_observations$colInd],
                file=observation_file_base)
}

#' Plot the reference samples
#'
#' Args:
#'    @param ref_data: Data to plot as references. Rows = Cells, Col = Genes
#'    @param ref_groups: Groups of references to plot together.
#'    @param col_pal: The color palette to use.
#'    @param contig_seps: Indices for line seperators of contigs
#'    @param file_base_name: Base of the file to used to make output file names.
#'    @param heatmap_widths: Width of heatmap column
#'
#' Returns:
#'    @return Void
plot_cnv_references <- function(ref_data,
                                ref_groups,
                                col_pal,
                                contig_seps,
                                file_base_name,
                                heatmap_widths){

    logging::loginfo("plot_cnv_references:Start")
    logging::loginfo(paste("Reference data size: Cells=",
                           ncol(ref_data),
                           "Genes=",
                           nrow(ref_data),
                           sep=" "))
    number_references <- ncol(ref_data)
    reference_ylab <- NA
    reference_data_file <- paste(file_base_name, "_references.txt", sep="")

    ref_seps <- c()
    # Handle only one reference
    # heatmap2 requires a 2 x 2 matrix, so with one reference
    # I just duplicate the row and hid the second name so it
    # visually looks like it is just taking up the full realestate.
    if(number_references == 1){
        ref_data <- cbind(ref_data, ref_data)
        names(ref_data) <- c("",names(ref_data)[1])
    }

    # Handle reference groups
    # If there is more than one reference group, visually break
    # up the groups with a row seperator. Also plot the rows in
    # order so the current groups are show and seperated.
    if(length(ref_groups) > 1){
        i_cur_idx <- 0
        order_idx <- c()
        for (ref_grp in ref_groups){
            i_cur_idx <- i_cur_idx + length(ref_grp)
            ref_seps <- c(ref_seps, i_cur_idx)
            order_idx <- c(order_idx, ref_grp)
        }
        ref_seps <- ref_seps[1:(length(ref_seps) - 1)]
        ref_data <- ref_data[, order_idx, drop=FALSE]
    }
    logging::loginfo(paste("plot_cnv_references:Number reference groups=",
                           length(ref_groups)),
                           sep=" ")
    ref_data <- t(ref_data)
    if (number_references > 20){
        # The reference labs can become clutered
        # Dynamically change labels given a certain number of labels.
        reference_ylab <- "References"
        row.names(ref_data) <- rep("", number_references)
    }
    # Print controls
    par(fig=c(0,1,0,1), new=TRUE)
    data_references <-GMD::heatmap.3(ref_data,
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
                   cexRow=0.8,
                   scale="none")
    #               rowsep=ref_seps,
    #               col=col_pal,
    #               # Seperate by contigs
    #               colsep=contig_seps,
    #               # Seperate by reference / not reference
    #               sepcolor="black",
    #               sepwidth=c(0.01,0.01),
    #               # Color by contigs
    #               lmat=rbind(c(2,1),c(4,5),c(3,0)),
    #               lhei=c(.25,.8,1.6),
    #               lwid=heatmap_widths)

    #par(fig=c(0,1,0,1), new=TRUE)
    #data_references <- gplots::heatmap.2(ref_data,
    #                             main=NA,
    #                             ylab=reference_ylab,
    #                             xlab=NA,
    #                             key=FALSE,
    #                             labCol=rep("", nrow(ref_data)),
    #                             notecol="black",
    #                             trace="none",
    #                             dendrogram="none",
    #                             Colv=FALSE,
    #                             Rowv=FALSE,
    #                             cexRow=0.8,
    #                             scale="none",
    #                             rowsep=ref_seps,
    #                             col=col_pal,
    #                             # Seperate by contigs
    #                             colsep=contig_seps,
    #                             # Seperate by reference / not reference
    #                             sepcolor="black",
    #                             sepwidth=c(0.01,0.01),
    #                             # Color by contigs
    #                             lmat=rbind(c(2,1),c(4,5),c(3,0)),
    #                             lhei=c(.25,.8,1.6),
    #                             lwid=heatmap_widths)

    # Write data to file
    logging::loginfo(paste("plot_cnv_references:Writing reference data to",
                           reference_data_file,
                           sep=" "))
    write.table(ref_data[data_references$rowInd,data_references$colInd],
                file=reference_data_file)
}

#' Return the indices of the rows that average above the cut off
#'
#' Args:
#'    @param data: Data to measure the average row and evaluate
#'                 against the cutoff. Row = Genes, Col = Cells.
#'    @param cutoff: Threshold to be above to be kept.
#'
#' Returns:
#'    @return Returns a vector of row indicies to keep (are above the cutoff).
above_cutoff <- function(data, cutoff){

    logging::loginfo(paste("::above_cutoff:Start", sep=""))
    average_gene <- log2(rowMeans( ( ( (2 ^ data) - 1) * 10 ), na.rm=TRUE) + 1 )
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
#'    @param data: Data (expression) matrix where the row names should be in
#'                 the row names of the genomic_position file.
#'    @param genomic_position: Data frame read in from the genomic position file
#'
#' Returns:
#'    @return Returns a matrix of expression in the order of the
#'            genomic_position file. NULL is returned if the genes in both
#'            data parameters do not match.
order_reduce <- function(data, genomic_position){

    logging::loginfo(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL)
    if (is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }
    # Reduce to genes in pos file
    keep_genes <- row.names(data)[which(row.names(data)
                                  %in% row.names(genomic_position))]

    # Set the chr to factor so the order can be arbitrarily set and sorted.
    chr_levels <- unique(genomic_position[[CHR]])
    genomic_position[[CHR]] <- factor(genomic_position[[CHR]],
                                   levels=chr_levels)
    if(length(keep_genes)){
        ret_results$expr <- data[keep_genes, , drop=FALSE]
        ret_results$order <- genomic_position[keep_genes, , drop=FALSE]
    }
    return(ret_results)
}

#' Remove values that are too close to the average and are considered noise.
#'
#' Args:
#'    @param smooth_matrix: A matrix of values, smoothed, and with average 
#'                          reference removed. Row = Genes, Col = Cells.
#'    @param threshold: The amount of difference a value must be from the
#'                      reference before the value can be kept and not 
#'                      removed as noise.
#' Returns:
#'    @return Denoised matrix
remove_noise <- function(smooth_matrix, threshold){

    logging::loginfo(paste("::remove_noise:Start.", sep=""))
    if (threshold > 0){
        smooth_matrix[abs(smooth_matrix) < threshold] <- 0
    }
    return(smooth_matrix)
}


#' Remove the tails of values of a specific chromosome.
#' The smooth_matrix values are expected to be in genomic order.
#' If the tail is too large and no contig will be left 1/3 of the
#' contig is left.
#'
#' Args:
#'    @param smooth_matrix: Smoothed values in genomic order.
#'                          Row = Genes, Col = Cells.
#'    @param chr: Indices of the chr in which the tails are to be removed.
#'    @param tail_length: Length of the tail to remove on both ends of the
#'                        chr indices.
#' Returns:
#'    @return Indices to remove.
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

#' Smooth a matrix by column using a simple moving average.
#' Tails of the averages use a window length that is truncated to
#' available data.
#'
#' Args:
#'    @param data: Data matrix to smooth. Row = Genes, Col = Cells.
#'    @param window_length: Length of window to use for the moving average.
#'        Should be a positive, odd integer.
#'
#' Returns:
#'    @return Matrix with columns smoothed with a simple moving average.
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

#' Helper function for smoothing the ends of a moving average.
#'
#' Args:
#'    @param obs_data: Data to smooth
#'    @param obs_tails: Length of the tail to smooth.
#' Returns:
#'    @return: Data smoothed.
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

#' Smooth vector of values over the given window length.
#'
#' Args:
#'    @param obs_data: Vector of data to smooth with a moving average.
#'    @param window_length: Length of the window for smoothing.
#'        Must be and odd, positive, integer.
#' Returns:
#'    @return Vector of values smoothed with a moving average.
smooth_window_helper <- function(obs_data, window_length){

    return(filter(obs_data, rep(1 / window_length, window_length), sides=2))
}


#' Measure an euclidean distance weighted by complexity.
#' It is assumed rows are clustered
#' Returns
#'    @return Distance object of Euclidean distance weighted by complexity.
sparseEuclid <- function(data){
    feature_count <- nrow(data)
    obs_complex <- apply(data, MARGIN=2, FUN=function(x) return(sum(x!=0)))
    per_complex <- obs_complex / feature_count
    out_matrix <- matrix(rep(NA,feature_count^2),ncol=feature_count)
    names(out_matrix) <- row.names(data)
    row.names(out_matrix) <- row.names(data)
    for( fi_one in 1:feature_count ){
        for( fi_two in fi_one:feature_count ){
            f1 <- data[fi_one,]
            f2 <- data[fi_two,]
            d <- sqrt(sum(per_complex * (f1 - f2)^2))
            out_matrix[fi_one, fi_two] <- d
            out_matrix[fi_two, fi_one] <- d
        }
    }
    return(as.dist(out_matrix))
}

# Load libraries
library(ape)
library("RColorBrewer", character.only=TRUE)
library(GMD)
library(gplots)
library(optparse)
library(logging)

# If called as source from commandline
if (identical(environment(),globalenv()) &&
    !length(grep("^source\\(", sys.calls()))){

    # Logging level choices
    C_LEVEL_CHOICES <- names(loglevels)
    # Visualization outlier thresholding and bounding method choices
    C_VIS_OUTLIER_CHOICES <- c("average_bound")

    # Command line arguments
    pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                      "--pdf pdf_file",
                                      "data_matrix genomic_positions"))

    pargs <- optparse::add_option(pargs, c("--color_safe"),
                        type="logical",
                        default=FALSE,
                        action="store_true",
                        dest="use_color_safe",
                        metavar="Color_Safe",
                        help=paste("To support the needs of those who see ",
                                   "colors differently, use this option to",
                                   "change the colors to a palette visibly ",
                                   "distinct to all color blindness. ",
                                   " [Default %default]"))

    pargs <- optparse::add_option(pargs, c("--cutoff"),
                        type="integer",
                        default=0,
                        action="store",
                        dest="cutoff",
                        metavar="Cutoff",
                        help=paste("A number >= 0 is expected. A cut off for",
                                   "the average expression of genes to be used",
                                   "for CNV inference. [Default %default]"))

    pargs <- optparse::add_option(pargs, c("--transform"),
                        type="logical",
                        default=FALSE,
                        action="store_true",
                        dest="log_transform",
                        metavar="LogTransform",
                        help=paste("Matrix is assumed to be Log2(TPM+1) ",
                                   "transformed. If instead it is raw TPMs ",
                                   "use this flag so that the data will be ",
                                   "transformed. [Default %default]"))

    pargs <- optparse::add_option(pargs, c("--log"),
                        type="character",
                        action="store",
                        default=NA,
                        dest="log_file",
                        metavar="Log",
                        help=paste("File for logging. If not given,",
                                   "logging will occur to console.",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--log_level"),
                        type="character",
                        action="store",
                        default="INFO",
                        dest="log_level",
                        metavar="LogLevel",
                        help=paste("Logging level. Valid choices are",
                                   paste(C_LEVEL_CHOICES,collapse=", "),
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--noise_filter"),
                        type="numeric",
                        default=0,
                        action="store",
                        dest="magnitude_filter",
                        metavar="Magnitude_Filter",
                        help=paste("A value must be atleast this much more or",
                                   "less than the reference to be plotted",
                                   "[Default %default]."))

    pargs <- optparse::add_option(pargs, c("--max_centered_expression"),
                        type="integer",
                        default=3,
                        action="store",
                        dest="max_centered_expression",
                        metavar="Max_centered_expression",
                        help=paste("This value and -1 * this value are used",
                                   "as the maximum value expression that can",
                                   "exist after centering data. If a value is",
                                   "outside of this range, it is truncated to",
                                   "be within this range [Default %default]."))

    pargs <- optparse::add_option(pargs, c("--obs_groups"),
                        type="character",
                        default=1,
                        action="store",
                        dest="num_obs_groups",
                        metavar="Number_of_observation_groups",
                        help=paste("Number of groups in which to break ",
                                   "the observations.",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--pdf"),
                        type="character",
                        action="store",
                        dest="pdf_file",
                        metavar="Output_PDF_Visualization",
                        help=paste("Output PDF for visualizing the analysis.",
                                   "[Default %default][REQUIRED]"))

    pargs <- optparse::add_option(pargs, c("--ref"),
                        type="character",
                        default=NULL,
                        action="store",
                        dest="reference_observations",
                        metavar="Input_reference_observations",
                        help=paste("Tab delimited characters are expected.",
                                   "Names of the subset of samples ( data's",
                                   "columns ) that should be used as",
                                   "references if not given, the average of",
                                   "all samples will be the reference.",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--ref_groups"),
                        type="character",
                        default=1,
                        action="store",
                        dest="num_groups",
                        metavar="Number_of_reference_groups",
                        help=paste("Number of groups in which to break ",
                                   "the reference observations.",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--steps"),
                        type="logical",
                        default=FALSE,
                        action="store_true",
                        dest="plot_steps",
                        metavar="plot_steps",
                        help=paste("Using this argument turns on plotting ",
                                   "intemediate steps. The plots will occur ",
                                   "in the same directory as the output pdf. ",
                                   "Please note this option increases the time",
                                   " needed to run [Default %default]"))

    pargs <- optparse::add_option(pargs, c("--vis_bound_method"),
                        type="character",
                        default=NA,
                        action="store",
                        dest="bound_method_vis",
                        metavar="Outlier_Removal_Method_Vis",
                        help=paste("Method to automatically detect and bound",
                                   "outliers. Used for visualizing. If both",
                                   "this argument and ",
                                   "--vis_bound_threshold are given, this will",
                                   "not be used. Valid choices are",
                                   paste(C_VIS_OUTLIER_CHOICES, collapse=", "),
                                   " [Default %default]"))

    pargs <- optparse::add_option(pargs, c("--vis_bound_threshold"),
                        type="character",
                        default=NA,
                        action="store",
                        dest="bound_threshold_vis",
                        metavar="Outlier_Removal_Threshold_Vis",
                        help=paste("Used as upper and lower bounds for values",
                                   "in the visualization. If a value is",
                                   "outside this bound it will be replaced by",
                                   "the closest bound. Should be given in",
                                   "the form of 1,1 (upper bound, lower bound)",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--window"),
                        type="integer",
                        default=101,
                        action="store",
                        dest="window_length",
                        metavar="Window_Lengh",
                        help=paste("Window length for the smoothing.",
                                   "[Default %default]"))

    pargs <- optparse::add_option(pargs, c("--tail"),
                        type="integer",
                        default=NA,
                        action="store",
                        dest="contig_tail",
                        metavar="contig_tail",
                        help=paste("Contig tail to be removed.",
                                   "[Default %default]"))

    args_parsed <- optparse::parse_args(pargs, positional_arguments=2)
    args <- args_parsed$options
    args["input_matrix"] <- args_parsed$args[1]
    args["gene_order"] <- args_parsed$args[2]

    # Check arguments
    args <- check_arguments(args)

    # Parse bounds
    bounds_viz <- c(NA,NA)
    if (!is.na(args$bound_threshold_vis)){
        bounds_viz <- as.numeric(unlist(strsplit(args$bound_threshold_vis,",")))
    }
    if (length(bounds_viz) != 2){
        error_message <- paste("Please use the correct format for the argument",
                               "--vis_bound_threshold . Two numbers seperated",
                               "by a comma is expected (lowerbound,upperbound)",
                               ". As an example, to indicate that outliers are",
                               "outside of -1 and 1 give the following.",
                               "--vis_bound_threshold -1,1")
        stop(error_message)
    }

    # Set up logging file
    logging::basicConfig(level=args$log_level)
    if (!is.na(args$log_file)){
        logging::addHandler(logging::writeToFile,
                             file=args$log_file,
                             level=args$log_level)
    }

    # Log the input parameters
    logging::loginfo(paste("::Input arguments. Start.")) 
    for (arg_name in names(args)){
        logging::loginfo(paste(":Iput_Argument:",arg_name,"=",args[[arg_name]],
                               sep="")) 
    }
    logging::loginfo(paste("::Input arguments. End.")) 

    # Manage inputs
    logging::loginfo(paste("::Reading data matrix.", sep=""))
    # Row = Genes/Features, Col = Cells/Observations
    expression_data <- read.table(args$input_matrix)
    logging::loginfo(paste("Original matrix dimensions (r,c)=",
                  paste(dim(expression_data), collapse=",")))

    # Default the gene order to the given order
    # If given a file read the user defined order.
    input_gene_order <- seq(1, nrow(expression_data), 1)
    if (args$gene_order != ""){
        input_gene_order <- read.table(args$gene_order, row.names=1)
        names(input_gene_order) <- c(CHR, START, STOP)
    }
    logging::loginfo(paste("::Reading gene order.", sep=""))
    logging::logdebug(paste(head(args$gene_order[1]), collapse=","))

    # Default the reference samples to all
    input_reference_samples <- colnames(expression_data)
    if (!is.null(args$reference_observations)){
        # This argument can be either a list of column labels
        # which is a comma delimited list of column labels
        # holding a comma delimited list of column labels
        refs <- args$reference_observations
        if (file.exists(args$reference_observations)){
            refs <- scan(args$reference_observations,
                         what="character",
                         quiet=TRUE)
            refs <- paste(refs, collapse=",")
        }
        # Split on comma
        refs <- unique(unlist(strsplit(refs, ",", fixed=FALSE)))
        # Remove multiple spaces to single spaces
        refs <- unique(unlist(strsplit(refs, " ", fixed=FALSE)))
        refs <- refs[refs != ""]
        if (length(refs) > 0){
            input_reference_samples <- refs
        }
    }

    # Make sure the given reference samples are in the matrix.
    if (length(input_reference_samples) !=
        length(intersect(input_reference_samples, colnames(expression_data)))){
        missing_reference_sample <- setdiff(input_reference_samples,
                                            colnames(expression_data))
        error_message <- paste("Please make sure that all the reference sample",
                              "names match a sample in your data matrix.",
                              "Attention to: ",
                              paste(missing_reference_sample, collapse=","))
        logging::logerror(error_message)
        stop(error_message)
    }

    # Order and reduce the expression to the genomic file.
    order_ret <- order_reduce(data=expression_data,
                                        genomic_position=input_gene_order)
    expression_data <- order_ret$expr
    input_gene_order <- order_ret$order
    if(is.null(expression_data)){
        error_message <- paste("None of the genes in the expression data",
                              "matched the genes in the reference genomic",
                              "position file. Analysis Stopped.")
        stop(error_message)
    }

    # Run CNV inference
    infer_cnv(data=expression_data,
                        gene_order=input_gene_order,
                        cutoff=args$cutoff,
                        reference_obs=input_reference_samples,
                        transform_data=args$log_transform,
                        window_length=args$window_length,
                        max_centered_threshold=args$max_centered_expression,
                        noise_threshold=args$magnitude_filter,
                        num_ref_groups=args$num_groups,
                        num_obs_groups=args$num_obs,
                        pdf_path=args$pdf_file,
                        plot_steps=args$plot_steps,
                        contig_tail=args$contig_tail,
                        color_safe=args$use_color_safe,
                        method_bound_vis=args$bound_method_vis,
                        lower_bound_vis=bounds_viz[1],
                        upper_bound_vis=bounds_viz[2])
}
