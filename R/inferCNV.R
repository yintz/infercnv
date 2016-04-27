#!/usr/bin/env Rscript


#' Remove the average of the genes of the reference observations from all 
#' observations' expression.
#'
#' Args:
#'    @param average_data: Matrix containing the data to remove average from
#'                         (this includes the reference observations).
#'                         Row = Genes, Col = Cells.
#'    @param reference_observations: Indices of reference observations.
#'                                   Only these are used in the average.
#'
#' Returns:
#'    @return Expression with the average gene expression in the reference 
#'            observations removed.
average_over_ref <- function(average_data,
                             ref_observations){

    logging::loginfo(paste("::average_over_ref:Start", sep=""))
    average_reference_obs <- average_data[, ref_observations]
    if (length(ref_observations) > 1){
        average_reference_obs <- rowMeans(average_data[, ref_observations],
                                          na.rm=TRUE)
    }
    return(list(data=sweep(average_data, MARGIN=1,
                           average_reference_obs, FUN="-"),
                averages=average_reference_obs))
}

get_adj_group_average <- function(data_adj,
                                  averages,
                                  ref_observations,
                                  ref_groups){

    logging::loginfo(paste("::get_adj_group_average:Start", sep=""))
    average_reference_obs <- data_adj[, ref_observations]
    ret_adj = list()
    adj_max <- rep(0,length(averages))
    adj_min <- rep(0,length(averages))
    for (ref_group in ref_groups){
        grp_average <- rowMeans(average_reference_obs[,ref_group], na.rm=TRUE)
        grp_diff <- grp_average - averages
        adj_max <- pmax(adj_max, grp_diff)
        adj_min <- pmin(adj_min, grp_diff)
    }
    return(list(max_adj=adj_max, min_adj=adj_min))
}

#' Split up reference observations in to k groups and return indices
#' for the different groups.
#'
#' Args:
#'    @param average_data: Matrix containing data. Row = Genes, Col = Cells.
#'    @param ref_obs: Indices of reference obervations.
#'    @param num_groups: The number of groups to partition nodes in.
#'
#' Returns:
#'    @return Returns a list of grouped reference observations given as
#'            vectors of groups.
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

    # Get grouping
    if (length(num_groups) == 1){
        num_groups <- unlist(num_groups)
        if (num_groups > ncol(average_data)){
            num_groups <- ncol(average_data)
        }

        # If only one group is needed, short-circuit and
        # return all indices immediately
        if((num_groups < 2)||
           (length(ref_obs) < 2)){
            ret_groups[[1]] <- ref_obs
            ret_order <- ref_obs
            return(ret_groups)
        }

        # Get HCLUST
        # Get reference observations only.
        average_reference_obs <- t(average_data[, ref_obs])
        hc <- hclust(dist(average_reference_obs))

        split_groups <- cutree(hc, k=num_groups)
        split_groups <- split_groups[hc$order]
        # Keep the sort of the hclust
        for(cut_group in unique(split_groups)){
            ret_groups[[cut_group]] <- hc$order[which(split_groups==cut_group)]
        }
    } else {
        ret_groups <- num_groups
    }
    return(ret_groups)
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
    # Warn that an average of the samples is used in the absence of
    # normal / reference samples
    if (is.null(arguments$reference_observations)){
        logging::logwarn(paste(":: --reference_observations: No reference ",
                      "samples were given, the average of the samples ",
                      "will be used.",
                      sep=""))
    }

    # Require the contig tail to be above 0
    if (is.na(arguments$contig_tail)){
        arguments$contig_tail <- (arguments$window_length - 1) / 2
    }
    if (arguments$contig_tail < 0){
        logging::logerror(paste(":: --tail: Please enter a value",
                                "greater or equal to zero for the tail.",
                                sep=""))
    }

    if (! is.na(as.integer(arguments$num_groups))){
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
            exit(10)
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
                exit(9)
            }
        }
        arguments$num_groups <- num_groups
    }
    return(arguments)
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
#'    @param reference_obs: Indices of the subset of samples (data's columns)
#'                          that should be used as references.
#'                          If not given, the average of all samples will 
#'                          be the reference.
#'    @param window_length: Length of the window for the moving average
#'                          (smoothing). Should be an odd integer.
#'    @param max_centered_threshold: The maximum value a a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value.
#'    @param noise_threshold: The minimum difference a value can be from the 
#'                            average reference in order for it not to be
#'                            removed as noise.
#'    @param pdf_path: The path to what to save the pdf as. The raw data is
#'                     also written to this path but with the extension .txt .
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
                      pdf_path,
                      contig_tail=(window_length - 1) / 2){

    logging::loginfo(paste("::infer_cnv:Start", sep=""))
    # Remove any gene without position information
    remove_by_position <- -1 * which(gene_order[2] + gene_order[3] == 0)
    if(length(remove_by_position)){
        gene_order <- gene_order[remove_by_position,]
        data <- data[remove_by_position,]
    }
    logging::loginfo(paste("::infer_cnv:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    logging::logdebug(paste("::infer_cnv:Removed indices:"))
    #logging::logdebug(paste(remove_by_position * -1, collapse=","))
    #write.table(data, file="01_remove_pos.txt")

    # Make sure data is log transformed + 1
    if (transform_data){
        data <- log2(data + 1)
    }

    # Reduce by cutoff
    keep_gene_indices <- above_cutoff(data, cutoff)
    if (!is.null(keep_gene_indices)){
        data <- data[keep_gene_indices,]
        gene_order <- gene_order[keep_gene_indices,]
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff, ",
                      "new dimensions (r,c) = ",
                      paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
        logging::logdebug(paste("::infer_cnv:Keeping indices.", sep=""))
        #logging::logdebug(paste(keep_gene_indices, collapse=","))
    } else {
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff.", sep=""))
        logging::logwarn(paste("::infer_cnv::No indicies left to keep.",
                      " Stoping."))
        stop()
    }
    #write.table(data, file="02_cutoff.txt")
    #write.table(gene_order, file="02b_cutoff.txt")

    # Order data by genomic region
    data <- data[with(gene_order, order(chr,start,stop)),]
    # Order gene order the same way
    chr_order <- gene_order[with(gene_order, order(chr,start,stop)),][1]
    gene_order <- NULL
    #write.table(data, file="03a_ordered.txt")
    #write.table(chr_order, file="03b_ordered_chr.txt")

    # Center data (automatically ignores zeros)
    data <- center_with_threshold(data, max_centered_threshold)
    logging::loginfo(paste("::infer_cnv:Outlier removal, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    #write.table(data, file="04_center.txt")

    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_length)
    data <- NULL
    logging::loginfo(paste("::infer_cnv:Smoothed data.", sep=""))

    # Split the reference data into groups if requested
    groups_ref <- split_references(average_data=data_smoothed,
                                   ref_obs=reference_obs,
                                   num_groups=num_ref_groups)
    logging::loginfo(paste("::infer_cnv:split_reference. ",
                           "found ",length(groups_ref)," reference groups.",
                          sep=""))

    # Remove average reference
    ret_average <- average_over_ref(average_data=data_smoothed,
                                      ref_observations=reference_obs)
    data_smoothed <- ret_average$data
    logging::loginfo(paste("::infer_cnv:Remove average, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    #write.table(data_smoothed, file="05_remove_average.txt")

    # Get adjustments for removing noise given multiple reference groups
    # Which may have different patterns
    average_adjustments <- get_adj_group_average(data_adj=data_smoothed,
                                                 averages=ret_average$averages,
                                                 ref_observations=reference_obs,
                                                 ref_groups=groups_ref)
    logging::loginfo(paste("::infer_cnv:get_adj_group_average. End"))

    # Remove Ends
    logging::logdebug(chr_order)
    for (chr in unlist(unique(chr_order))){
        logging::loginfo(paste("::infer_cnv:Remove tail contig ",
                               chr,
                               ".",
                               sep=""))
        data_smoothed <- remove_tails(data_smoothed,
                                      which(chr_order == chr),
                                      contig_tail)
    }
    logging::loginfo(paste("::infer_cnv:Remove ends, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    #write.table(data_smoothed, file="05_remove_ends.txt")

    # Remove noise
    data_smoothed <- remove_noise(ref=reference_obs,
                                  smooth_matrix=data_smoothed,
                                  threshold=noise_threshold,
                                  max_adj=average_adjustments$max_adj,
                                  min_adj=average_adjustments$min_adj)
    logging::loginfo(paste("::infer_cnv:Remove moise, ",
                           "new dimensions (r,c) = ",
                           paste(dim(data_smoothed), collapse=","),
                           " Total=", sum(data_smoothed),
                           " Min=", min(data_smoothed),
                           " Max=", max(data_smoothed),
                           ".", sep=""))
    #write.table(data_smoothed, file="06_remove_noise.txt")

    # Plot and write data
    logging::loginfo(paste("::infer_cnv:Drawing plots to file:",
                  pdf_path, sep=""))
    logging::loginfo(paste("::infer_cnv:Current data dimensions (r,c)=",
                           paste(dim(data_smoothed), collapse=","), sep=""))
    plot_cnv(plot_data=data_smoothed,
             contigs=paste(as.vector(as.matrix(chr_order))),
             reference_idx=reference_obs,
             ref_groups=groups_ref,
             pdf_path=pdf_path)
    logging::loginfo(paste("::infer_cnv:Writing final data to ",
                  paste(pdf_path, ".txt", sep=""), sep=""))
    write.table(data_smoothed, file=paste(pdf_path, ".txt", sep=""))
}

#' Plot the matrix as a heatmap.
#' Clustering is on observation only, gene position is preserved.
#'
#' Args:
#'    @param plot_data: Data matrix to plot (columns are observations).
#'    @param pdf_path: Path to save pdf file.
#'
#' Returns:
#'    @return No return
plot_cnv <- function(plot_data,
                     contigs,
                     reference_idx,
                     ref_groups,
                     pdf_path){

    logging::loginfo(paste("::plot_cnv:Start", sep=""))
    logging::logdebug(paste("::plot_cnv:Current data dimensions (r,c)=",
                            paste(dim(plot_data), collapse=","),
                            " Total=", sum(plot_data),
                            " Min=", min(plot_data),
                            " Max=", max(plot_data),
                            ".", sep=""))

    # Contigs
    unique_contigs <- unique(contigs)
    ct.colors <- colorRampPalette(brewer.pal(12,"Set3"))(length(unique_contigs))
    names(ct.colors) <- unique_contigs

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
    pdf(pdf_path,
        useDingbats=FALSE,
        width=10,
        height=7.5,
        paper="USr")

    # Plot heatmap
    par(mar=c(.5,.5,.5,.5))

    # Plot Observation Samples
    # Remove observation col names, too many to plot
    # Will try and keep the reference names
    # They are more informative anyway
    obs_data <- plot_data
    if (!is.null(ref_idx)){
        obs_data <- plot_data[,-1 * ref_idx,drop=FALSE]
        #TODO there will always be more obs, flip this with ref.
        if (ncol(obs_data) == 1){
                plot_data <- cbind(obs_data, obs_data)
                names(obs_data) <- c("", names(obs_data)[1])
        }
    }
    obs_data <- t(obs_data)
    row.names(obs_data) <- rep("", nrow(obs_data))

    par(fig=c(0,1,0,1), new=FALSE)
    heatmap.2(obs_data,
        main="Copy Number Variation Inference",
        ylab="Observations (Cells)",
        xlab="Genomic Region",
        key=TRUE,
        labCol=contig_labels,
        notecol="black",
        density.info="histogram",
        denscol="blue",
        trace="none",
        dendrogram="row",
        Colv=FALSE,
        cexRow=0.8,
        scale="none",
        #col="cm.colors",
        col=colorRampPalette(c("purple","white","orange"))(n=100),
        # Seperate by contigs
        colsep=col_sep,
        # Seperate by reference / not reference
        sepcolor="black",
        sepwidth=c(0.01,0.01),
        # Color by contigs
        ColSideColors=ct.colors[contigs],
        # Position heatmap elements
        lmat=rbind(c(5,4),c(0,1),c(3,2)),
        lhei=c(1.9,0.1,4),
        lwid=heatmap_widths)
    obs_data <- NULL

    # Plot Reference Samples
    if(!is.null(ref_idx)){
        plot_data <- plot_data[,ref_idx, drop=FALSE]
        number_references <- ncol(plot_data)
        reference_ylab <- NA
        ref_seps = c()
        # heatmap2 requires a 2 x 2 matrix, so with one reference
        # I just duplicate the row and hid the second name so it
        # visually looks like it is just taking up the full realestate.
        if(number_references == 1){
            plot_data <- cbind(plot_data, plot_data)
            names(plot_data) <- c("",names(plot_data)[1])
        }
        # If there is more than one reference group, visually break
        # up the groups with a row seperator. Also plot the rows in
        # order so the current groups are show and seperated.
        if(length(ref_groups) > 1){
            i_cur_idx <- 0
            order_idx <- c()
            for(ref_grp in ref_groups){
                i_cur_idx <- i_cur_idx + length(ref_grp)
                ref_seps <- c(ref_seps, i_cur_idx)
                order_idx <- c(order_idx, ref_grp)
            }
            ref_seps = ref_seps[1:(length(ref_seps) - 1)]
            plot_data <- plot_data[,order_idx, drop=FALSE]
        }
        plot_data <- t(plot_data)
        if (number_references > 20){
            # The reference labs can become clutered
            # Dynamically change labels given a certain number of labels.
            reference_ylab <- "References"
            row.names(plot_data) <- rep("", number_references)
        }
        # Print controls
        par(fig=c(0,1,0,1), new=TRUE)
        heatmap.2(plot_data,
            main=NA,
            ylab=reference_ylab,
            xlab=NA,
            key=FALSE,
            labCol=rep("", nrow(plot_data)),
            notecol="black",
            trace="none",
            dendrogram="none",
            Colv=FALSE,
            Rowv=FALSE,
            cexRow=0.8,
            scale="none",
            rowsep=ref_seps,
            #col="cm.colors",
            col=colorRampPalette(c("purple","white","orange"))(n=100),
            # Seperate by contigs
            colsep=col_sep,
            # Seperate by reference / not reference
            sepcolor="black",
            sepwidth=c(0.01,0.01),
            # Color by contigs
            lmat=rbind(c(2,1),c(4,5),c(3,0)),
            lhei=c(.25,.8,1.6),
            lwid=heatmap_widths)
    }
    dev.off()
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
    average_gene <- rowMeans( ( (2 ^ data) - 1), na.rm=TRUE)
    logging::loginfo(paste("::infer_cnv:Averages (counts).", sep=""))
    #logging::logdebug(paste(average_gene, collapse=","))
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
#'    @param genomic_position:
#'
#' Returns:
#'    @return Returns a matrix of expression in the order of the
#'            genomic_position file. NULL is returned if the genes in both
#'            data parameters do not match.
order_reduce <- function(data, genomic_position){

    logging::loginfo(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL)
    if(is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }
    keep_genes <- row.names(data)[which(row.names(data)
                                  %in% row.names(genomic_position))]
    if(length(keep_genes)){
        ret_results$expr <- data[keep_genes,]
        ret_results$order <- genomic_position[keep_genes,]
    }
    return(ret_results)
}

#' Remove values that are too close to the average and are considered noise.
#'
#' Args:
#'    @param ref: Indices of the samples / observation that are reference.
#'        These will not be removed or changed.
#'    @param smooth_matrix: A matrix of values, smoothed, and with average 
#'                          reference removed. Row = Genes, Col = Cells.
#'    @param threshold: The amount of difference a value must be from the
#'                      reference before the value can be kept and not 
#'                      removed as noise.
#' Returns:
#'    @return Denoised matrix
remove_noise <- function(ref, smooth_matrix, threshold,
                         max_adj, min_adj){

    logging::loginfo(paste("::remove_noise:Start.", sep=""))
    if (is.null(max_adj)){
        max_adj <- rep(0, nrows(smooth_matrix))
    }
    if (is.null(min_adj)){
        min_adj <- rep(0, nrows(smooth_matrix))
    }
    if (length(ref) == nrow(smooth_matrix)){
        ref <- c()
    }
    denoised_matrix <- smooth_matrix
    # Adjust the thresholds and turn into vectors
    max_adj <- max_adj + threshold
    max_adj <- matrix(max_adj, nrow=nrow(smooth_matrix),
                      ncol=ncol(smooth_matrix), byrow=FALSE)
    min_adj <- min_adj - threshold
    min_adj <- matrix(min_adj, nrow=nrow(smooth_matrix),
                      ncol=ncol(smooth_matrix), byrow=FALSE)
    # Adjust for max value
    denoised_matrix[max_adj > denoised_matrix & denoised_matrix > 0] <- 0
    denoised_matrix[min_adj < denoised_matrix & denoised_matrix < 0] <- 0
    #denoised_matrix[which(abs(denoised_matrix) < threshold)] <- 0
    if (length(ref) > 0){
        denoised_matrix[, ref] <- smooth_matrix[, ref]
    }
    return(denoised_matrix)
}

#' Remove the tails of values of a specific chromosome.
#' The smooth_matrix values are expected to be in genomic order.
#'
#' Args:
#'    @param smooth_matrix: Smoothed values in genomic order.
#'                          Row = Genes, Col = Cells.
#'    @param chr: Indices of the chr in which the tails are to be removed.
#'    @param tail_length: Length of the tail to remove on both ends of the
#'                        chr indices.
#' Returns:
#'    @return Matrix with tails of chr set to 0.
remove_tails <- function(smooth_matrix, chr, tail_length){

    logging::loginfo(paste("::remove_tails:Start.", sep=""))
    if (tail_length < 1){
        return(smooth_matrix)
    }
    if (length(chr) < (tail_length * 2)){
         smooth_matrix[chr,] <- 0
    } else {
         chr_length <- length(chr)
         smooth_matrix[chr[1:tail_length],] <- 0
         smooth_matrix[chr[( ( chr_length + 1 ) - tail_length ):
                       chr_length],] <- 0
    }
    return(smooth_matrix)
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
        data_sm[row_end,] <- data_end[row_end,]
        data_sm[end_bound, ] <- data_end[end_bound, ]
    }
    return(data_sm)
}

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

    #logging::logdebug(paste("::smooth_window_helper:Start.", sep=""))
    return(filter(obs_data, rep(1 / window_length, window_length), sides=2))
}


# If called as source from commandline
if (identical(environment(),globalenv()) &&
    !length(grep("^source\\(", sys.calls()))){

    library("RColorBrewer", character.only=TRUE)
    library(gplots)
    library(optparse)
    library(logging)

    # Logging level choices
    C_LEVEL_CHOICES <- names(loglevels)

    # Command line arguments
    pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                      "--pdf pdf_file",
                                      "data_matrix genomic_positions"))

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
                        help=paste("Matrix is assumed to be Log2+1 ",
                                   "transformed. If instead it is raw counts ",
                                   "use this flag so that the data will be ",
                                   "transformed. [Default %default]"))

    #pargs <- optparse::add_option(pargs, c("--delim"),
    #                    type="character",
    #                    default=" ",
    #                    action="store",
    #                    dest="delim",
    #                    metavar="Expression_Delimiter",
    #                    help=paste("Delimiter for expression matrix.",
    #                               "[Default %default][REQUIRED]"))

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
                        type="integer",
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

    # Set up logging file
    logging::basicConfig(level=args$log_level)
    if (!is.na(args$log_file)){
        optparse::addHandler(optparse::writeToFile,
                             file=args$log_file,
                             level=args$log_level)
    }

    # Manage inputs
    logging::loginfo(paste("::Reading data matrix.", sep=""))
    # Row = Genes/Features, Col = Cells/Observations
    expression_data <- read.table(args$input_matrix)
    #write.table(expression_data, file="0_data.txt")
    logging::loginfo(paste("Original matrix dimensions (r,c)=",
                  paste(dim(expression_data), collapse=",")))

    # Default the gene order to the given order
    # If given a file read and sort and order.
    input_gene_order <- seq(1, nrow(expression_data), 1)
    if (args$gene_order != ""){
        input_gene_order <- read.table(args$gene_order, row.names=1)
        names(input_gene_order) <- c("chr", "start", "stop")
    }
    logging::loginfo(paste("::Reading gene order.", sep=""))
    logging::logdebug(paste(head(args$gene_order[1]), collapse=","))
    #write.table(input_gene_order, file="0a_data.txt")

    # Default the reference samples to all
    input_reference_samples <- colnames(expression_data)
    if (!is.null(args$reference_observations)){
        input_reference_samples <- unique(
                                       unlist(
                                           strsplit(
                                               args$reference_observations,
                                               ",",
                                               fixed=FALSE)))
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
    #write.table(expression_data, file="0c_order_reduce.txt")
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
                        pdf_path=args$pdf_file,
                        contig_tail=args$contig_tail)
}
