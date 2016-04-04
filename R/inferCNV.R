#!/usr/bin/env Rscript



#' Remove the average of the genes of the reference observations from all 
#' observations' expression.
#'
#' Args:
#'    @param average_data: Matrix containing the data to remove average from
#'                         (this includes the reference observations).
#'    @param reference_observations: Indices of reference observations.
#'                                   Only these are used in the average.
#'
#' Returns:
#'    @return Expression with the average gene expression in the reference 
#'            observations removed.
average_over_ref <- function(average_data, ref_observations){

    logging::logdebug(paste("::average_over_ref:Start", sep=""))
    average_reference_obs <- average_data[, ref_observations]
    if (length(ref_observations) > 1){
        average_reference_obs <- rowMeans(average_data[, ref_observations],
                                          na.rm=TRUE)
    }
    return(sweep(average_data, MARGIN=1, average_reference_obs, FUN="-"))
}

#' Center data and threshold (both negative and postive values)
#'
#' Args:
#'    @param center_data: Matrix to center
#'    @param threshold: Values will be required to be with -/+1 * 
#'                      threshold after centering.
#' Returns:
#'    @return Centered and thresholded matrix
center_with_threshold <- function(center_data, threshold){

    logging::logdebug(paste("::center_with_threshold:Start", sep=""))
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

    logging::logdebug(paste("::check_arguments:Start", sep=""))
    # Require an input matrix
    #     if (!("input_matrix" %in% names(arguments)) ||
    #          (arguments$input_matrix == "")){
    #        logging::logerror(paste(":: --data_matrix: Please enter a path ",
    #                                "to the data matrix. This file should be ",
    #                                "comma delimited (rows: genes by ",
    #                                "columns: samples).",
    #                                sep=""))
    #    }
    # Require the name of a output pdf file
    if (!( "pdf_file" %in% names(arguments)) || (arguments$pdf_file == "")){
        logging::logerror(paste(":: --pdf: Please enter a file path to ",
                                "save the heatmap.",
                                 sep=""))
    }
    # Require a gene order file
    #    if (arguments$gene_order != ""){
    #        logging::logerror(paste(":: --gene_order: A file was given for ",
    #                                "the gene order, genes will be reordered ",
    #                                "according to this file.",
    #                                sep=""))
    #    }
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
        logwarn(paste(":: --reference_observations: No reference samples were",
                      " given, the average of the samples will be used.",
                      sep=""))
    }
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
infer_cnv <- function(data, gene_order, cutoff, reference_obs,
                      window_length, max_centered_threshold,
                      noise_threshold, pdf_path){

    logging::logdebug(paste("::infer_cnv:Start", sep=""))

    # Remove any gene without position information
    remove_by_position <- -1 * which(gene_order[2] + gene_order[3] == 0)
    gene_order <- gene_order[remove_by_position,]
    data <- data[remove_by_position,]
    logging::loginfo(paste("::infer_cnv:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),".", sep=""))
    logging::logdebug(paste("::infer_cnv:Removed indices:"))
    logging::logdebug(paste(remove_by_position * -1, collapse=","))

    # Reduce by cutoff
    keep_gene_indices <- infercnv::above_cutoff(data, cutoff)
    if (!is.null(keep_gene_indices)){
        data <- data[keep_gene_indices,]
        gene_order <- gene_order[keep_gene_indices,]
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff, ",
                      "new dimensions (r,c) = ",
                      paste(dim(data), collapse=","), sep=""))
        logging::logdebug(paste("::infer_cnv:Keeping indices.", sep=""))
        logging::logdebug(paste(keep_gene_indices, collapse=","))
    } else {
        logging::loginfo(paste("::infer_cnv:Reduce by cutoff.", sep=""))
        logwarn(paste("::infer_cnv::No indicies left to keep.",
                      " Stoppping."))
        stop()
    }

    # Order data by genomic region
    data <- data[with(gene_order, order(chr,start,stop)),]
    # Order gene order the same way
    chr_order <- gene_order[with(gene_order, order(chr,start,stop)),][1]
    gene_order <- NULL

    # Center data (automatically ignores zeros)
    data <- infercnv::center_with_threshold(data, max_centered_threshold)
    logging::loginfo(paste("::infer_cnv:Outlier removal.", sep=""))

    # Smooth the data with gene windows
    data_smoothed <- infercnv::smooth_window(data, window_length)
    data <- NULL
    logging::loginfo(paste("::infer_cnv:Smoothed data.", sep=""))

    # Remove average reference
    data_smoothed <- infercnv::average_over_ref(average_data=data_smoothed,
                                            ref_observations=reference_obs)
    logging::loginfo(paste("::infer_cnv:Remove average.", sep=""))

    # Remove Ends
    logging::loginfo(chr_order)
    for (chr in unlist(unique(chr_order))){
        logging::loginfo(paste("::infer_cnv:Remove tail contig ",
                               chr,
                               ".",
                               sep=""))
        data_smoothed <- infercnv::remove_tails(data_smoothed,
                                                which(chr_order == chr),
                                                (window_length - 1) / 2)
    }

    # Remove noise
    data_smoothed <- infercnv::remove_noise(reference_obs,
                                  data_smoothed,
                                  noise_threshold)
    logging::loginfo(paste("::infer_cnv:Remove noise.", sep=""))

    # Plot and write data
    logging::loginfo(paste("::infer_cnv:Drawing plots to file:",
                  pdf_path, sep=""))
    logging::loginfo(paste("::infer_cnv:Current data dimensions (r,c)=",
                           paste(dim(data_smoothed), collapse=","), sep=""))
    infercnv::plot_cnv(data_smoothed,
                       paste(as.vector(as.matrix(chr_order))),
                       reference_obs,
                       pdf_path)
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
plot_cnv <- function(plot_data, contigs, reference_idx, pdf_path){

    logging::loginfo(paste("::plot_cnv:Start", sep=""))
    logging::logdebug(paste("::plot_cnv:Current data dimensions (r,c)=",
                   paste(dim(plot_data), collapse=","), sep=""))
    logging::logdebug(paste(head(plot_data[1]), collapse=","))

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

    # Define data for for row conditions, will be used to create a colorstack
    row.conditions <- contigs
    names(row.conditions) <- contigs

    # Column seperation by contig and label axes with only one instance of name
    contig_tbl <- table(contigs)
    col_sep <- cumsum(contig_tbl)
    col_sep <- col_sep[-1 * length(col_sep)]
    contig_labels <- c()
    for(contig_name in names(contig_tbl)){
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
    plot_data <- t(plot_data)

    # Plot Observation Samples
    obs_data <- plot_data
    if(!is.null(ref_idx)){
        obs_data <- plot_data[-1 * ref_idx,,drop=FALSE]
        #TODO there will always be more obs, flip this with ref.
        if(nrow(obs_data) == 1){
                plot_data <- rbind(obs_data, obs_data)
                row.names(obs_data) <- c("", row.names(obs_data)[1])
        }
    }

    par(fig=c(0,1,0,1), new=FALSE)
    heatmap.2(obs_data,
        main="Copy Number Variation Inference",
        ylab="Cell",
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
        col="cm.colors",
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
        # heatmap2 requires a 2 x 2 matrix, so with one reference
        # I just duplicate the row and hid the second name so it
        # visually looks like it is just taking up the full realestate.
        plot_data <- plot_data[ref_idx,, drop=FALSE]
        if(nrow(plot_data) == 1){
            plot_data <- rbind(plot_data, plot_data)
            row.names(plot_data) <- c("",row.names(plot_data)[1])
        }
        # Print controls
        par(fig=c(0,1,0,1), new=TRUE)
        heatmap.2(plot_data,
            main=NA,
            ylab=NA,
            xlab=NA,
            key=FALSE,
            labCol=rep("", nrow(plot_data)),
            notecol="black",
            trace="none",
            dendrogram="none",
            Colv=FALSE,
            cexRow=0.8,
            scale="none",
            col="cm.colors",
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
#'                 against the cutoff.
#'    @param cutoff: Threshold to be above to be kept.
#'
#' Returns:
#'    @return Returns a vector of row indicies to keep (are above the cutoff).
above_cutoff <- function(data, cutoff){

    logging::logdebug(paste("::above_cutoff:Start", sep=""))
    average_gene <- log2(rowMeans( ( (2 ^ data) - 1), na.rm=TRUE) + 1)
    logging::loginfo(paste("::infer_cnv:Averages.", sep=""))
    logging::logdebug(paste(average_gene, collapse=","))
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
    keep_genes <- which(row.names(genomic_position) %in%
                        row.names(data))
    if(length(keep_genes)){
        return(data[keep_genes,])
    } else {
        return(NULL)
    }
}

#' Remove values that are too close to the average and are considered noise.
#'
#' Args:
#'    @param ref: Indices of the samples / observation that are reference.
#'        These will not be removed or changed.
#'    @param smooth_matrix: A matrix of values, smoothed, and with average 
#'                          reference removed.
#'    @param threshold: The amount of difference a value must be from the
#'                      reference before the value can be kept and not 
#'                      removed as noise.
#' Returns:
#'    @return Denoised matrix
remove_noise <- function(ref, smooth_matrix, threshold){

    logging::logdebug(paste("::remove_noise:Start.", sep=""))
    if (length(ref) == ncol(smooth_matrix)){
        ref <- c()
    }
    denoised_matrix <- smooth_matrix
    denoised_matrix[which(abs(denoised_matrix) < threshold)] <- 0
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
#'    @param chr: Indices of the chr in which the tails are to be removed.
#'    @param tail_length: Length of the tail to remove on both ends of the
#'                        chr indices.
#' Returns:
#'    @return Matrix with tails of chr set to 0.
remove_tails <- function(smooth_matrix, chr, tail_length){

    logging::logdebug(paste("::remove_tails:Start.", sep=""))
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
#' Tail so of the averages usea window length that is truncated to
#' available data.
#'
#' Args:
#'    @param data: Data matrix to smooth (columns=observations)
#'    @param window_length: Length of window to use for the moving average.
#'        Should be a positive, odd integer.
#'
#' Returns:
#'    @return Matrix with columns smoothed with a simple moving average.
smooth_window <- function(data, window_length){

    if (window_length < 2){
        return(data)
    }
    if (window_length > nrow(data)){
        return(data)
    }
    logging::logdebug(paste("::smooth_window:Start.", sep=""))
    data_sm <- data.frame(matrix(rep( NA, nrow(data) * ncol(data)),
                          nrow=nrow(data),
                          ncol=ncol(data)))
    data_sm <- apply(data,
                     2,
                     infercnv::smooth_window_helper,
                     window_length=window_length)
    tails <- (window_length - 1) / 2
    for (obs in 1:ncol(data)){
        for (tail_end in 1:tails){
            bounds <- tail_end - 1
            end_tail <- nrow(data_sm) - bounds
            data_sm[tail_end, obs] <- mean(data[, obs][(tail_end - bounds):
                                                       (tail_end + bounds)],
                                                       na.rm=TRUE)
            data_sm[end_tail, obs] <- mean(data[, obs][(end_tail - bounds):
                                                       (end_tail + bounds)],
                                                       na.rm=TRUE)
        }
    }
    return(data_sm)
}

#' Smooth vector of values over the given window length.
#'
#' Args:
#'    @param row_data: Vector of data to smooth with a moving average.
#'    @param window_length: Length of the window for smoothing.
#'        Must be and odd, positive, integer.
#' Returns:
#'    @return Vector of values smoothed with a moving average.
smooth_window_helper <- function(row_data, window_length){

    logging::logdebug(paste("::smooth_window_helper:Start.", sep=""))
    return(filter(row_data, rep(1 / window_length, window_length), sides=2))
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

    #pargs <- optparse::add_option(pargs, c("--data_matrix"),
    #                    type="character",
    #                    action="store",
    #                    dest="input_matrix",
    #                    metavar="Input_data_matrix",
    #                    help=paste("A file path is expected. Expression",
    #                               "matrix (row:genes X column:samples),",
    #                               "assumed to be log2(TPM+1)."))

    #pargs <- optparse::add_option(pargs, c("--gene_order"),
    #                    type="character",
    #                    default="",
    #                    action="store",
    #                    dest="gene_order",
    #                    metavar="Input_gene_order",
    #                    help=paste("A file path is expected. Tab delimited",
    #                               "file of chromosome, start, and stop for",
    #                               "each gene in the data matrix."))

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

    pargs <- optparse::add_option(pargs, c("--window"),
                        type="integer",
                        default=101,
                        action="store",
                        dest="window_length",
                        metavar="Window_Lengh",
                        help=paste("Window length for the smoothing.",
                                   "[Default %default]"))

    args_parsed <- optparse::parse_args(pargs, positional_arguments=2)
    args <- args_parsed$options
    args["input_matrix"] <- args_parsed$args[1]
    args["gene_order"] <- args_parsed$args[2]

    # Check arguments
    infercnv::check_arguments(args)

    # Set up logging file
    optparse::basicConfig(level=args$log_level)
    if (!is.na(args$log_file)){
        optparse::addHandler(optparse::writeToFile,
                             file=args$log_file,
                             level=args$log_level)
    }

    # Manage inputs
    logging::loginfo(paste("::Reading data matrix.", sep=""))
    expression_data <- read.table(args$input_matrix)
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
    expression_data <- infercnv::order_reduce(data=expression_data,
                                   genomic_position=input_gene_order)
    if(is.null(expression_data)){
        error_message <- paste("None of the genes in the expression data",
                              "matched the genes in the reference genomic",
                              "position file. Analysis Stopped.")
        stop(error_message)
    }

    # Run CNV inference
    infercnv::infer_cnv(data=expression_data,
                        gene_order=input_gene_order,
                        cutoff=args$cutoff,
                        reference_obs=input_reference_samples,
                        window_length=args$window_length,
                        max_centered_threshold=args$max_centered_expression,
                        noise_threshold=args$magnitude_filter,
                        pdf_path=args$pdf_file)
}
