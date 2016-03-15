#!/usr/bin/env Rscript

library(RColorBrewer)
library(gplots)
library(optparse)
library(logging)

# Constants
C_PROGRAM_NAME = "scCNV"
C_LEVEL_CHOICES = names(loglevels)

#' Remove the average of the genes of the reference observations from all observations' expression.
#'
#' Args:
#'    @param average_data: Matrix containing the data to remove average from (this includes the reference observations.
#'    @param reference_observations: Indices of reference observations. Only these are used in the average.
#'
#' Returns:
#'    @return Expression with the average gene expression in the reference observations removed.
average_over_reference <- function(average_data, reference_observations){

    logdebug(paste(C_PROGRAM_NAME, "::average_over_reference:Start", sep=""))
    average_reference_obs <- average_data[, reference_observations[1]]
    if (length(reference_observations) > 1){
        average_reference_obs <- rowMeans(average_data[, reference_observations], na.rm=TRUE)
    }
    return(sweep(average_data, MARGIN=1, average_reference_obs, FUN="-"))
}

#' Center data and threshold (both negative and postive values)
#'
#' Args:
#'    @param center_data: Matrix to center
#'    @param threshold: Values will be required to be with -/+1 * threshold after centering.
#' Returns:
#'    @return Centered and thresholded matrix
center_with_threshold <- function(center_data, threshold){

    logdebug(paste(C_PROGRAM_NAME, "::center_with_threshold:Start", sep=""))
    # Center data (automatically ignores zeros)
    center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
    # Cap values between threshold and -threshold and recenter
    center_data[center_data > threshold] <- threshold
    center_data[center_data < (-1 * threshold)] <- -1 * threshold
    center_data <- center_data - rowMeans(center_data, na.rm=TRUE)
    return(center_data)
}

#' Check arguments and make sure the user input meet certain additional requirements.
#'
#' Args:
#'    @param arguments: Parsed arguments from user.
#' Returns:
#'    @return No return.
check_arguments <- function(arguments){

    logdebug(paste(C_PROGRAM_NAME, "::check_arguments:Start", sep=""))
    # Require an input matrix
    #     if (!("input_matrix" %in% names(arguments)) || (arguments$input_matrix == "")){
    #        logerror(paste(C_PROGRAM_NAME,
    #                       ":: --data_matrix: Please enter a path to the data matrix. ",
    #                       "This file should be comma delimited (rows: genes by columns: samples).",
    #                       sep=""))
    #    }
    # Require the name of a output pdf file
    if (!( "pdf_file" %in% names(arguments)) || (arguments$pdf_file == "")){
        logerror(paste(C_PROGRAM_NAME,
                       ":: --pdf: Please enter a file path to save the heatmap.",
                       sep=""))
    }
    # Require a gene order file
    #    if (arguments$gene_order != ""){
    #        logerror(paste(C_PROGRAM_NAME,
    #                       ":: --gene_order: A file was given for the gene order, ",
    #                       "genes will be reordered according to this file.",
    #                       sep=""))
    #    }
    # Require the cut off to be above 0
    if (arguments$cutoff < 0){
        logerror(paste(C_PROGRAM_NAME,
                       ":: --cutoff: Please enter a value greater or equal to zero for the cut off.",
                       sep=""))
    }
    # Require the logging level to be one handled by logging
    if (!(arguments$log_level %in% C_LEVEL_CHOICES)){
        logerror(paste(C_PROGRAM_NAME,
                       ":: --log_level: PLease use a logging level given here: ",
                       C_LEVEL_CHOICES,
                       collapse=",", sep=""))
    }
    # Warn that an average of the samples is used in the absence of normal / reference samples
    if (is.na(arguments$reference_observations)){
        logwarn(paste(C_PROGRAM_NAME,
                      ":: --reference_observations: No reference samples were given, ",
                      "the average of the samples will be used.",
                      sep=""))
    }
}

#' Infer CNV changes given a matrix of RNASeq counts.
#' Output a pdf and matrix of final values.
#'
#' Args:
#'    @param data: Expression matrix (genes X samples), assumed to be log2(TPM+1) .
#'    @param gene_order: Ordering of the genes (data's rows) according to their genomic location
#'        To include all genes use 0.
#'    @param cutoff: Cut-off for the average expression of genes to be used for CNV inference.
#'    @param reference_obs: Indices of the subset of samples ( data's columns ) that should be used as references
#'      If not given, the average of all samples will be the reference.
#'    @param window_length: Length of the window for the moving average (smoothing).
#'        Should be an odd integer.
#'    @param max_centered_threshold: The maximum value a a value can have after centering. Also sets a lower bound of -1 * this value.
#'    @param noise_threshold: The minimum difference a value can be from the average reference in order for it not to be removed as noise.
#'    @param pdf_path: The path to what to save the pdf as. The raw data is also written to this path but with the extension .txt .
#'
#' Returns:
#'    @return No return.
infer_cnv <- function(data, gene_order, cutoff, reference_obs,
                      window_length, max_centered_threshold, noise_threshold, pdf_path){

    logdebug(paste(C_PROGRAM_NAME, "::infer_cnv:Start", sep=""))

    # Remove any gene without position information
    remove_by_position <- -1 * which(gene_order[2] + gene_order[3] == 0)
    gene_order <- gene_order[remove_by_position, ]
    data <- data[remove_by_position, ]
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Reduction from positional data, new dimensions (r,c) = ",
                  paste(dim(data), collapse=","),".", sep=""))
    logdebug(paste(C_PROGRAM_NAME, "::infer_cnv:Removed indices:"))
    logdebug(paste(remove_by_position*-1, collapse=","))

    # Reduce by cutoff
    keep_gene_indices <- above_cutoff(data, cutoff)
    if (!is.null(keep_gene_indices)){
        data <- data[keep_gene_indices, ]
        gene_order <- gene_order[keep_gene_indices, ]
        loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Reduce by cutoff, new dimensions (r,c) = ",
                      paste(dim(data), collapse=","), sep=""))
        logdebug(paste(C_PROGRAM_NAME, "::infer_cnv:Keeping indices.", sep=""))
        logdebug(paste(keep_gene_indices, collapse=","))
    } else {
        loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Reduce by cutoff.", sep=""))
        logwarn(paste(C_PROGRAM_NAME, "::infer_cnv::NO indicies left to keep. Stoppping."))
        stop()
    }

    # Order data by genomic region
    data <- data[with(gene_order, order(chr,start,stop)),]
    # Order gene order the same way
    chr_order <- gene_order[with(gene_order, order(chr,start,stop)),][1]
    gene_order = NULL

    # Center data (automatically ignores zeros)
    data <- center_with_threshold(data, max_centered_threshold)
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:After outlier removal.", sep=""))

    # Smooth the data with gene windows
    data_smoothed <- smooth_window(data, window_length)
    print(data_smoothed)
    data = NULL
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Smoothed data.", sep=""))

    # Remove average reference
    data_smoothed <- average_over_reference(average_data=data_smoothed,
                                            reference_observations=reference_obs)
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Remove average.", sep=""))

    # Remove Ends
    loginfo(chr_order)
    for (chr in unlist(unique(chr_order))){
        loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Remove tail contig ", chr, ".", sep=""))
        data_smoothed <- remove_tails(data_smoothed, which(chr_order==chr), (window_length - 1) / 2)
    }

    # Remove noise
    data_smoothed <- remove_noise(reference_obs, data_smoothed, noise_threshold)
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Remove noise.", sep=""))

    # Plot and write data
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Drawing plots to file:", pdf_path, sep=""))
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Current data dimensions (r,c) = ",
                  paste(dim(data_smoothed), collapse=","), sep=""))
    plot_cnv(data_smoothed, paste(as.vector(as.matrix(chr_order))), pdf_path)
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Writing final data to",
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
plot_cnv <- function(plot_data, contigs, pdf_path){

    print(contigs)
    # Plot heatmap 
    # Cluster CNV patterns
    # TODO change to using Sam's tool.
    loginfo(paste(C_PROGRAM_NAME, "::plot_cnv:Start", sep=""))
    logdebug(paste(C_PROGRAM_NAME, "::plot_cnv:Current data dimensions (r,c) = ", 
                   paste(dim(plot_data), collapse=","), sep=""))
    #plot_data <- t(plot_data)
    logdebug(paste(head(plot_data[1]), collapse=","))

    # Heatmaps using CompoHeatMap
    # Rows observations, Columns CHR
    # Observation clustering
    pdf(pdf_path, useDingbats=FALSE)
    # Contigs
    ## Define the colors corresponding to contigs
    ## Make sure there are no spaces in the name
    unique_contigs <- unique(contigs)
    ct.colors <- colorRampPalette(brewer.pal(8,"Accent"))(length(unique_contigs))
    #names(ct.colors) <- unique_contigs
    #row.names(plot_data) <- contigs

    # Row seperation based on reference
    ref_blocks <- unique(unlist(tapply(reference_idx, c(0,cumsum(diff(refernce_idx)!=1)),range)))
    ref_sep = setdiff( c(ref_blocks-1, ref_blocks+1), ref_blocks)

    ## Define data for row cluster size => will be used to create a vertical barplot
    clu.sz <- as.integer(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
    names(clu.sz) <- contigs
    clu.sz.colors <- color.palette(c("white", "purple"), 20)(20)

    ## Define data for for row conditions => will be used to create a colorstack
    row.conditions <- contigs
    names(row.conditions) <- contigs

    ## Test building each plot component and composing them all in one command.
    ## The supporting components of the plot can be accessed in the returned list.
    #p.complete.l=create.gg.hmap.w.barps(hmap.dat=plot_data,
    #    barp.dat.v=clu.sz,
    #    colorstack.dat.l=list(row.conditions),
    #    colorstack.fill.l=list(ct.colors),
    #    widths=c(2,1,27),
    #    plot.title="Copy Number Variation Inference",
    #    leg.title.l=c("Contig"),
    #    print.plot=FALSE,
    #    ## heatmap-specific arguments follow
    #    dend.c.ord=NULL,
    #    dend.r.ord=NULL,
    #    dend.r=FALSE,
    #    dend.c=FALSE,
    #    cut.frac.r.h=NULL,
    #    cut.frac.c.h=NULL,
    #    hmap.col=get.hmap.col(range.val=range(plot_data),mid.val=0),
    #    leg.title="Centered RNA-Seq Intensity",
    #    norm.rows=FALSE)

    heatmap.2(plot_data,
        main="Copy Number Variation Inference",
        notecol="black",
        density.info="histogram",
        trace="none",
        dendrogram="row",
        ColV=FALSE,
        scale="none",
        col="heat.colors",
        # Seperate by contigs
        colsep=table(contigs),
        # Seperate by reference / not reference
        rowsep=ref_sep,
        sepcolor="white",
        sepwdith=c(0.5,0.5),
        # Color by contigs
        ColSideColors=ct.colors[ct.colors])
    dev.off()
}

#' Return the indices of the rows that average above the cut off
#'
#' Args:
#'    @param data: Data to measure the average row and evaluate against the cutoff.
#'    @param cutoff: Threshold to be above to be kept.
#'
#' Returns:
#'    @return Returns a vector of row indicies to keep (are above the cutoff).
above_cutoff <- function(data, cutoff){

    logdebug(paste(C_PROGRAM_NAME, "::above_cutoff:Start", sep=""))
    average_gene <- log2(rowMeans(((2^data) - 1), na.rm=TRUE) + 1)
    loginfo(paste(C_PROGRAM_NAME, "::infer_cnv:Averages.", sep=""))
    logdebug(paste(average_gene, collapse=","))
    # Find averages above a certain threshold
    indicies <- which(average_gene > cutoff)
    if (length(indicies) > 0){
        return(indicies)
    } else {
        return(NULL)
    }
}

#' Remove values that are too close to the average and are considered noise.
#'
#' Args:
#'    @param ref: Indices of the samples / observation that are reference.
#'        These will not be removed or changed.
#'    @param smooth_matrix: A matrix of values, smoothed, and with average reference removed.
#'    @param threshold: The amount of difference a value must be from the reference before
#'        the value can be kept and not removed as noise.
#' Returns:
#'    @return Denoised matrix
remove_noise <- function(ref, smooth_matrix, threshold){

    logdebug(paste(C_PROGRAM_NAME, "::remove_noise:Start.", sep=""))
    if (length(ref) == ncol(smooth_matrix)){
        ref <- c()
    }
    denoised_matrix <- smooth_matrix
    denoised_matrix[which(abs(denoised_matrix) < threshold)] <- 0
    if (length(ref) > 0){
        denoised_matrix[, ref] = smooth_matrix[, ref]
    }
    return(denoised_matrix)
}

#' Remove the tails of values of a specific chromosome.
#' The smooth_matrix values are expected to be in genomic order.
#'
#' Args:
#'    @param smooth_matrix: Smoothed values in genomic order.
#'    @param chr: Indices of the chr in which the tails are to be removed.
#'    @param tail_length: Length of the tail to remove on both ends of the chr indices.
#' Returns:
#'    @return Matrix with tails of chr set to 0.
remove_tails <- function(smooth_matrix, chr, tail_length){

    logdebug(paste(C_PROGRAM_NAME, "::remove_tails:Start.", sep=""))
    if (tail_length < 1 ){
        return(smooth_matrix)
    }
    if (length(chr) < (tail_length * 2)){
         smooth_matrix[chr, ] <- 0
    } else {
         chr_length <- length(chr)
         smooth_matrix[chr[1:tail_length], ] <- 0
         smooth_matrix[chr[((chr_length+1) - tail_length):chr_length], ] <- 0
    }
    return(smooth_matrix)
}

#' Smooth a matrix by column using a simple moving average.
#' Tail so of the averages usea window length that is truncated to available data.
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
    logdebug(paste(C_PROGRAM_NAME, "::smooth_window:Start.", sep=""))
    data_smoothed <- data.frame(matrix(rep( NA, nrow(data) * ncol(data)), nrow=nrow(data), ncol=ncol(data)))
    data_smoothed <- apply(data, 2, smooth_window_helper, window_length=window_length)
    tails <- (window_length - 1) / 2
    len_data <- nrow(data)
    for (obs in 1:ncol(data)){
        for (tail_end in 1:tails){
            bounds <- tail_end - 1
            end_tail <- nrow( data_smoothed ) - bounds
            data_smoothed[tail_end, obs] <- mean(data[, obs][(tail_end - bounds):(tail_end + bounds)], na.rm=TRUE)
            data_smoothed[end_tail, obs] <- mean(data[, obs][(end_tail - bounds):(end_tail + bounds)], na.rm=TRUE)
        }
    }
    return(data_smoothed)
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

    logdebug(paste(C_PROGRAM_NAME, "::smooth_window_helper:Start.", sep=""))
    return(filter(row_data, rep(1 / window_length, window_length), sides=2))
}


# Calls if called by commandline but not sourced
if (identical(environment(),globalenv()) &&
    !length(grep("^source\\(", sys.calls()))){

    # Command line arguments
    pArgs <- OptionParser(usage="%prog [options] --pdf pdf_file data_matrix genomic_positions")

    pArgs <- add_option(pArgs, c("--cutoff"),
                        type="integer",
                        default=0,
                        action="store",
                        dest="cutoff",
                        metavar="Cutoff",
                        help="A number >= 0 is expected. A cut off for the average expression of genes to be used for CNV inference. [Default %default]")

    #pArgs <- add_option(pArgs, c("--data_matrix"),
    #                    type="character",
    #                    action="store",
    #                    dest="input_matrix",
    #                    metavar="Input_data_matrix",
    #                    help="A file path is expected. Expression matrix (row:genes X column:samples), assumed to be log2(TPM+1).")

    #pArgs <- add_option(pArgs, c("--gene_order"),
    #                    type="character",
    #                    default="",
    #                    action="store",
    #                    dest="gene_order",
    #                    metavar="Input_gene_order",
    #                    help="A file path is expected. Tab delimited file of chromosome, start, and stop for each gene in the data matrix.")

    pArgs <- add_option(pArgs, c("--log"),
                        type="character",
                        action="store",
                        default=NA,
                        dest="log_file",
                        metavar="Log",
                        help="File for logging. If not given, logging will occur to console. [Default %default]")

    pArgs <- add_option(pArgs, c("--log_level"),
                        type="character",
                        action="store",
                        default="INFO",
                        dest="log_level",
                        metavar="LogLevel",
                        help=paste("Logging level. Valid choices are ", C_LEVEL_CHOICES, "[Default %default]", collapse=","))

    pArgs <- add_option(pArgs, c("--noise_filter"),
                        type="integer",
                        default=0,
                        action="store",
                        dest="magnitude_filter",
                        metavar="Magnitude_Filter",
                        help="A value must be atleast this much more or less than the reference to be plotted [Default %default].")

    pArgs <- add_option(pArgs, c("--max_centered_expression"),
                        type="integer",
                        default=3,
                        action="store",
                        dest="max_centered_expression",
                        metavar="Max_centered_expression",
                        help="This value and -1 * this value are used as the maximum value expression that can exist after centering data. If a value is outside of this range, it is truncated to be within this range [Default %default].")

    pArgs <- add_option(pArgs, c("--pdf"),
                        type="character",
                        action="store",
                        dest="pdf_file",
                        metavar="Output_PDF_Visualization",
                        help="Output PDF for visualizing the analysis. [Default %default][REQUIRED]")

    pArgs <- add_option(pArgs, c("--ref"),
                        type="character",
                        default=NA,
                        action="store",
                        dest="reference_observations",
                        metavar="Input_reference_observations",
                        help="Tab delimited integers are expected. Indices of the subset of samples ( data's columns ) that should be used as references if not given, the average of all samples will be the reference. [Default %default]")

    pArgs <- add_option(pArgs, c("--window"),
                        type="integer",
                        default=101,
                        action="store",
                        dest="window_length",
                        metavar="Window_Lengh",
                        help="Window length for the smoothing. [Default %default]")

    args_parsed <- parse_args(pArgs, positional_arguments=2)
    args <- args_parsed$options
    args["input_matrix"] <- args_parsed$args[1]
    args["gene_order"] <- args_parsed$args[2]

    # Check arguments
    check_arguments(args)

    # Set up logging file
    basicConfig(level=args$log_level)
    if (!is.na(args$log_file)){
        addHandler(writeToFile, file=args$log_file, level=args$log_level)
    }

    # Manage inputs
    loginfo(paste(C_PROGRAM_NAME, "::Reading data matrix.", sep=""))
    expression_data <- read.table(args$input_matrix)
    loginfo(paste(C_PROGRAM_NAME, "Original matrix dimensions (r,c) = ", paste(dim(expression_data), collapse=",")))

    # Default the gene order to the given order
    # If given a file read and sort and order.
    input_gene_order <- seq(1, nrow(expression_data), 1)
    if (args$gene_order != ""){
        input_gene_order <- read.table(args$gene_order, col.names=c("chr", "start", "stop"))
    }
    loginfo(paste(C_PROGRAM_NAME, "::Reading gene order.", sep=""))
    logdebug(paste(head(args$gene_order[1]), collapse=","))

    # Default the reference samples to all
    input_reference_samples <- colnames(expression_data)
    if (!is.na(args$reference_observations)){
        input_reference_samples <- unique(unlist(split(args$reference_observations,",")))
    }

    # Make sure the given reference samples are in the matrix.
    if (length(input_reference_samples) != length(intersect(input_reference_samples, colnames(expression_data)))){
        missing_reference_sample <- setdiff(input_reference_samples, colnames(expression_data))
        error_message = paste(C_PROGRAM_NAME,
                              "Please make sure that all the reference sample names match a sample in your data matrix.\nAttention to: ",
                              paste(missing_reference_sample, collapse=","))
        logerror(error_message)
        stop(error_message)
    }

    # Run CNV inference
    infer_cnv(data=expression_data,
              gene_order=input_gene_order,
              cutoff=args$cutoff,
              reference_obs=input_reference_samples,
              window_length=args$window_length,
              max_centered_threshold=args$max_centered_expression,
              noise_threshold=args$magnitude_filter,
              pdf_path=args$pdf_file)
}
