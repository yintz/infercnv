#!/usr/bin/env Rscript


options(error = function() traceback(2))

#' The infercnv Class
#'
#' An infercnv object encapsulates the expression data and gene chromosome ordering information
#' that is leveraged by infercnv for data exploration.  The infercnv object is passed among the
#' infercnv data processing and plotting routines.
#'
#' Slots in the infercnv object include:
#'
#' @slot expr.data  <matrix>  the count or expression data matrix
#'
#' @slot gene_order  <data.frame> chromosomal gene order
#'
#' @slot reference_grouped_cell_indices <list>  mapping [['group_name']] to c(cell column indices) for reference (normal) cells
#'
#' @slot observation_grouped_cell_indices <list> mapping [['group_name']] to c(cell column indices) for observation (tumor) cells
#' 

infercnv <- methods::setClass(
                         "infercnv",
                         slots = c(
                             expr.data = "matrix",
                             gene_order= "data.frame",
                             reference_grouped_cell_indices = "list",
                             observation_grouped_cell_indices = "list") )




#' @title CreateInfercnvObject
#'
#' Creation of an infercnv object. This requires the following inputs:
#'     raw_counts_matrix : the matrix file of genes (rows) vs. cells (columns) containing the raw counts
#'     gene_order_file : data file containing the positions of each gene along each chromosome in the genome.
#'     annotations_file : a description of the cells, indicating the cell type classifications
#'     ref_group_names : a vector containing the classifications of the reference (normal) cells to use for infering cnv
#'
#' A more detailed description of each input is provided below:
#'
#' The raw_counts_matrix:
#'
#'           MGH54_P16_F12 MGH53_P5_C12 MGH54_P12_C10 MGH54_P16_F02 MGH54_P11_C11  ...
#' DDX11L1     0.0000000     0.000000      0.000000      0.000000     0.0000000
#' WASH7P      0.0000000     2.231939      7.186235      5.284944     0.9650009
#' FAM138A     0.1709991     0.000000      0.000000      0.000000     0.0000000
#' OR4F5       0.0000000     0.000000      0.000000      0.000000     0.0000000
#' OR4F29      0.0000000     0.000000      0.000000      0.000000     0.0000000
#' ...
#'
#' The gene_order_file, contains chromosome, start, and stop position for each gene, tab-delimited:
#'
#'          chr  start   stop
#' DDX11L1 chr1  11869  14412
#' WASH7P  chr1  14363  29806
#' FAM138A chr1  34554  36081
#' OR4F5   chr1  69091  70008
#' OR4F29  chr1 367640 368634
#' OR4F16  chr1 621059 622053
#' ...
#' 
#' The annotations_file, containing the cell name and the cell type classification, tab-delimited.
#'
#'             V1                   V2
#' 1 MGH54_P2_C12 Microglia/Macrophage
#' 2 MGH36_P6_F03 Microglia/Macrophage
#' 3 MGH53_P4_H08 Microglia/Macrophage
#' 4 MGH53_P2_E09 Microglia/Macrophage
#' 5 MGH36_P5_E12 Oligodendrocytes (non-malignant)
#' 6 MGH54_P2_H07 Oligodendrocytes (non-malignant)
#' ...
#' 179  93_P9_H03 malignant
#' 180 93_P10_D04 malignant
#' 181  93_P8_G09 malignant
#' 182 93_P10_B10 malignant
#' 183  93_P9_C07 malignant
#' 184  93_P8_A12 malignant
#' ...
#'
#'
#' and the ref_group_names vector might look like so:  c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")

CreateInfercnvObject <- function(raw_counts_matrix, gene_order_file, annotations_file, ref_group_names, delim="\t") {
    
    # input expression data
    raw.data <- read.table(raw_counts_matrix, sep=delim, header=TRUE, row.names=1, check.names=FALSE)    

    # get gene order info
    gene_order <- read.table(gene_order_file, header=FALSE, row.names=1, sep="\t")
    names(gene_order) <- c(C_CHR, C_START, C_STOP)

    # read annotations file
    input_classifications <- read.table(annotations_file, header=FALSE, row.names=1, sep=delim, stringsAsFactors=FALSE)

    # make sure all reference samples are accounted for:
    if (! all( rownames(input_classifications) %in% colnames(raw.data)) ) {
        
        missing_cells <- rownames(input_classifications)[ ! ( rownames(input_classifications) %in% colnames(raw.data) ) ]

        error_message <- paste("Please make sure that all the annotated cell ",
                               "names match a sample in your data matrix. ",
                               "Attention to: ",
                               paste(missing_cells, collapse=","))
        stop(error_message)
    }

    # restrict expression data to the annotated cells.
    raw.data <- raw.data[,colnames(raw.data) %in% rownames(input_classifications)]
        
    # reorder cell classifications according to expression matrix column names
    input_classifications <- input_classifications[order(match(row.names(input_classifications), colnames(raw.data))), , drop=FALSE]

    
    # get indices for reference cells
    ref_group_cell_indices = list()
    for (name_group in ref_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        cell_names =  rownames(input_classifications)[cell_indices]
        ref_group_cell_indices[[ name_group ]] <- cell_indices
    }
    
    # rest of the cells are the 'observed' set.
    all_group_names <- unique(input_classifications[,1])
    obs_group_names <- setdiff(all_group_names, ref_group_names)
        
    # extract the genes indicated in the gene ordering file:
    order_ret <- .order_reduce(data=raw.data, genomic_position=gene_order)

    num_genes_removed = dim(raw.data)[1] - dim(order_ret$exp)[1]

    if (num_genes_removed > 0) {
        flog.info(paste("num genes removed taking into account provided gene ordering list: ",
                        num_genes_removed,
                        " = ",
                        num_genes_removed / dim(raw.data)[1] * 100,
                        "% removed.", sep=""))
    }
    
    raw.data <- order_ret$expr
    input_gene_order <- order_ret$order
    
    if(is.null(raw.data)) {
        error_message <- paste("None of the genes in the expression data",
                               "matched the genes in the reference genomic",
                               "position file. Analysis Stopped.")
        stop(error_message)
    }
    
    # define groupings according to the observation annotation names
    

    obs_group_cell_indices = list()
    for (name_group in obs_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        cell_names = rownames(input_classifications)[cell_indices]
        obs_group_cell_indices[[ name_group ]] <- cell_indices
    }

    raw.data = as.matrix(raw.data)
    
    object <- new(
        Class = "infercnv",
        expr.data = raw.data, #simple copy for now
        gene_order = input_gene_order,
        reference_grouped_cell_indices = ref_group_cell_indices,
        observation_grouped_cell_indices = obs_group_cell_indices)
    
    return(object)
}


# Order the data and subset the data to data in the genomic position file.
#
# Args:
# @param data Data (expression) matrix where the row names should be in
#                 the row names of the genomic_position file.
# @param genomic_position Data frame read in from the genomic position file
#
# @return Returns a matrix of expression in the order of the
#            genomic_position file. NULL is returned if the genes in both
#            data parameters do not match.
#

.order_reduce <- function(data, genomic_position){
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
    chr_levels <- unique(genomic_position[[C_CHR]])
    ret_results$order[[C_CHR]] <- factor(ret_results$order[[C_CHR]],
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


