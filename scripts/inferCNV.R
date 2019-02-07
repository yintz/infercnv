#!/usr/bin/env Rscript

# To use inferCNV via command-line interface, first install inferCNV per Wiki,
# then run a command like the following:
#
# ./scripts/inferCNV.R \
#     --output_dir=output_cli \
#     --ref_groups="Microglia/Macrophage","Oligodendrocytes (non-malignant)" \
#     --annotations_file=example/oligodendroglioma_annotations_downsampled.txt \
#     example/oligodendroglioma_expression_downsampled.counts.matrix \
#     example/gencode_downsampled.txt

# Load libraries
library(ape)
library("RColorBrewer", character.only=TRUE)
library(GMD)
library(optparse)
library(logging)
if (!require('fastcluster')) {
    warning("fastcluster library not available, using the default hclust method instead.")
}
library(infercnv)

# Logging level choices
C_LEVEL_CHOICES <- names(loglevels)
# Visualization outlier thresholding and bounding method choices
C_VIS_OUTLIER_CHOICES <- c("average_bound")
C_REF_SUBTRACT_METHODS <- c("by_mean", "by_quantiles")
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

CHR = "chr"
START = "start"
STOP = "stop"

logging::basicConfig(level='INFO') #initialize to info setting.

#' Check arguments and make sure the user input meet certain
#' additional requirements.
#'
#' Args:
#'    @param arguments: Parsed arguments from user.
#'
#' Returns:
#'    @return: Updated arguments
check_arguments <- function(arguments){

    logging::loginfo(paste("::check_arguments:Start", sep=""))
    # Require the name of a output pdf file
    if ( (!( "output_dir" %in% names(arguments))) || (arguments$output_dir == "") || (is.na(arguments$output_dir)) ) {
        logging::logerror(paste(":: --output_dir: Please enter a file path to ",
                                "save the heatmap.",
                                 sep=""))

        stop("error, no --output_dir")
    }

    # Require the cut off to be above 0
    if (arguments$cutoff < 0){
        logging::logerror(paste(":: --cutoff: Please enter a value",
                                "greater or equal to zero for the cut off.",
                                sep=""))

        stop("error, no --cutoff")
    }

    # Require the logging level to be one handled by logging
    if (!(arguments$log_level %in% C_LEVEL_CHOICES)){
        logging::logerror(paste(":: --log_level: Please use a logging level ",
                                "given here: ", C_LEVEL_CHOICES,
                                collapse=",", sep=""))
        stop("error, not recognizing log level")
    }

    # Require the visualization outlier detection to be a correct choice.
    if (!(arguments$bound_method_vis %in% C_VIS_OUTLIER_CHOICES)){
        logging::logerror(paste(":: --vis_bound_method: Please use a method ",
                                "given here: ", C_VIS_OUTLIER_CHOICES,
                                collapse=",", sep=""))
        stop("error, must specify acceptable --vis_bound_method")
    }

    if (! (arguments$ref_subtract_method %in% C_REF_SUBTRACT_METHODS) ) {
        logging::logerror(paste(":: --ref_subtract_method: acceptable values are: ",
                                paste(C_REF_SUBTRACT_METHODS, collapse=","), sep="") )
        stop("error, must specify acceptable --ref_subtract_method")
    }


    if (! (arguments$hclust_method %in% C_HCLUST_METHODS) ) {
        logging::logerror(paste(":: --hclust_method: acceptable values are: ",
                                paste(C_HCLUST_METHODS, collapse=","), sep="") )
        stop("error, must specify acceptable --hclust_method")
    }

    # # Warn that an average of the samples is used in the absence of
    # # normal / reference samples
    # if (is.null(arguments$reference_observations)){
    #     logging::logwarn(paste(":: --reference_observations: No reference ",
    #                   "samples were given, the average of the samples ",
    #                   "will be used.",
    #                   sep=""))
    # }

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

        stop(980)
    }

    # if (! is.na(suppressWarnings(as.integer(arguments$name_ref_groups)))){
    #     arguments$name_ref_groups <- list(as.integer(arguments$name_ref_groups))
    # } else {

    if (! is.na(arguments$name_ref_groups)) {
        arguments$name_ref_groups <- unlist(strsplit(arguments$name_ref_groups,","))
    }
    else {
        if(!is.null(arguments$num_ref_groups)) {
            logging::logerror(paste("::  cannot use --num_ref_groups without",
                                    "using --ref_groups as the average of all",
                                    "cells is used."))
            stop(978)
        }
    }
    return(arguments)
}

# Command line arguments
pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--output_dir directory",
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

pargs <- optparse::add_option(pargs, c("--contig_lab_size"),
                              type="integer",
                              action="store",
                              default=1,
                              dest="contig_label_size",
                              metavar="Contig_Label_Size",
                              help=paste("Used to increase or decrease the text labels",
                                         "for the X axis (contig names).",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--cutoff"),
                              type="numeric",
                              default=0,
                              action="store",
                              dest="cutoff",
                              metavar="Cutoff",
                              help=paste("A number >= 0 is expected. A cut off for",
                                         "the average expression of genes to be used",
                                         "for CNV inference (use the value before log2 transformation). [Default %default]"))

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

pargs <- optparse::add_option(pargs, c("--log_file"),
                              type="character",
                              action="store",
                              default=NA,
                              dest="log_file",
                              metavar="Log",
                              help=paste("File for logging. If not given,",
                                         "logging will occur to console.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--delim"),
                              type="character",
                              action="store",
                              default="\t",
                              dest="delim",
                              metavar="Delimiter",
                              help=paste("Delimiter for reading expression matrix",
                                        " and writing matrices output.",
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

pargs <- optparse::add_option(pargs, c("--output_dir"),
                              type="character",
                              action="store",
                              dest="output_dir",
                              metavar="Output_Directory",
                              help=paste("Output directory for analysis products.",
                                         "[Default %default][REQUIRED]"))

pargs <- optparse::add_option(pargs, c("--ref"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="reference_observations",
                              metavar="Input_reference_observations",
                              help=paste("Tab delimited characters are expected.",
                                         "Names of the subset each sample ( data's",
                                         "columns ) is part of.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--num_ref_groups"),
                              type="integer",
                              default=NULL,
                              action="store",
                              dest="num_ref_groups",
                              metavar="Number_of_reference_groups",
                              help=paste("Define a number of groups to",
                                         "make automatically by unsupervised",
                                         "clustering. This ignores annotations",
                                         "within references, but does not",
                                         "mix them with observations.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--ref_groups"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="name_ref_groups",
                              metavar="Name_of_reference_groups",
                              help=paste("Names of groups from --ref table whose cells",
                                         "are to be used as reference groups.",
                                         "[REQUIRED]"))

pargs <- optparse::add_option(pargs, c("--ref_subtract_method"),
                              type="character",
                              default="by_mean",
                              action="store",
                              dest="ref_subtract_method",
                              metavar="Reference_Subtraction_Method",
                              help=paste("Method used to subtract the reference values from the observations. Valid choices are: ",
                                         paste(C_REF_SUBTRACT_METHODS, collapse=", "),
                                         " [Default %default]"))

pargs <- optparse::add_option(pargs, c("--hclust_method"),
                              type="character",
                              default="complete",
                              action="store",
                              dest="hclust_method",
                              metavar="Hierarchical_Clustering_Method",
                              help=paste("Method used for hierarchical clustering of cells. Valid choices are: ",
                                         paste(C_HCLUST_METHODS, collapse=", "),
                                         " [Default %default]"))



pargs <- optparse::add_option(pargs,c("--obs_cluster_contig"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="clustering_contig",
                              metavar="Clustering_Contig",
                              help=paste("When clustering observation samples, ",
                                         "all genomic locations are used unless ",
                                         "this option is given. The expected value ",
                                         "is one of the contigs (Chr) in the genomic ",
                                         "positions file (case senstive). All genomic ",
                                         "positions will be plotted but only the given ",
                                         "contig will be used in clustering / group ",
                                         "creation."))

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
                              default="average_bound",
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

pargs <- optparse::add_option(pargs, c("--title"),
                              type="character",
                              default="Copy Number Variation Inference",
                              action="store",
                              dest="fig_main",
                              metavar="Figure_Title",
                              help=paste("Title of the figure.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--title_obs"),
                              type="character",
                              default="Observations (Cells)",
                              action="store",
                              dest="obs_main",
                              metavar="Observations_Title",
                              help=paste("Title of the observations matrix Y-axis.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--title_ref"),
                              type="character",
                              default="References (Cells)",
                              action="store",
                              dest="ref_main",
                              metavar="References_Title",
                              help=paste("Title of the references matrix Y-axis (if used).",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--save"),
                              type="logical",
                              action="store_true",
                              default=FALSE,
                              dest="save",
                              metavar="save",
                              help="Save workspace as infercnv.Rdata")

pargs <- optparse::add_option(pargs, c("--ngchm"),
                              type="logical",
                              action="store_true",
                              default=FALSE,
                              dest="ngchm",
                              metavar="NextGen_HeatMap",
                              help=paste("Create a Next Generation Clustered Heat Map"))

pargs <- optparse::add_option(pargs, c("--path_to_shaidyMapGen"),
                              type="character",
                              action="store",
                              default=NULL,
                              dest="path_to_shaidyMapGen",
                              metavar="Path_To_ShaidyMapGenp",
                              help=paste("This is the pathway to the java application ShaidyMapGen.jar.",
                                    "If this is not assigned, then an enviornmental variable that ",
                                    "contains the "))

pargs <- optparse::add_option(pargs, c("--gene_symbol"),
                              type="character",
                              action="store",
                              default=NULL,
                              dest="gene_symbol",
                              metavar="Gene_Symbol",
                              help=paste("The labeling type used to represent the genes in the expression",
                                   "data. This needs to be passed in order to add linkouts to the ",
                                   "genes. Possible gene label types to choose from are specified on",
                                   "the broadinstitute/inferCNV wiki and bmbroom/NGCHM-config-biobase."))

pargs <- optparse::add_option(pargs, c("--annotations_file"),
                              type="character",
                              action="store",
                              default=NULL,
                              dest="annotations_file",
                              metavar="Annotations_File",
                              help=paste("Path to file describing the cells, indicating the cell type classifications"))


args_parsed <- optparse::parse_args(pargs, positional_arguments=2)
args <- args_parsed$options
args["input_matrix"] <- args_parsed$args[1]
args["gene_order"] <- args_parsed$args[2]

# Check arguments
args <- check_arguments(args)


# Make sure the output directory exists
if(!file.exists(args$output_dir)){
    dir.create(args$output_dir)
}

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
    logging::loginfo(paste(":Input_Argument:",arg_name,"=",args[[arg_name]],
                           sep=""))
}
logging::loginfo(paste("::Input arguments. End."))

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=args$input_matrix,
                                    annotations_file=args$annotations_file,
                                    delim=args$delim,
                                    gene_order_file=args$gene_order,
                                    ref_group_names=args$name_ref_groups)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=args$output_dir,
                             cluster_by_groups=T,
                             plot_steps=F,
                             mask_nonDE_genes=T,
                             include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
                             )