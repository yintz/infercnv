#!/usr/bin/env Rscript


# Load libraries
library(ape)
library("RColorBrewer", character.only=TRUE)
library(GMD)
library(optparse)
library(logging)
library(infercnv)

# Logging level choices
C_LEVEL_CHOICES <- names(loglevels)
# Visualization outlier thresholding and bounding method choices
C_VIS_OUTLIER_CHOICES <- c("average_bound")
C_REF_SUBTRACT_METHODS <- c("by_mean", "by_quantiles")

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
                                C_REF_SUBTRACT_METHODS, collapse=",", sep="") )
        stop("error, must specify acceptable --ref_subtract_method")
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

        stop(980)
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
                              help=paste("Indicies of groups in which to group",
                                         "the reference samples; or a number",
                                         "of groups to make automatically",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--ref_subtract_method"),
                              type="character",
                              default="by_mean",
                              action="store",
                              dest="ref_subtract_method",
                              metavar="Reference_Subtraction_Method",
                              help=paste("Method used to subtract the reference values from the observations. Valid choices are",
                                         paste(C_REF_SUBTRACT_METHODS, collapse=", "),
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

# Manage inputs
logging::loginfo(paste("::Reading data matrix.", sep=""))
# Row = Genes/Features, Col = Cells/Observations
expression_data <- read.table(args$input_matrix, sep=args$delim, header=T, row.names=1)
logging::loginfo(paste("Original matrix dimensions (r,c)=",
                 paste(dim(expression_data), collapse=",")))

# Read in the gen_pos file
input_gene_order <- seq(1, nrow(expression_data), 1)
if (args$gene_order != ""){
    input_gene_order <- read.table(args$gene_order, header=F, row.names=1, sep="\t")
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
    # Normalize names with make.names so they are treated
    # as the matrix column names
    refs <- make.names(refs)
    if (length(refs) > 0){
        input_reference_samples <- refs
    }
    logging::logdebug(paste("::Reference observations set to: ", input_reference_samples, collapse="\n"))
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
    logging::logdebug(paste("::colnames(expression_data): ", colnames(expression_data), collapse="\n"))
    logging::logerror(error_message)
    stop(error_message)
}

# Order and reduce the expression to the genomic file.
order_ret <- infercnv::order_reduce(data=expression_data,
                                    genomic_position=input_gene_order)
expression_data <- order_ret$expr
input_gene_order <- order_ret$order
if(is.null(expression_data)){
    error_message <- paste("None of the genes in the expression data",
                           "matched the genes in the reference genomic",
                           "position file. Analysis Stopped.")
    stop(error_message)
}


if (args$save) {
    logging::loginfo("Saving workspace")
    save.image("infercnv.Rdata")
}

# Run CNV inference
ret_list = infercnv::infer_cnv(data=expression_data,
                               gene_order=input_gene_order,
                               cutoff=args$cutoff,
                               reference_obs=input_reference_samples,
                               transform_data=args$log_transform,
                               window_length=args$window_length,
                               max_centered_threshold=args$max_centered_expression,
                               noise_threshold=args$magnitude_filter,
                               num_ref_groups=args$num_groups,
                               out_path=args$output_dir,
                               k_obs_groups=args$num_obs,
                               plot_steps=args$plot_steps,
                               contig_tail=args$contig_tail,
                               method_bound_vis=args$bound_method_vis,
                               lower_bound_vis=bounds_viz[1],
                               upper_bound_vis=bounds_viz[2],
                               ref_subtract_method=args$ref_subtract_method)

# Log output
logging::loginfo(paste("::infer_cnv:Writing final data to ",
                       file.path(args$output_dir,
                       "expression_pre_vis_transform.txt"), sep="_"))
# Output data before viz outlier
write.table(ret_list["PREVIZ"], sep=args$delim,
            file=file.path(args$output_dir,
                       "expression_pre_vis_transform.txt"))
# Output data after viz outlier
write.table(ret_list["VIZ"], sep=args$delim,
            file=file.path(args$output_dir,
                       "expression_post_viz_transform.txt"))
logging::loginfo(paste("::infer_cnv:Current data dimensions (r,c)=",
                       paste(dim(ret_list[["VIZ"]]), collapse=","), sep=""))

logging::loginfo(paste("::infer_cnv:Drawing plots to file:",
                           args$output_dir, sep=""))


if (args$save) {
    logging::loginfo("Saving workspace")
    save.image("infercnv.Rdata")
}


if (args$plot_steps) {
    logging::loginfo("See results from each stage plotted separately")
}  else {
      
    infercnv::plot_cnv(plot_data=ret_list[["VIZ"]],
                       contigs=ret_list[["CONTIGS"]],
                       k_obs_groups=args$num_obs,
                       reference_idx=ret_list[["REF_OBS_IDX"]],
                       ref_contig=args$clustering_contig,
                       contig_cex=args$contig_label_size,
                       ref_groups=ret_list[["REF_GROUPS"]],
                       out_dir=args$output_dir,
                       color_safe_pal=args$use_color_safe,
                       title=args$fig_main,
                       obs_title=args$obs_main,
                       ref_title=args$ref_main)
    
}

