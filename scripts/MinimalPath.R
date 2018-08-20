#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
# Command line arguments
pargs <- optparse::OptionParser()
pargs <- optparse::add_option(pargs, c("--file_dir"),
                              type="character",
                              action="store",
                              dest="file_dir",
                              metavar="File_Directory",
                              help=paste("Path to the expression data files created by inferCNV.",
                                         "[Default %default][REQUIRED]"))
pargs <- optparse::add_option(pargs, c("--genes"),
                              type="character",
                              action="store",
                              dest="genes",
                              metavar="Genes",
                              default = NULL,
                              help=paste("Desired Gene to plot.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--one_ref"),
                              type="logical",
                              action="store",
                              dest="one_ref",
                              default = FALSE,
                              metavar="ONE_REFERENCE",
                              help=paste("Option to combine mulitple references into one.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--contig"),
                              type="character",
                              action="store",
                              dest="contig",
                              default = NULL,
                              metavar="CHOOSE_CHROMOSOME",
                              help=paste("Chromosome in which the median density will be plotted.",
                                         "[Default %default]"))
args_parsed <- optparse::parse_args(pargs)


#' @title Plot Gene Distributions Along InferCNV Steps
#' @description Takes in expression data created by inferCNV and plots the 
#'              distribution of desired genes expression or entire chromosome 
#'              along each step in inferCNV. Plots are output as a PDF file. 
#' @param files_dir (string) Path to where the expression data files produced by inferCNV are located.
#' @param genes (string or vector) Gene(s) of interest.
#' @param one_ref (logical) Option to have references combined into one if TRUE or separated into different references if FALSE. Default value is FALSE.
#' @param contig (string) The label for the chromosome that will be plotted. 
#' @return 
#' PDF File with density plots showing gene expression levels along each step in inferCNV. 
#' @export
# Requires:
#   gridExtra, grid, data.table, ggplot2
MinimalPath <- function (files_dir,
                         genes = NULL,
                         one_ref = FALSE,
                         contig = NULL){
    
    # ----------------------Select Gene & Add Cell Types-----------------------------------------------------------------------
    create_gene_eset <- function(plot_data, 
                                 gene,
                                 reference_idx, 
                                 ref_groups) {
        # Function to create identify references 
        #		Input: 
        #				plot_data: (data.table) expression data, rows = Genes, columns = Cells.
        #				gene: (string) Gene of interest.
        #				reference_idx: (vector) Vector containing cell line reference identifiers.
        #				ref_groups: (list) List containing index locations of reference groups.
        #       Output:
        #               Data frame with expression data and annotation information
        
        cell_ids <- colnames(plot_data)
        # label observed or Reference(s)
        ## find what index values are not in ref_groups 
        cell_type <- replace(cell_ids, !(1:length(cell_ids) %in% unlist(ref_groups)), paste("Observed"))
        #check
        if (any(!(input_reference_samples %in% cell_type))){
            stop("reference samples not in cell_type")
        }
        ## Label the references based on index locations 
        if (one_ref == FALSE) {
            ## Label the references based on index locations 
            if (length(ref_groups) > 1) {
                for(i in 1:length(ref_groups)){ 
                    cell_type <- replace(cell_type, ref_groups[[i]], paste("Reference",toString(i),sep = "")) 
                }
            } else {
                for(i in 1:length(ref_groups)){ 
                    cell_type <- replace(cell_type, ref_groups[[1]], paste("Reference"))
                }
            }
        } else {
            cell_type <- replace(cell_type, unlist(ref_groups), paste("Reference"))
        }
        
        # check if all reference cells are in cell type 
        obs_cells <- which(cell_type != "Observed")
        if (!(all(reference_idx %in% cell_ids[obs_cells]))){
            missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
            error_message <- paste("Error: Not all references are accounted for.",
                                   "Make sure the reference names match the names in the data.\n",
                                   "Check the following reference cell lines: ", 
                                   paste(missing_refs, collapse = ","))
            stop(error_message)
        }
        #filter the expression set to only include the gene of interest 
        if (!missing(gene)) {
            new_plot_data <- plot_data[gene,] 
            filtered_plot_data_t <- t(new_plot_data)
            # add the labled cell annotation information
            final_eset <- data.frame(filtered_plot_data_t, cell_type,
                                     check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
            colnames(final_eset) <- c("exp","cell_type")
        } else {
            filtered_plot_data_t <- t(plot_data )
            # add the labled cell annotation information
            final_eset <- data.table(filtered_plot_data_t, cell_type = cell_type, check.names = FALSE, stringsAsFactors = FALSE)
        }
        return(final_eset)
    }
    
    # ----------------------Check Pathways----------------------------------------------------------------------------------------
    # Error handling 
    ## check if pathway exists
    env <- list.files(files_dir, pattern="infercnv.*.Rdata", full.names=TRUE)
    
    if (any(!c(file.exists(files_dir,env)))){
        path_error <- c(files_dir,env)[which(!file.exists(files_dir,env))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the following pathways are correct: ")
        stop(error_message, paste(path_error, collapse = ", "))
    }
    # ----------------------Load Environment and variables----------------------------------------------------------------------------------------
    if ( length(env) == 0 ) {
        stop("No enviornmental files found")
    }
    for (i in env) {
        load(file = i)
    }
    ref_groups = ret_list$REF_GROUPS
    reference_idx = ret_list$REF_OBS_IDX
    
    # ----------------------Load Data----------------------------------------------------------------------------------------
    filenames <- list.files(files_dir, pattern="*.pdf.txt", full.names=TRUE)
    if (file.exists(paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))){
        logging::loginfo("Loading Rdata file.")
        load(paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))
    } else {
        logging::loginfo("Reading data files.")
        plot_data <- lapply(filenames , function(x) { 
            read.table(file = x, sep=delim, header=TRUE, row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
        })
        save(plot_data, file = paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))
    }
    # ----------------------Create output directory----------------------------------------------------------------------------------------   
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    
    # ----------------------Create graph title----------------------------------------------------------------------------------------
    concat <- function(files){ 
        # Create the titles for each graph:
        #       Gets the names of the files and capitalizes first letter of each word 
        #       input: the paths to file 
        #       output: vector of graph titles 
        temp1 <- sapply(strsplit(files, "/"), tail, n = 1)
        temp2 <- sapply(strsplit(temp1, ".pdf.txt"),head,n=1)
        temp3 <- strsplit(temp2, "_")
        sapply(temp3, function(x){
            tmp <- paste0(toupper(substring(x[-1], 1,1)), substring(x[-1], 2), collapse=" ")
            paste(x[1], tmp)
        })
    }
    titles <- concat(filenames)
    
    if (!is.null(genes)){
        if (is.character(genes)) {
            genes <- unlist(strsplit(genes, split = ","))
        }
        # check to see if the gene is in the data 
        logging::loginfo(paste0("Checking the following genes exists in data: ", paste(genes, collapse = ", ")))
        missing_gene <- lapply(plot_data, function(x){
            if (!(all(genes %in% row.names(x)))){
                genes[which(!(genes %in% row.names(x)))]
            } 
        })
        if (!is.null(unlist(missing_gene))) {
            stop(paste0("Error: The following gene(s) are not found in expression data: ", paste(unique(unlist(missing_gene)), collapse = ", ")))
        }
        
        
        
        # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
        for (gene in genes){
            logging::loginfo(paste0("Plotting ",gene))
            # run create_gene_eset function on each file to gather the expression data 
            data_list <- list()
            for (q in 1:length(plot_data)){
                q<-q
                data_list[[q]] <- create_gene_eset(plot_data = plot_data[[q]], 
                                                   gene = gene, 
                                                   reference_idx = reference_idx, 
                                                   ref_groups = ref_groups)
            }
            # check that all cell id's match 
            if(length(unique(lapply(data_list, function(x) row.names(x))))!=1){
                stop("Cell id's don't match between plotting data.")
            }
            
            plot_list <- list()
            for (i in 1:length(data_list)) {
                plot_list[[i]] <-
                    assign(paste0("plot",i), 
                           ggplot(na.omit(data_list[[i]]), aes(x = exp, fill = cell_type, colour = cell_type)) +
                               geom_density(alpha = 0.3, adjust = 3/5) +
                               labs(title=titles[i]) +
                               theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 9),
                                     panel.background = element_blank(),
                                     legend.position = "none",
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank())+
                               scale_x_continuous(breaks = pretty( data_list[[i]]$exp, n = 5)))
            }
            # Generate output file name with the desired output directory 
            ## File name is the name of the gene 
            output_file_name <- paste0(gene, ".pdf")
            pdf(file.path(out_directory,output_file_name), onefile = FALSE, height = 9, width = 9)
            pull_legend<-function(plot){
                plot_params <- ggplot_gtable(ggplot_build(plot))
                # get which grob in list is the legend called "guide=box" 
                legend_idx <- which(sapply(plot_params$grobs, function(x) x$name) == "guide-box")
                legend <- plot_params$grobs[[legend_idx]]
                return(legend)}
            temp <- plot_list[[1]] + 
                theme(legend.position = "bottom") + 
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(temp)
            
            arrange_plots <- gridExtra::arrangeGrob( grobs = plot_list,
                                                     ncol=3, # use apply and get to get variables for each step instead
                                                     top    = grid::textGrob(gene,
                                                                             gp = grid::gpar(fontsize=20,font=2)),
                                                     left   = grid::textGrob("Density",
                                                                             gp = grid::gpar(fontsize=15), rot = 90),
                                                     bottom = grid::textGrob("Expression Levels",
                                                                             gp = grid::gpar(fontsize=15)))
            # create the pdf and plot it on a grid 
            gridExtra::grid.arrange(arrange_plots,
                                    plot_legend, nrow = 2, ncol = 1, heights = c(10,1))
            dev.off()
        }
    }
    if (!is.null(contig)){
        data_list <- list()
        for (q in 1:length(plot_data)){
            q<-q
            data_list[[q]] <- create_gene_eset(plot_data = plot_data[[q]],
                                               reference_idx = reference_idx, 
                                               ref_groups = ref_groups)
        }
        location <- data.table(input_gene_order, "Gene" = row.names(input_gene_order))
        plot_list <- list()
        for (i in 1:length(data_list)) {
            currentGenes <- intersect(colnames(data_list[[i]]),location[chr == contig,]$Gene)
            geneMedians <- as.data.table(data_list[[i]])[,c("cell_type",currentGenes), with = FALSE][,lapply(.SD,median), by = cell_type]
            df <- melt(geneMedians, id = "cell_type")
            colnames(df) <- c("cell_type", "Gene", "Median")
            plot_list[[i]] <- assign(paste0("plot",i), 
                                     ggplot(df, aes(x = Median, fill = cell_type, colour = cell_type)) + 
                                         geom_density(alpha = 0.3, adjust = 3/5) +
                                         labs(title=titles[i]) +
                                         theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 9),
                                               panel.background = element_blank(),
                                               legend.position = "none",
                                               axis.title.x = element_blank(),
                                               axis.title.y = element_blank())+
                                         scale_x_continuous(breaks = pretty(df$Median, n = 5)))
        }
        pull_legend<-function(plot){
            plot_params <- ggplot_gtable(ggplot_build(plot))
            # get which grob in list is the legend called "guide=box" 
            legend_idx <- which(sapply(plot_params$grobs, function(x) x$name) == "guide-box")
            legend <- plot_params$grobs[[legend_idx]]
            return(legend)
        }
        temp <- plot_list[[1]] + 
            theme(legend.position = "bottom") +
            guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
        plot_legend <- pull_legend(temp)
        
        output_file_name <- paste0(contig,".pdf")
        pdf(file.path(out_directory,output_file_name), onefile = FALSE, height = 9, width = 9)
        # create the pdf and plot it on a grid 
        grid.arrange(arrangeGrob( grobs = plot_list,
                                  ncol=3, # use apply and get to get variables for each step instead
                                  top    = textGrob(contig,
                                                    gp = gpar(fontsize=20,font=2)),
                                  left   = textGrob("Density",
                                                    gp = gpar(fontsize=15), rot = 90),
                                  bottom = textGrob("Expression Levels",
                                                    gp = gpar(fontsize=15))),
                     plot_legend, nrow = 2, ncol = 1, heights = c(10,1))
    }
}


if (!is.null(args)) {
    MinimalPath(files_dir = args_parsed$file_dir, 
                genes = args_parsed$genes,
                one_ref = args_parsed$one_ref,
                contig = args_parsed$contig)
}
