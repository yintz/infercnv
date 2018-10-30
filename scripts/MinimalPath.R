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
#' 
#' @return 
#' PDF File with density plots showing gene expression levels along each step in inferCNV. 
#' @export
# Requires:
#   gridExtra, grid, data.table, ggplot2


configure_data <- function(infercnv_obj,
                           gene,
                           cluster_by_groups,
                           one_ref,
                           contig){
    # Function to identify references and observed cell lines within expresion dataframe 
    #		Input: 
    #				infercnv_obj: (S4) infercnv object 
    #				gene: (string) Gene of interest.
    #               one_ref: (boolean) If True, will combine references into one 
    #       Output:
    #               Data frame with expression data and annotation information
    
    # get the cell ID names in order, (reference, observed)
    row_order <- c(colnames(infercnv_obj@expr.data))
    cell_type <- replace(row_order, 1:length(row_order) %in% unlist(infercnv_obj@observation_grouped_cell_indices), paste("Observed"))
    
    ## Label the references based on index locations 
    if (one_ref == FALSE || length(infercnv_obj@reference_grouped_cell_indices) > 1) {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, infercnv_obj@reference_grouped_cell_indices[[i]], paste("Reference",toString(i),sep = "")) 
        }
    } else {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, 1:length(cell_type) %in% unlist(infercnv_obj@reference_grouped_cell_indices), paste("Reference"))
        }
    }
    
    #filter the expression set to only include the gene of interest 
    if (!missing(gene)) {
        new_plot_data <- infercnv_obj@expr.data[gene,] 
        # add the labled cell annotation information
        final_data <- data.frame(new_plot_data, cell_type,
                                 check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
        colnames(final_data) <- c("exp","cell_type")
    } 
    if (!missing(contig)){
        # check if contig in data
        if (contig %in% infercnv_obj@gene_order[[1]]){
            final_data <- data.frame(t(infercnv_obj@expr.data), cell_type,
                                     check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
        } else {
            error_message <- paste("Contig selected is not found in the data.\n",
                                   "Check to make sure the correct contig label was selected.")
            stop(error_message)
        }
    }
    
    # check if all reference cells are in cell type 
    ref_cells <- which(cell_type != "Observed")
    reference_ids <- colnames(infercnv_obj@expr.data)[unlist(infercnv_obj@reference_grouped_cell_indices)]
    if (!(all(reference_ids %in% colnames(infercnv_obj@expr.data)[ref_cells]))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
        error_message <- paste("Error: Not all references are accounted for.",
                               "Make sure the reference names match the names in the data.\n",
                               "Check the following reference cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    # check if all observed cells are in cell type 
    observed_ids <-colnames(infercnv_obj@expr.data)[unlist(infercnv_obj@observation_grouped_cell_indices)]
    obs_cells <- which(cell_type == "Observed")
    if (!(all(observed_ids %in% colnames(infercnv_obj@expr.data)[obs_cells]))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
        error_message <- paste("Error: Not all observed cell lines are accounted for.",
                               "Make sure the cell id names match the names in the data.\n",
                               "Check the following oberved cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    return(final_data)
}

concat <- function(files){ 
    # Function to Create the titles for each graph:
    #   Gets the names of the files and capitalizes first letter of each word 
    #       input: the paths to file 
    #       output: vector of graph titles 
    sub_titles <- sapply(files, function(x){
        unlist(strsplit(basename(x),split="\\."))[1]
    })
    temp <- strsplit(sub_titles, "_")
    sapply(temp, function(x){
        tmp <- paste0(toupper(substring(x[-1], 1,1)), substring(x[-1], 2), collapse=" ")
        paste(x[1], tmp)
    })
}

pull_legend<-function(plot){
    plot_params <- ggplot_gtable(ggplot_build(plot))
    # get which grob in list is the legend called "guide=box" 
    legend_idx <- which(sapply(plot_params$grobs, function(x) x$name) == "guide-box")
    legend <- plot_params$grobs[[legend_idx]]
    return(legend)
}


MinimalPath <- function (files_dir,
                         genes = NULL,
                         one_ref = FALSE,
                         contig = NULL){
    
    # ----------------------Check Pathways----------------------------------------------------------------------------------------
    # Error handling 
    ## check if pathway exists
    obj_files <- list.files(files_dir, pattern="*infercnv_obj", full.names=TRUE)
    
    if (any(!c(file.exists(files_dir, obj_files)))){
        path_error <- c(files_dir,obj_files)[which(!file.exists(files_dir,obj_files))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the following pathways are correct: ")
        stop(error_message, paste(path_error, collapse = ", "))
    }
    # ----------------------Load infercnv objects----------------------------------------------------------------------------------------
    if ( length(obj_files) == 0 ) {
        stop("No enviornmental files found")
    }
    for (i in obj_files) {
        load(file = i)
    }

    obj_id <- lapply(obj_files, function(x){
        load(x)})
    final_obj <- get(tail(obj_id,n=1)[[1]])
    
    ref_groups = names(final_obj@reference_grouped_cell_indices)
    reference_idx = unlist(final_obj@reference_grouped_cell_indices)
    
    # ----------------------Create output directory----------------------------------------------------------------------------------------   
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    
    # ----------------------Create graph title----------------------------------------------------------------------------------------
    titles <- concat(obj_files)
    
    if (!is.null(genes)){
        if (is.character(genes)) {
            genes <- unlist(strsplit(genes, split = ","))
        }
        # check to see if the gene is in the data 
        logging::loginfo(paste0("Checking the following genes exists in data: ", paste(genes, collapse = ", ")))
        missing_gene <- lapply(obj_id, function(x){
            if (!(all(genes %in% row.names(get(x)@expr.data)))){
                genes[which(!(genes %in% row.names(get(x)@expr.data)))]
            }
        })
        if (!is.null(unlist(missing_gene))) {
            stop(paste0("Error: The following gene(s) are not found in expression data: ", paste(unique(unlist(missing_gene)), collapse = ", ")))
        }
        
        
        
        # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
        for (gene in genes){
            logging::loginfo(paste0("Plotting ",gene))
            # run configure_data function on each file to gather the expression data 
            data_list <- list()
            for (q in 1:length(obj_id)){
                q<-q
                data_list[[q]] <- configure_data(infercnv_obj = get(obj_id[[q]]), 
                                                   gene = gene,
                                                   one_ref = one_ref)
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
            temp <- plot_list[[1]] + 
                theme(legend.position = "bottom") + 
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(temp)
            
            arrange_plots <- gridExtra::arrangeGrob( grobs = plot_list,
                                                     ncol=3, # use apply and get to get variables for each step instead
                                                     top    = grid::textGrob(gene,
                                                                             gp = grid::gpar(fontsize=20,font=2)),
                                                     left   = grid::textGrob("Density",
                                                                             gp = grid::gpar(fontsize=15), rot=90),
                                                     bottom = grid::textGrob("Expression Levels",
                                                                             gp = grid::gpar(fontsize=15)))
            # create the pdf and plot it on a grid 
            gridExtra::grid.arrange(arrange_plots, plot_legend, 
                                    nrow = 2, ncol = 1, 
                                    heights = c(10,1))
            dev.off()
        }
    }
    
    if (!is.null(contig)){
        data_list <- list()
        for (q in 1:length(obj_id)){
            q<-q
            data_list[[q]] <- configure_data(infercnv_obj = get(obj_id[[q]]),
                                             one_ref = one_ref,
                                             contig = contig)
        }
        
        plot_list <- list()
        for (i in 1:length(data_list)) {
            current_order <- get(obj_id[[i]])@gene_order
            current_genes <- row.names(current_order[which(current_order[1] == contig),])
            
            geneMedians <- as.data.table(data_list[[i]])[,c("cell_type",current_genes), with = FALSE][,lapply(.SD,median), by = cell_type]
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
