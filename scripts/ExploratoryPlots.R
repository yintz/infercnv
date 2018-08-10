#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
pargs <- optparse::OptionParser()
pargs <- optparse::add_option(pargs, c("--file_dir"),
                              type="character",
                              action="store",
                              dest="file_dir",
                              metavar="File_Directory",
                              help=paste("Path to the expression data files created by inferCNV.",
                                         "[Default %default][REQUIRED]"))
args_parsed <- optparse::parse_args(pargs)
    # ----------------------Select Gene & Add Cell Types-----------------------------------------------------------------------
create_gene_eset <- function(ploting_data,
                             files_dir,
                             reference_idx, 
                             ref_groups) {
    cell_ids <- colnames(ploting_data)
    # label observed or Reference(s)
    ## find what index values are not in ref_groups 
    cell_type <- replace(cell_ids, !(1:length(cell_ids) %in% unlist(ref_groups)), paste("Observed"))
    #check
    if (any(!(reference_idx %in% cell_type))){
        stop("reference samples not in cell_type")
    }
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

    filtered_plot_data_t <- t(ploting_data )
    # add the labled cell annotation information
    final_eset <- data.table(filtered_plot_data_t, cell_type = cell_type, check.names = FALSE, stringsAsFactors = FALSE)
    return(final_eset)
}
    
CreatePlots <- function(files_dir){
    # ----------------------Load Environment and variables----------------------------------------------------------------------------------------
    inferCNVenv <- c("infercnv.preprocess.Rdata", "infercnv.processed.ward.D.Rdata")
    env <- paste(files_dir,inferCNVenv, sep = .Platform$file.sep)
    if (all(!(file.exists(env)))){
        path_error <- c(inferCNVenv)[which(!file.exists(inferCNVenv))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the enviorments created by inferCNV are in the following pathway: ")
        stop(error_message, paste(path_error, collapse = ", "))
    } else {
        for (i in env) {
            load(file = i)
        }
    }
    # ----------------------Load Data----------------------------------------------------------------------------------------
    filenames <- list.files(files_dir, pattern="*.pdf.txt", full.names=TRUE)
    reference_idx = ret_list$REF_OBS_IDX
    ref_groups = ret_list$REF_GROUPS
    if (file.exists(paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))) {
        message("Loading Rdata file.")
        load(paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))
    } else {
        message("Reading data files.")
        plot_data <- lapply(filenames , function(x) { 
            read.table(file = x, header=TRUE, row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
        })
        save(plot_data, file = paste(files_dir,"MinimalPathData.Rdata", sep = .Platform$file.sep))
    }
    data_list <- list()
    for (q in 1:length(plot_data)){
        q<-q
        data_list[[q]] <- create_gene_eset(ploting_data = plot_data[[q]], 
                                           files_dir = path,
                                           reference_idx = reference_idx, 
                                           ref_groups = ref_groups)
    }
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
    #### get chromosome seperator locations
    chrIdx <- cumsum(table(ret_list$CONTIGS))
    dat <- data_list[[9]]
    geneOrder <- head(colnames(dat),n=-1)
    
    # ----------------------plot 2 differnet cells  -----------------------------------------------------------------------
    message("Plotting a single reference cell against a observed cell.")
    pdf(paste(out_directory,"SINGLE.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    obs_dat <-  dat[cell_type == "Observed"]
    ref_dat <- dat[cell_type != "Observed"]
    twoCells <- rbind(obs_dat[1,],ref_dat[1,])
    singleData <- melt(twoCells, id = "cell_type")
    colnames(singleData) <- c("cell_type", "Gene", "Expression")
    singlePlot <- ggplot(singleData, aes(x = Gene, y = Expression, group = cell_type)) +
                    geom_line(aes(color = cell_type)) +
                    labs(title="Expression Levels Between One Reference And Observed Cell") + 
                    theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                          panel.background = element_blank(),
                          panel.grid.major.x = element_line(colour = "grey"),
                          legend.position = "bottom",
                          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
                    scale_x_discrete(limits = geneOrder,
                                     breaks = geneOrder[chrIdx],
                                     labels = c(names(chrIdx))) 
    print(singlePlot)
    dev.off()
    
    # ----------------------plot all values  -----------------------------------------------------------------------
    message("Plotting all expression values.")
    pdf(paste(out_directory,"ALL_VALUES.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    allData <- melt(dat, id = "cell_type")
    colnames(allData) <- c("cell_type", "Gene", "Expression")
    allPlot <- ggplot(allData, aes(x = Gene, y = Expression, group = cell_type)) +
                geom_line(aes(color = cell_type)) +
                labs(title="Expression Levels For All Genes") + 
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "bottom",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = c(names(chrIdx))) 
    print(allPlot)
    dev.off()
    
    # ----------------------PLOT MEAN VALUES-----------------------------------------------------------------------
    message("Plotting mean expression values.")
    pdf(paste(out_directory,"MEAN.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    geneMeans <- dat[,lapply(.SD,mean), by = cell_type]
    meanData <- melt(geneMeans, id = "cell_type")
    colnames(meanData) <- c("cell_type", "Gene", "Expression")
    meanPlot <- ggplot(meanData, aes(x = Gene, y = Expression, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Mean Expression") + 
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = c(names(chrIdx))) 
    print(meanPlot)
    dev.off()
    # ----------------------MEDIAN VALUES -----------------------------------------------------------------------
    message("Plotting median expression values.")
    pdf(paste(out_directory,"MEDIAN.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    geneMedians <- dat[,lapply(.SD,median), by = cell_type]
    medianData <- melt(geneMedians, id = "cell_type")
    colnames(medianData) <- c("cell_type", "Gene", "Expression")
    medianPlot <- ggplot(medianData, aes(x = Gene, y = Expression, group = cell_type)) +
                geom_line(aes(color = cell_type)) +
                labs(title="Median Expression") +
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "bottom",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = c(names(chrIdx))) 
    print(medianPlot)
    dev.off()
    
    # ----------------------VAR VALUES -----------------------------------------------------------------------
    message("Plotting Variance in expression within genes.")
    pdf(paste(out_directory,"Variance.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    geneVar <- dat[,lapply(.SD,var), by = cell_type]
    varData <- melt(geneVar, id = "cell_type")
    colnames(varData) <- c("cell_type", "Gene", "Variance")
    varPlot <- ggplot(varData, aes(x = Gene, y = Variance, group = cell_type)) +
                geom_line(aes(color = cell_type)) +
                labs(title="Variance In Expression Within Genes") +
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "bottom",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = c(names(chrIdx))) 
    print(varPlot)
    dev.off()
    
    # ----------------------MAD VALUES -----------------------------------------------------------------------
    message("Plotting MAD of gene expression.")
    pdf(paste(out_directory,"MAD.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    geneMAD <- dat[,lapply(.SD,mad), by = cell_type]
    MADData <- melt(geneMAD, id = "cell_type")
    colnames(MADData) <- c("cell_type", "Gene", "MAD")
    madPlot <- ggplot(MADData, aes(x = Gene, y = MAD, group = cell_type)) +
                geom_line(aes(color = cell_type)) +
                labs(title="MAD Expression") +
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "bottom",
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = c(names(chrIdx))) 
    print(madPlot)
    dev.off()
    
    
    # ----------------------KS.test between genes -----------------------------------------------------------------------
    message("Plotting KS and P-values.")
    # seperate the observed and reference data 
    obs_dat <- dat[cell_type == "Observed"]
    ref_dat <- dat[cell_type != "Observed"]
    obs_dat <- as.data.frame(obs_dat)[,1:(length(obs_dat[1,])-1)]
    ref_dat <- as.data.frame(ref_dat)[,1:(length(ref_dat[1,])-1)] 
    
    ksValues <- lapply(geneOrder,function(i,d1,d2){
        ks_values <- ks.test(d1[,i],d2[,i])
        c(KS = ks_values$statistic, p.value = ks_values$p.value)}, obs_dat, ref_dat)
    names(ksValues) <- geneOrder
    ksData <- t(data.frame(ksValues, check.names = FALSE) )
    ksData <- transform(ksData, KS.D = as.numeric(KS.D),
                         p.value = as.numeric(p.value))
    ksData$Gene <- row.names(ksData)
    #plot KS values
    ksPlot <- ggplot(data = ksData, aes(x = Gene, group = 1)) +
        geom_line(aes(y = KS.D), color = "#FF9999") +
        ylab(label="KS Value") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = c(names(chrIdx))) 
    #plot p-values
    pValPlot <- ggplot(data = ksData, aes(x = Gene, group = 1)) +
        geom_line(aes(y = p.value), color = "skyblue2") +
        ylab(label="P-Value") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = c(names(chrIdx))) 
    
    
    arrange_plots <- gridExtra::arrangeGrob( ksPlot, pValPlot,
                                             ncol=1, # use apply and get to get variables for each step instead
                                             top    = grid::textGrob("KS Values With Associated P-Values",
                                                                     gp = grid::gpar(fontsize=20,font=2)),
                                             bottom = grid::textGrob("Gene",
                                                                     gp = grid::gpar(fontsize=15)))
    # create the pdf and plot it on a grid 
    pdf(paste(out_directory,"KS.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    gridExtra::grid.arrange(arrange_plots, nrow = 2, ncol = 1, heights = c(10,1))
    dev.off()
}

if (!is.null(args)) {
    CreatePlots(files_dir = args_parsed$file_dir)
}
