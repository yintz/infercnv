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
pargs <- optparse::add_option(pargs, c("--cluster_by_groups"),
                              type="logical",
                              action="store",
                              dest="cluster_by_groups",
                              default = "TRUE",
                              metavar="Cluster_By_Groups",
                              help=paste("Option to cluster the observed cells by unique type.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--one_reference"),
                              type="logical",
                              action="store",
                              dest="one_reference",
                              default = "FALSE",
                              metavar="One_Reference",
                              help=paste("Option to treat all references as one.",
                                         "[Default %default]"))
args_parsed <- optparse::parse_args(pargs)

pull_legend<-function(plot){
    plot_params <- ggplot_gtable(ggplot_build(plot))
    # get which grob in list is the legend called "guide=box" 
    legend_idx <- which(sapply(plot_params$grobs, function(x) x$name) == "guide-box")
    legend <- plot_params$grobs[[legend_idx]]
    return(legend)
}
# ----------------------Select Gene & Add Cell Types-----------------------------------------------------------------------
create_gene_eset <- function(plotting_data,
                             reference_idx, 
                             ref_groups,
                             cluster_by_groups,
                             obs_annotations_groups,
                             obs_annotations_names,
                             one_ref
) {
    cell_ids <- colnames(plotting_data)
    cell_type <- colnames(plotting_data)
    # label observed or Reference(s)
    ## find what index values are not in ref_groups 
    obs_idx <- c(1:length(cell_ids))[-unlist(ref_groups)]
    if (cluster_by_groups == TRUE){
        cell_type[obs_idx] <- sapply(obs_annotations_groups, function(x){ obs_annotations_names[x] })
    } else {
        cell_type <- replace(cell_ids, !(1:length(cell_ids) %in% unlist(ref_groups)), paste("Observed"))
    }
    
    
    #check
    if (any(!(reference_idx %in% cell_type))){
        stop("reference samples not in cell_type")
    }
    ## Label the references based on index locations 
    if (length(ref_groups) > 1 & one_ref == FALSE) {
        for(i in 1:length(ref_groups)){ 
            cell_type <- replace(cell_type, ref_groups[[i]], paste("Reference",toString(i),sep = "")) 
        }
    } else {
        for(i in 1:length(ref_groups)){ 
            cell_type <- replace(cell_type, unlist(ref_groups), paste("Reference"))
        }
    }
    # check if all reference cells are in cell type 
    ref_cells <- which(!(cell_type %in% c(obs_annotations_names)))
    if (!(all(reference_idx %in% cell_ids[ref_cells]))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
        error_message <- paste("Error: Not all references are accounted for.",
                               "Make sure the reference names match the names in the data.\n",
                               "Check the following reference cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    
    filtered_plot_data_t <- t(plotting_data )
    # add the labled cell annotation information
    final_eset <- data.table(filtered_plot_data_t, cell_type = cell_type, check.names = FALSE, stringsAsFactors = FALSE)
    return(final_eset)
}

CreatePlots <- function(files_dir,
                        cluster_by_groups = TRUE,
                        one_reference = FALSE){
    # ----------------------Load Environment and variables----------------------------------------------------------------------------------------
    print(files_dir)
    env <- list.files(files_dir, pattern="infercnv.*.Rdata", full.names=TRUE)
    if (all(!(file.exists(env)))){
        path_error <- c(env)[which(!file.exists(env))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the enviorments created by inferCNV are in the following pathway: ")
        stop(error_message, paste(path_error, collapse = ", "))
    } else {
        for (i in env) {
            load(file = i)
        }
    }
    # ----------------------Load Data----------------------------------------------------------------------------------------
    reference_idx <- ret_list$REF_OBS_IDX
    ref_groups <- ret_list$REF_GROUPS
    filenames <- list.files(files_dir, pattern="*.pdf.txt", full.names=TRUE)
    
    if (length(filenames)>0){
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
            obs_idx <- c(1:ncol(plot_data[[q]]))[-unlist(ref_groups)]
            data_list[[q]] <- create_gene_eset(plotting_data = plot_data[[q]], 
                                               reference_idx = reference_idx, 
                                               ref_groups = ref_groups,
                                               cluster_by_groups = cluster_by_groups,
                                               obs_annotations_groups = obs_annotations_groups,
                                               obs_annotations_names = obs_annotations_names,
                                               one_ref = one_reference
            )
        }
        dat <- data_list[[9]]
    } else {
        message("loading data from environmnet")
        obs_idx <- c(1:ncol(data.frame(ret_list$VIZ)))[-unlist(ref_groups)]
        data_list <- create_gene_eset(plotting_data = data.frame(ret_list$VIZ),
                                      reference_idx = reference_idx, 
                                      ref_groups = ref_groups,
                                      cluster_by_groups = cluster_by_groups,
                                      obs_annotations_groups = obs_annotations_groups,
                                      obs_annotations_names = obs_annotations_names,
                                      one_ref = one_reference
        )
        dat <- data_list
    }
    # create directory
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
    #### get chromosome seperator locations
    contigOrder <- levels(input_gene_order$chr)
    counts <- table(factor(ret_list$CONTIGS),useNA = "no")[contigOrder]
    cumulativeCounts <- cumsum(c(1,counts[-which(is.na(names(counts)))]))
    chrIdx <- head(cumulativeCounts,n=-1L)
    geneOrder <- head(colnames(dat),n=-1)
    # create PDF 
    pdf(paste(out_directory,"ExploratoryPlots.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 6, width = 9)
    # ----------------------plot 2 differnet cells  -----------------------------------------------------------------------
    message("Plotting a single reference cell against a observed cell.")
    obs_dat <-  dat[obs_idx]
    ref_dat <- dat[-obs_idx]
    twoCells <- rbind(obs_dat[1,],ref_dat[1,])
    singleData <- melt(twoCells, id = "cell_type")
    colnames(singleData) <- c("cell_type", "Gene", "Expression")
    singlePlot <- ggplot(singleData, aes(x = Gene, y = Expression, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Expression Of One Reference And One Observed Cell") + 
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 10),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(singlePlot)
    
    
    # ----------------------plot all values  -----------------------------------------------------------------------
    # message("Plotting all expression values.")
    # #pdf(paste(out_directory,"ALL_VALUES.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 9, width = 9)
    # allData <- melt(dat, id = "cell_type")
    # colnames(allData) <- c("cell_type", "Gene", "Expression")
    # allPlot <- ggplot(allData, aes(x = Gene, y = Expression, group = cell_type)) +
    #     geom_line(aes(color = cell_type)) +
    #     labs(title="Expression Levels For All Genes") + 
    #     theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
    #           panel.background = element_blank(),
    #           panel.grid.major.x = element_line(colour = "grey"),
    #           legend.position = "bottom",
    #           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    #     scale_x_discrete(limits = geneOrder,
    #                      breaks = geneOrder[chrIdx],
    #                      labels = names(cumulativeCounts[-1]))
    #print(allPlot)
    #dev.off()
    # ----------------------PLOT MAX VALUES-----------------------------------------------------------------------
    # max_val <- dat[,lapply(.SD,max), by = cell_type]
    # maxData <- melt(max_val, id = "cell_type")
    # colnames(maxData) <- c("cell_type","Gene","Max_Value")
    # ggplot(maxData, aes(x = Gene, y = Max_Value, group = cell_type)) +
    #     geom_line(aes(color = cell_type)) +
    #     labs(title="Max Expression") + 
    #     theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
    #           panel.background = element_blank(),
    #           panel.grid.major.x = element_line(colour = "grey"),
    #           legend.position = "bottom",
    #           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    #     scale_x_discrete(limits = geneOrder,
    #                      breaks = geneOrder[chrIdx],
    #                      labels = names(cumulativeCounts[-1]))
    
    # ----------------------PLOT MEAN VALUES-----------------------------------------------------------------------
    message("Plotting mean expression values.")
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
                         labels = names(cumulativeCounts[-1]))
    print(meanPlot)
    
    # ----------------------MEDIAN VALUES -----------------------------------------------------------------------
    message("Plotting median expression values.")
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
                         labels = names(cumulativeCounts[-1]))
    print(medianPlot)
    
    # ----------------------VAR VALUES -----------------------------------------------------------------------
    message("Plotting Variance in expression within genes.")
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
                         labels = names(cumulativeCounts[-1]))
    print(varPlot)
    
    # ----------------------MAD VALUES -----------------------------------------------------------------------
    message("Plotting MAD of gene expression.")
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
                         labels = names(cumulativeCounts[-1]))
    print(madPlot)
    
    
    # ----------------------KS.test between genes -----------------------------------------------------------------------
    message("Plotting KS and P-values.")
    # seperate the observed and reference data 
    obs_dat <-  dat[obs_idx]
    ref_dat <- dat[-obs_idx]
    obs_dat <- as.data.frame(obs_dat)[,1:(ncol(obs_dat)-1)]
    ref_dat <- as.data.frame(ref_dat)[,1:(ncol(ref_dat)-1)] 
    
    ksValues <- lapply(geneOrder,function(i,d1,d2){
        ks_values <- ks.test(d1[,i],d2[,i])
        c(KS = ks_values$statistic, p.value = (-10*log10(ks_values$p.value)))}, obs_dat, ref_dat)
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
                         labels = names(cumulativeCounts[-1]))
    #plot p-values
    pValPlot <- ggplot(data = ksData, aes(x = Gene, group = 1)) +
        geom_line(aes(y = p.value), color = "skyblue2") +
        ylab(label="Phred Score (-10log10(P-Value))") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    
    # create the pdf and plot it on a grid 
    KSPLOT <- gridExtra::arrangeGrob( ksPlot, pValPlot,
                                      ncol=1, # use apply and get to get variables for each step instead
                                      top    = grid::textGrob("KS Values With Associated P-Values",
                                                              gp = grid::gpar(fontsize=20,font=2)),
                                      bottom = grid::textGrob("Gene",
                                                              gp = grid::gpar(fontsize=15)))
    gridExtra::grid.arrange(KSPLOT, ncol = 1, heights = c(6))
    
    # ----------------------t.test between genes -----------------------------------------------------------------------
    message("Plotting t-test.")
    # seperate the observed and reference data 
    plotTtest <- function(dat){
        obs_dat <- dat[obs_idx]
        ref_dat <- dat[-obs_idx]
        obs_dat <- as.data.frame(obs_dat)[,1:( ncol(obs_dat)-1)]
        ref_dat <- as.data.frame(ref_dat)[,1:(ncol(ref_dat)-1)] 
        
        ttestValues <- data.frame(t(sapply(geneOrder,function(i,d1,d2){
            values <- (t.test(d1[,i],d2[,i]))
            ttest <- c(t.score = values$statistic, # t value the calculated difference represented in units of standard error (sigma)
                       t.pvalue = values$p.value,
                       confidence = values$conf.int) 
        }, obs_dat, ref_dat)))
        ttestValues$Gene <- row.names(ttestValues)
        ttestValues$t.fdr <- p.adjust(ttestValues$t.pvalue, method = "BH")
        ttestValues$t.pvalue <- -10*log10(ttestValues$t.pvalue)
        ttestValues$t.fdr <- -10*log10(ttestValues$t.fdr)
        
        ttest_score <- ggplot(data = ttestValues, aes(x = Gene)) +
            geom_line(aes(y = t.score.t, color = "T Value", group = 1)) +
            #geom_line(aes(y = confidence1, color = "blue", group = 1)) +
            labs(title="T-Value") +
            theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                  panel.background = element_blank(),
                  panel.grid.major.x = element_line(colour = "grey"),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
            scale_x_discrete(limits = geneOrder,
                             breaks = geneOrder[chrIdx],
                             labels = names(cumulativeCounts[-1]))
        ttest_pval <- ggplot(data = ttestValues, aes(x = Gene)) +
            geom_line(aes(y = t.pvalue, color = "P-Value", group = 1)) +
            geom_line(aes(y = t.fdr, color = "Adjust p-Value", group = 1)) +
            labs(title="P-Value And Adjusted P-Value") +
            theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                  panel.background = element_blank(),
                  panel.grid.major.x = element_line(colour = "grey"),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
            scale_x_discrete(limits = geneOrder,
                             breaks = geneOrder[chrIdx],
                             labels = names(cumulativeCounts[-1]))
        TTESTPLOT <- gridExtra::arrangeGrob( ttest_score, ttest_pval,
                                             ncol=1, # use apply and get to get variables for each step instead
                                             top    = grid::textGrob("t-Test Values With Associated P-Values",
                                                                     gp = grid::gpar(fontsize=20,font=2)))
        gridExtra::grid.arrange(TTESTPLOT, ncol = 1, heights = c(6))
        
    }
    print(plotTtest(dat))
    
    #adjacent covarience 
    # adjCov <- function(dat){
    #     n = table(dat$cell_type)[unique(dat$cell_type)]
    #     covData <- list()
    #     for (i in 1:length(n)){
    #         print(names(n[i]))
    #         cell <- as.data.frame(dat[cell_type == names(n[i]),])
    #         covData[[i]] <- data.frame(apply(cell,2, function(j){
    #             cov(as.numeric(j[1:(n[i]-1)]), as.numeric(j[2:n[i]]))
    #             }))
    #         covData[[i]]$Gene <- row.names(covData[[i]])
    #         covData[[i]]$cell_type <- names(n[i])
    #         colnames(covData[[i]]) <- c("Covariation","Gene", "cell_type")
    #     }
    #     for (i in 1:length(covData)){
    #         covData[[1]] <- rbind(covData[[1]],covData[[i]])
    #     }
    # covPlot <- ggplot(data = covData[[1]], aes(x = Gene, y = Covariation, colour = cell_type, fill = cell_type)) +
    #     geom_point(alpha = .3) +
    #     scale_x_discrete(limits = geneOrder,
    #                      breaks = geneOrder[chrIdx],
    #                      labels = names(cumulativeCounts[-1]))+
    #     ylim(-.1,.1)
    # return(covPlot)
    # }
    # 
    # TTESTPLOT <- gridExtra::arrangeGrob( 
    #                                      adjCov(data_list[[7]]),
    #                                      adjCov(data_list[[8]]),
    #                                      adjCov(data_list[[9]]),
    #                                      ncol=1, # use apply and get to get variables for each step instead
    #                                      top    = grid::textGrob("Adjacent Covariation",
    #                                                              gp = grid::gpar(fontsize=20,font=2)),
    #                                      bottom = grid::textGrob("Gene",
    #                                                              gp = grid::gpar(fontsize=15)))
    # gridExtra::grid.arrange(TTESTPLOT, ncol = 1, heights = c(6))
    
    
    linearTest <- function(dat, one_reference=one_reference){
        plotting <- function(dat, one_reference){
            # run linear models
            test <- lapply(1:(ncol(dat)-1), function(i){
                summary(lm(unlist(dat[,i,with =FALSE])~ dat$cell_type))
            })
            coef <- data.frame(t(sapply(1:length(test), function(x){
                temp <- test[[x]]$coefficients
                estimates <- temp[,1]
            })))
            pval <- data.frame(t(sapply(1:length(test), function(x){
                temp <- test[[x]]$coefficients
                estimates <- temp[,4]
            })))
            
            coef$Gene <- colnames(dat)[-ncol(dat)]
            coefData <- melt(coef, id = "Gene")
            colnames(coefData) <- c("Gene", "cell_type", "Coefficient")
            coefPlot <- ggplot(data = coefData, aes(x = Gene, y = Coefficient, colour = cell_type, fill = cell_type,  group = length(unique(dat$cell_type)))) +
                geom_line(alpha = .3) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = names(cumulativeCounts[-1])) +
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "none",
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
            
            pval$Gene <- colnames(dat)[-ncol(dat)]
            pvalData <- melt(pval, id = "Gene")
            colnames(pvalData) <- c("Gene", "cell_type", "pVal")
            lmPvalPlot <- ggplot(data = pvalData, aes(x = Gene, y = -10*log10(pVal), colour = cell_type, fill = cell_type, group = length(unique(dat$cell_type)))) +
                geom_line(alpha = .3) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = names(cumulativeCounts[-1]))+
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "none",
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
            lmPlot_list <- list(coefPlot,lmPvalPlot)
            return(lmPlot_list)
        }
        lmPlots <- c()
        if (length(ref_groups) > 1 & one_reference == FALSE){
            dat$cell_type <- factor(dat$cell_type)
            types <- length(unique(dat$cell_type))
            refs <- sapply(1:length(ref_groups), function(i){paste("Reference",toString(i),sep = "")}) 
            names(refs) <- refs
            
            plotrefs <- lapply(refs, function(i){
                contr.treatment(types, base = which(sort(unique(dat$cell_type)) == i))})
            for (i in 1:length(plotrefs)){
                contrasts(dat$cell_type) <- plotrefs[[i]] 
                tempPlot <- plotting(dat, one_reference = one_reference)
                lmPlots <- c(lmPlots, tempPlot)}
            LMPLOT <- gridExtra::arrangeGrob( 
                grobs = lmPlots,
                ncol=1, # use apply and get to get variables for each step instead
                top    = grid::textGrob("Linear Regression Model Estimates and P-Values",
                                        gp = grid::gpar(fontsize=15,font=2)),
                bottom = grid::textGrob("Gene",
                                        gp = grid::gpar(fontsize=15)))
            plot_with_legend <- lmPlots[[1]] + theme(legend.position = "bottom") +
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(plot_with_legend)
            print(gridExtra::grid.arrange(LMPLOT,
                                          plot_legend, ncol = 1, heights = c(7,1)))
        } else {
            # need to make the reference the intercept 
            dat$cell_type <- factor(dat$cell_type)
            types <- length(unique(dat$cell_type))
            contrasts(dat$cell_type) <- contr.treatment(types, base = which(sort(unique(dat$cell_type)) == "Reference"))
            lmPlots <- plotting(dat, one_reference = one_reference)
            LMPLOT <- gridExtra::arrangeGrob( 
                grobs = lmPlots,
                ncol=1, # use apply and get to get variables for each step instead
                top    = grid::textGrob("Linear Regression Model Estimates and P-Values",
                                        gp = grid::gpar(fontsize=15,font=2)),
                bottom = grid::textGrob("Gene",
                                        gp = grid::gpar(fontsize=15)))
            plot_with_legend <- lmPlots[[1]] + 
                theme(legend.position = "bottom") +
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(plot_with_legend)
            print(gridExtra::grid.arrange(LMPLOT, plot_legend, ncol = 1, heights = c(10,1)))
        }
    }
    linearTest(dat, one_reference = one_reference)
    
    
    dev.off()
}

if (!is.null(args)) {
    CreatePlots(files_dir = args_parsed$file_dir,
                cluster_by_groups = args_parsed$cluster_by_groups,
                one_reference = args_parsed$one_reference)
}
