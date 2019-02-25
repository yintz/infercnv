#!/usr/bin/env Rscript
#

##################################
# create MCMC_inferCNV S4 object #
##################################

#' MCMC_inferCNV class
#' 
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior 
#' probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by 
#' inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster and each individual
#' cell line in the CNV. 
#' 
#' @aliases infercnv-package
#'
#' @slot bugs_model BUGS model. 
#' @slot sig fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line 
#' @slot mu Mean values to be used for determining the distribution of each cell line 
#' @slot group_id ID's given to the cell clusters.
#' @slot cell_gene List containing the Cells and Genes that make up each CNV.
#' @slot mcmc Simulation output from sampling.
#' @slot cnv_probabilities Probabilities of each CNV belonging to a particular state from 0 (least likely)to 1 (most likely). 
#' @slot cell_probabilities Probabilities of each cell being in a particular state, from 0 (least likely)to 1 (most likely).
#' @slot args Input arguments given by the user 
#' @slot cnv_regions ID for each CNV found by the HMM
#' @slot States States that are identified and (depending on posterior MCMC input methods) modified.
#'
#'
#'
#'
#' @return Returns a MCMC_inferCNV_obj
#' @export
MCMC_inferCNV <- setClass("MCMC_inferCNV", slots = c(bugs_model = "character",
                                                     sig = "numeric",
                                                     mu = "numeric",
                                                     group_id = "integer",
                                                     cell_gene = "list",
                                                     mcmc = "list",
                                                     cnv_probabilities = "list",
                                                     cell_probabilities = "list",
                                                     args = "list",
                                                     cnv_regions = "factor",
                                                     States = "ANY"),
                          contains = "infercnv")

##########################
# Command line arguments #
##########################
pargs <- optparse::OptionParser()
pargs <- optparse::add_option(pargs, c("-f", "--infercnv_dir"),
                              type="character",
                              action="store",
                              dest="file_dir",
                              metavar="File_Directory",
                              help=paste("Path to files created by inferCNV.",
                                         "[Default %default][REQUIRED]"))
pargs <- optparse::add_option(pargs, c("-m", "--model"),
                              type="character",
                              action="store",
                              dest="model_file",
                              metavar="Model_File_Path",
                              help=paste("Path to the BUGS Model file.",
                                         "[Default %default][REQUIRED]"))
pargs <- optparse::add_option(pargs, c("-p","--parallel"),
                              type="character",
                              action="store",
                              dest="CORES",
                              default = NULL,
                              metavar="Number_of_Cores",
                              help=paste("Option to run parallel by specifying the number of cores to be used.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("-o","--out_dir"),
                              type="character",
                              action="store",
                              dest="output_dir",
                              default = NULL,
                              metavar="Output_Directory",
                              help=paste("Option to set the output directory to save the RDS object.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("-M","--method"),
                              type="character",
                              action="store",
                              dest="postMcmcMethod",
                              default = NULL,
                              metavar="Posterior_MCMC_Method",
                              help=paste("What actions to take after finishing the MCMC.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("-x","--plot"),
                              type="logical",
                              action="store",
                              dest="plotingProbs",
                              default = TRUE,
                              metavar="Plot_Probabilities",
                              help=paste("Plot the posterior probabilites for each CNV and each cell line in each cnv.",
                                         "[Default %default]"))

# Function to run the mixture model for given expression data 
# Runs on each cnv seperately. Then seperate by tumor subgroup.
#       input: 
#               gene_exp            : Gene expression data
#               MCMC_inferCNV_obj   : MCMC_inferCNV object 
#       return:
#               samples     : Results of the sampleing process 
#
run_gibb_sampling <- function( 	gene_exp, 
                                MCMC_inferCNV_obj
){
    if(is.null(ncol(gene_exp))){
        gene_exp <- data.frame(gene_exp)
    } 
    C = ncol(gene_exp)
    G = nrow(gene_exp)
    futile.logger::flog.info(paste("Cells: ",C))
    futile.logger::flog.info(paste("Genes: ",G))
    quiet=FALSE
    # make Data list for model 
    data <- list(
        'C' = C,                                 # number of cell lines 
        'G' = G,                                 # number of genes 
        'gexp' = gene_exp,                       # expression data
        'sig' = MCMC_inferCNV_obj@sig,           # fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line 
        'mu' = MCMC_inferCNV_obj@mu              # Mean values to be used for determining the distribution of each cell line 
    )
    # set initial values for each cell line begining states 
    inits <- list(
        list(epsilon = rep(1, C)),
        list(epsilon = rep(2, C)),
        list(epsilon = rep(3, C)),
        list(epsilon = rep(4, C)),
        list(epsilon = rep(5, C)),
        list(epsilon = rep(6, C))
    )
    # Create the model for rjags
    model <- rjags::jags.model(textConnection(MCMC_inferCNV_obj@bugs_model), 
                        data=data, 
                        inits=inits, # (Initialization) optional specification of initial values in the form of a list or a function
                        n.chains=6,  # the number of parallel chains for the model
                        n.adapt=500, # the number of iterations for adaptation (burn in)
                        quiet=FALSE)
    update(model, 200, progress.bar=ifelse(quiet,"none","text"))
    # run the rjags model 
    ## set the parameters to return from sampling 
    parameters <- c('theta', 'epsilon')#, 'gamma')
    samples <- rjags::coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(quiet,"none","text"))
    return(samples)
}
##########################################################################################################################################################

# Function to plot the probability for each cell line of being in a particular state 
plot_cell_prob <- function(df, title){
    df$mag = c(1:6)
    long_data <- reshape::melt(df, id = "mag")
    ggplot2::ggplot(long_data, aes(x = variable, y = value, fill = as.factor(mag)))+
        geom_bar(stat="identity", width = 1) +
        coord_flip() +
        theme(
            panel.grid = element_blank(), panel.background = element_blank(),panel.border = element_blank(),
            axis.text=element_text(size=20),
            plot.title = element_text(hjust = 0.5,size = 22),
            #legend.position = "none",
            legend.position="bottom",
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18))+
        labs(title = title) +
        #fill = "CNV States") + 
        xlab("Cell") +
        ylab("Probability")+
        scale_x_discrete(breaks =seq(1, ncol(df), 9))
}

# Function for total CNV probaility of belonging to each state using THETA prior 
cnv_prob <- function(combined_samples) {
    thetas <- combined_samples[,grepl('theta', colnames(combined_samples))]
    #print(paste("Thetas: ", dim(thetas)))
    return(thetas)
}

# Function for each individule cell probabilities, marginalize over the EPSILONS
cell_prob <- function(combined_samples) {
    epsilons <- combined_samples[,grepl('epsilon', colnames(combined_samples))]
    #print(paste("Epsilons: ", dim(epsilons)))
    epsilon_state_frequencies <- apply(as.data.frame(epsilons), 2, function(x) table(factor(x, levels=1:6)))
    cell_probs <- epsilon_state_frequencies/colSums(epsilon_state_frequencies)
    return(cell_probs)
}

## Fucntion to Plot the probability of each state for a CNV 
plot_cnv_prob <- function(df){
    means <- as.data.frame(colMeans(df))
    means$state <- c(1:6)
    colnames(means) <- c("Probability", "State")
    ggplot2::ggplot(data = means, aes(y = Probability, x= State, fill = as.factor(State))) +
        geom_bar(stat = "identity")
}


#################
# main function #
#################
#' @title inferCNVBayesNet: Run Bayesian Network Mixture Model To Obtain Posterior Probabilities For HMM Predicted States
#' 
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior 
#' probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by 
#' inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster and each individual
#' cell line in the CNV. 
#'
#' @param infercnv_dir Location of the directory of the inferCNV outputs.
#' @param model Path to the BUGS Model file.
#' @param parallel Option to run parallel by specifying the number of cores to be used.
#' @param out_dir (string) Path to where the output file should be saved to.
#' @param method What actions to take after finishing the MCMC.
#'

inferCNVBayesNet <- function( 
                              file_dir,
                              model_file,
                              CORES = NULL,
                              output_dir) {
    
    ################
    # CHECK INPUTS #
    ################
    if (!file.exists(file_dir)){
        error_message <- paste("Cannot find the supplied directory location for the infercnv output.",
                               "Please supply teh correct path for the output.")
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    if (!file.exists(model_file)){
        error_message <- paste("Cannot find the model file.",
                               "Please supply the correct path for the model file.")
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    if (!is.null(CORES)){
        if (as.integer(CORES) > parallel::detectCores()){
            error_message <- paste("Too many cores previded. The following system has ",detectCores(), " cores.",
                                   "Please select an appropriate amount.")
            futile.logger::flog.error(error_message)
            stop(error_message)
        } 
    }
    
    if(output_dir != "." & !file.exists(output_dir)){
        # create the output directory 
        dir.create(file.path(output_dir))
        futile.logger::flog.info(paste("Creating the following Directory: ", out_dir))
    }
    args_parsed <- list("file_dir" = file_dir,
                      "model_file" = model_file,
                      "CORES" = CORES,
                      "output_dir"= output_dir)
    #################################
    # LOAD DATA & INITIALIZE OBJECT #
    #################################
    
    ## create the S4 object
    MCMC_inferCNV_obj <- new("MCMC_inferCNV")
    MCMC_inferCNV_obj <- initializeObject(MCMC_inferCNV_obj, args_parsed)
    MCMC_inferCNV_obj <- getStates(MCMC_inferCNV_obj)
    
    #############
    # MEAN & SD #
    #############
    ## Get mean and sd of expression in the predicted cnv areas 
    MCMC_inferCNV_obj <- MeanSD(MCMC_inferCNV_obj)
    
    
    # check and print the number of genes in each cnv 
    number_of_genes <- sapply(cellGene(MCMC_inferCNV_obj), function(i) length(i$Genes))
    print(paste("Number of genes for each CNV: ", paste(number_of_genes, sep = " ",collapse = " ")))
    # check the lengths of each cell group 
    sapply(cellGene(MCMC_inferCNV_obj), function(x) length(x$Cells))
    
    ################################
    # Run MCMC Sampling            #
    ################################
    MCMC_inferCNV_obj <- runMCMC(MCMC_inferCNV_obj)
    
    
    ########
    # Plot #
    ########
    
    ## Chromosomes in which the cnv's are located in 
    # chr_per_cnv <- lapply(cell_gene,function(x){
    #     genes <- x$Genes
    #     infercnv_obj@gene_order[which(row.names(infercnv_obj@gene_order) %in% genes),]$chr
    # })
    
    ## Plot the resuls 
    if(args_parsed$plotingProbs == TRUE){
        plotProbabilities(MCMC_inferCNV_obj)
        postProbNormal(MCMC_inferCNV_obj)
    }
}

##########################
# Command Line Arguments #
##########################
# if (!is.null(args)){
#     # Set Constants 
#     args_parsed <- optparse::parse_args(pargs)
#     file_dir <- args_parsed$file_dir
#     model_file <- args_parsed$model_file
#     CORES <- args_parsed$CORES
#     output_dir <- args_parsed$output_dir
#     inferCNVBayesNet(file_dir,
#                      model_file,
#                      CORES,
#                      output_dir)
# } 