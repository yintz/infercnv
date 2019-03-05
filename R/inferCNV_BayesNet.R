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
#' 
# Requires:
#	infercnv, rjags, ggplot2, parallel, futile.logger, reshape
## build off of the present S4 object inferCNV_obj to add more slots 
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




#############
# Accessors #
#############

#' Access the values for cellGene 
#' 
#' This function returns the list of values in cellGene
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return A list.
#'
#' @exportMethod cellGene
#' @rdname cellGene-method
#' 
setGeneric(name = "cellGene", 
           def = function(obj) standardGeneric("cellGene"))
setMethod(f = "cellGene", 
          signature = "MCMC_inferCNV", 
          definition=function(obj) obj@cell_gene)


#######################
# Object Manipulation #
#######################
#' 
#' Get the cell Mean and Standard Deviation for identified cnv regions 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod MeanSD
#' @rdname MeanSD-method
#' 
setGeneric(name="MeanSD",
           def=function(obj)
               { standardGeneric("MeanSD") }
)

#' @rdname MeanSD-method
#' @aliases MeanSD
#' 
setMethod(f="MeanSD",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              gene_expr_by_cnv = .get_gene_expr_by_cnv(obj@.hspike)
              cnv_mean_sd = .get_gene_expr_mean_sd_by_cnv(gene_expr_by_cnv)
              cnv_sd <- cbind(lapply(cnv_mean_sd,function(x){x$sd}))
              cnv_mean <- cbind(lapply(cnv_mean_sd,function(x){x$mean}))
              ## Sort so in order of {x0,...,x1,..,x3} and get into a vector format
              obj@mu <- unlist(cbind(cnv_mean[sort(row.names(cnv_mean)),]))
              obj@sig <-  unlist(cbind(cnv_sd[sort(row.names(cnv_sd)),]))
              obj@sig <- 1/(obj@sig^2)
              if (obj@args$quietly == FALSE) {
                  print(paste("Means: ", obj@mu, collapse = ""))
                  print(paste("Sig: ", obj@sig, collapse = ""))
              }
              return(obj)
          }
)


#' Create a list that holds Genes and Cells for each separate identified CNV 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod getGenesCells
#' @rdname getGenesCells-method
#' 
setGeneric(name="getGenesCells",
           def=function(obj, pred_cnv_genes_df, cell_groups_df)
               { standardGeneric("getGenesCells") }
)

#' @rdname getGenesCells-method
#' @aliases getGenesCells
#' 
setMethod(f="getGenesCells",
          signature="MCMC_inferCNV",
          definition=function(obj, pred_cnv_genes_df, cell_groups_df)
          {
              ## list that holds Genes and Cells for each separate identified CNV 
              obj@cell_gene <- lapply(obj@cnv_regions,function(i) {
                  # subset the data to get the rows for the current CNV 
                  current_cnv <- pred_cnv_genes_df[which(i == pred_cnv_genes_df$gene_region_name),]
                  # get the index for the genes that are in each cnv 
                  genes <- current_cnv$gene
                  # pred_cnv_genes_df[which(pred_cnv_genes_df$gene_region_name %in% i),]$gene
                  gene_idx <- which(row.names(obj@expr.data) %in% genes)
                  
                  # get the index for the cells that are in each cnv 
                  sub_cells <- unique(current_cnv$cell_group_name)
                  cells_idx <- which(colnames(obj@expr.data) %in% cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
                  return(list("cnv_regions" = i, "Genes" = gene_idx, "Cells" = cells_idx))
              })
              return(obj)
          }
)




#' Initialize the MCMC_inferCNV_obj object 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod initializeObject
#' @rdname initializeObject-method
#' 
setGeneric(name="initializeObject",
           def=function(obj, args_parsed, infercnv_obj, HMM_obj)
               { standardGeneric("initializeObject") }
)

#' @rdname initializeObject-method
#' @aliases initializeObject
#' 
setMethod(f="initializeObject",
          signature="MCMC_inferCNV",
          definition=function(obj, args_parsed, infercnv_obj, HMM_obj)
          {
              futile.logger::flog.info(paste("Initializing new MCM InferCNV Object."))
              files <- list.files(args_parsed$file_dir, full.names = TRUE)
              
              # Validate the inferCNV Object 
              validate_infercnv_obj(infercnv_obj)
              
              ## create the S4 object
              obj <- MCMC_inferCNV(infercnv_obj)
              ## add the command line arguments 
              obj@args <- args_parsed
              
              ## Load the files for cnv predictions 
              cell_groups_PATH <- files[grep(files, pattern = "_HMM_preds.cell_groupings")]
              pred_cnv_genes_PATH <- files[grep(files, pattern = "_HMM_preds.pred_cnv_genes.dat")]
              # cell_groups_df <- read.table(cell_groups_PATH, header = T)
              cell_groups_df <- read.csv(cell_groups_PATH, sep = "\t", header = T, check.names = FALSE)
              pred_cnv_genes_df <- read.table(pred_cnv_genes_PATH, header = T, check.names = FALSE)
              
              # cnv region id's 
              obj@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
              futile.logger::flog.info(paste("Total CNV's: ", length(obj@cnv_regions)))
              
              ## Load Mixture Model File
              futile.logger::flog.info(paste("Loading BUGS Model."))
              obj@bugs_model <- readChar(obj@args$model_file,file.info(obj@args$model_file)$size)
              
              ## list that holds Genes and Cells for each separate identified CNV 
              obj <- getGenesCells(obj, pred_cnv_genes_df, cell_groups_df)
              
              # Create numerical ids for each subgroup of cells 
              ## group name ids
              cell_group_id <- unique(pred_cnv_genes_df$cell_group_name)
              
              ## set numerical id's for cell groups and set values in a vector for cell positions in the matrix 
              group_id <- rep(NA, max(unlist(obj@observation_grouped_cell_indices)))
              lapply(1:length(cell_group_id), function(i) {
                  ## cells in the cluster group 
                  cells <- cell_groups_df[cell_groups_df$cell_group_name %in% cell_group_id[i],]$cell
                  ## set the numerical id in the vector 
                  group_id[which(colnames(obj@expr.data) %in% cells)] <<- i
              })
              obj@group_id <- group_id
              
              return(obj)
          }
)


#' Get the state values from the inferCNV HMM object 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod getStates
#' @rdname getStates-method
#' 
setGeneric(name="getStates",
           def=function(obj, HMM_obj)
               { standardGeneric("getStates") }
)

#' @rdname getStates-method
#' @aliases getStates
#' 
setMethod(f="getStates",
          signature="MCMC_inferCNV",
          definition=function(obj, HMM_obj)
          {
              # Add the HMM defined states
              x <- HMM_obj@expr.data
              obj@States <- x
              return(obj)
          }
)


#' Set the probabilities for each CNV belonging to each state as well as probability of each cell belonging to a states 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod getProbabilities
#' @rdname getProbabilities-method
#' 
setGeneric(name="getProbabilities",
           def=function(obj)
               { standardGeneric("getProbabilities") }
)

#' @rdname getProbabilities-method
#' @aliases getProbabilities
#' 
setMethod(f="getProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              ## List holding state probabilities for each CNV 
              cnv_probabilities <- list()
              ## List for combining the chains in each simulation 
              combined_samples <- list()
              ## list holding the frequency of epsilon values for each cell line 
              ##  for each cnv region and subgroup 
              cell_probabilities <- list()
              
              for(j in 1:length(obj@mcmc)){
                  # combine the chains 
                  obj@mcmc[[j]] <- do.call(rbind, obj@mcmc[[j]])
                  # run function to get probabilities 
                  ## Thetas 
                  cnv_probabilities[[j]] <- cnv_prob(obj@mcmc[[j]])
                  ## Epsilons
                  cell_probabilities[[j]] <- cell_prob(obj@mcmc[[j]])
              }
              
              obj@cnv_probabilities <- cnv_probabilities
              obj@cell_probabilities <- cell_probabilities
              return(obj)
          }
)

#' Run simulations in Parallel
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod withParallel
#' @rdname withParallel-method
#' 
setGeneric(name="withParallel",
           def=function(obj)
               { standardGeneric("withParallel") }
)

#' @rdname withParallel-method
#' @aliases withParallel
#' 
setMethod(f="withParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Running Sampling Using Parallel with ",obj@args$CORES,"Cores"))
              obj@mcmc <- parallel::mclapply(1:length(obj@cell_gene), function(i){ 
                  if (obj@args$quietly == FALSE) {
                      futile.logger::flog.info(paste("Sampleing Number: ", i))
                  }
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- obj@group_id[obj@cell_gene[[i]]$Cells] # subset the tumor ids for the cells wanted
                      gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes, obj@cell_gene[[i]]$Cells]
                      return(run_gibb_sampling(gene_exp, obj))
                  } else {
                      return(list(NULL)) 
                  }
              },
              mc.cores=as.integer(obj@args$CORES))
              return(obj)
          }
)

#' Run simulations in Non-Parallel mode 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod nonParallel
#' @rdname nonParallel-method
#' 
setGeneric(name="nonParallel",
           def=function(obj)
               { standardGeneric("nonParallel") }
)

#' @rdname nonParallel-method
#' @aliases nonParallel
#' 
setMethod(f="nonParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Running Gibbs sampling in Non-Parallel Mode."))
              # Iterate over the CNV's and run the Gibbs sampling.
              obj@mcmc <- lapply(1:length(obj@cell_gene), function(i){
                  if (obj@args$quietly == FALSE) {
                      print(obj@args$quietly)
                      futile.logger::flog.info(paste("Sample Number: ", i))
                  }
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- obj@group_id[ obj@cell_gene[[i]]$Cells ] # subset the tumor ids for the cells wanted
                      gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes, obj@cell_gene[[i]]$Cells]
                      return(run_gibb_sampling(gene_exp, obj))
                  } else {
                      return(list(NULL))
                  }
              })
              return(obj)
          }
)

#' Run simulations and remove CNV's that have a probability of being normal above a set thresholld. 
#' This removes possible false posotives identified by the HMM. 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod removeCNV
#' @rdname removeCNV-method
#' 
setGeneric(name="removeCNV",
           def=function(obj)
               { standardGeneric("removeCNV") }
)

#' @rdname removeCNV-method
#' @aliases removeCNV
#' 
setMethod(f="removeCNV",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Mean values of the probability distribution of the CNV states p(CNV == {states 1:6})
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              if (any(cnv_means[3,] > obj@args$BayesMaxPNormal)){
                  remove_cnv <- which(cnv_means[3,] > obj@args$BayesMaxPNormal)
                  futile.logger::flog.info(paste("Removing ",length(remove_cnv), " CNV(s) identified by the HMM."))
                  
                  if (obj@args$quietly == FALSE) { print("CNV's being removed have the following posterior probabilities of being a normal state: ") }
                  
                  lapply(remove_cnv, function(i) {
                      if (obj@args$quietly == FALSE) {
                          print( paste(obj@cell_gene[[i]]$cnv_regions, ", Genes: ", length(obj@cell_gene[[i]]$Genes), " Cells: ", length(obj@cell_gene[[i]]$Cells)) )
                        }
                      ## Change the states to normal states
                      obj@States[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- 3
                  })
                  ## Remove the CNV's from the following matrices 
                  obj@cell_gene <- obj@cell_gene[-remove_cnv]
                  obj@cell_probabilities <- obj@cell_probabilities[-remove_cnv]
                  obj@cnv_probabilities <- obj@cnv_probabilities[-remove_cnv]
                  cnv_means <- cnv_means[,-remove_cnv]
              }
              
              # Write the state probabilities for each CNV to a table. 
              ## set column names to the CNV ID
              cnv_regions <- sapply(obj@cell_gene, function(i) { as.character(i$cnv_regions) })
              colnames(cnv_means) <- cnv_regions
              ## set row names to the states 1:6
              row.names(cnv_means) <- c(sprintf("State:%s",1:6))
              write.table(cnv_means,file = "CNV_State_Probabilities.dat", col.names = TRUE, row.names=TRUE, quote=FALSE, sep="\t")
              return(obj)
          }
)

#' Run simulations and remove cells from cnv's that are predicted to be normal 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod removeCells
#' @rdname removeCells-method
#' 
setGeneric(name="removeCells",
           def=function(obj)
               { standardGeneric("removeCells") }
)

#' @rdname removeCells-method
#' @aliases removeCells
#' 
setMethod(f="removeCells",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              if (any(do.call(cbind, obj@cell_probabilities)[3,] > obj@args$BayesMaxPNormal)){
                  lapply(1:length(obj@cell_probabilities), function(i) {
                      idx <- which(obj@cell_probabilities[[i]][3,] > obj@args$BayesMaxPNormal)
                      if(length(idx) > 0){
                          ## change the states to normal states
                          obj@States[ obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells[ idx ] ] <<- 3
                          ## remove these cells from the cnv
                          obj@cell_gene[[i]]$Cells <<- obj@cell_gene[[i]]$Cells[- idx]
                      }
                  })
                  # recursively run again 
                  obj <- runMCMC(obj)
              }
              return(obj)
          }
)

#' Run simulations using rjags.
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod runMCMC
#' @rdname runMCMC-method
#' 
setGeneric(name="runMCMC",
           def=function(obj)
               { standardGeneric("runMCMC") }
)

#' @rdname runMCMC-method
#' @aliases runMCMC
#' 
setMethod(f="runMCMC",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Run MCMC
              if(is.null(obj@args$CORES)){
                  obj <- nonParallel(obj)
              } else {
                  obj <- withParallel(obj)
              }
              
              # Get the probability of of each cell line and complete CNV belonging to a specific state
              obj <- getProbabilities(obj)
              
              if(!(is.null(obj@args$postMcmcMethod))){
                  if(obj@args$postMcmcMethod == "removeCNV"){
                      obj <- removeCNV(obj)
                  } else {
                      obj <- removeCells(obj)
                  }
              }
              
              return(obj)
          }
)


##########################
# Plotting functions #
##########################
#' Get the probability of each cnv being a normal state and plot these probabilities. 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod postProbNormal
#' @rdname postProbNormal-method
#' 
setGeneric(name="postProbNormal",
           def=function(obj)
               { standardGeneric("postProbNormal") }
)

#' @rdname postProbNormal-method
#' @aliases postProbNormal
#' 
setMethod(f="postProbNormal",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # get probability of the cnv's belonging to each state 
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              # Adjust the probabilities so greater probability corresponds to less likely to be normal 
              normal_prob <- 1 - cnv_means[3,]
              obj@expr.data[,] <- 0
              lapply(1:length(normal_prob), function(i) { 
                  ## change the states to normal states
                  obj@expr.data[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- normal_prob[i]
              })
              infercnv::plot_cnv(infercnv_obj          = obj,
                                 #k_obs_groups         = 4,
                                 #cluster_by_groups    = cluster_by_groups,
                                 title                 = sprintf("NormalProbabilities"),
                                 output_filename       = file.path(file.path(obj@args$out_dir),"infercnv.NormalProbabilities"),
                                 write_expr_matrix     = FALSE,
                                 x.center              = 0,
                                 x.range               = c(0,1)
              )
          }
)

#' Plots the probability for each cnv belonging to a specific state and the probability of 
#' each cell line belonging to a specific states. 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @return obj The MCMC_inferCNV_obj S4 object.
#' 
#' @exportMethod plotProbabilities
#' @rdname plotProbabilities-method
#' 
setGeneric(name="plotProbabilities",
           def=function(obj)
               { standardGeneric("plotProbabilities") }
)

#' @rdname plotProbabilities-method
#' @aliases plotProbabilities
#' 
setMethod(f="plotProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Plotting
              ## plots the probability of each cell line being a particular state
              ## plots the probability of a cnv being a particular state 
              
              # Plot the probabilities of epsilons 
              ep <- function(df){df[,grepl('epsilon', colnames(df))]}
              epsilons <- lapply(obj@mcmc, function(x) ep(x))
              
              pdf(file = file.path(file.path(obj@args$out_dir),"cellProbs.pdf"), onefile = TRUE)
              lapply(1:length(obj@cell_probabilities), function(i){ 
                  print(plot_cell_prob(as.data.frame(obj@cell_probabilities[[i]]), as.character(obj@cell_gene[[i]]$cnv_regions)))
              })
              dev.off()
              
              ## Plot the probability of each state for a CNV 
              pdf(file = file.path(file.path(obj@args$out_dir),"cnvProbs.pdf"), onefile = TRUE)
              lapply(obj@cnv_probabilities,function(x){
                  print(plot_cnv_prob(x))
              })
              dev.off()
          }
)


#' Return the InferCNV Object with the new adjucted CNV's 
#' 
#' Returns Infercnv Object 
#' 
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return An inferCNV object
#'
#' @exportMethod returningInferCNV
#' @rdname returningInferCNV-method
#' 
setGeneric(name = "returningInferCNV", 
           def = function(obj, infercnv_obj) 
               { standardGeneric("returningInferCNV") }
)
setMethod(f = "returningInferCNV", 
          signature = "MCMC_inferCNV", 
          definition=function(obj, infercnv_obj) {
              NewStates <- obj@States
              infercnv_obj@expr.data <- NewStates
              return(infercnv_obj)
          }
)



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
                              dest="out_dir",
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
    if (MCMC_inferCNV_obj@args$quietly == FALSE) {
        futile.logger::flog.info(paste("Cells: ",C))
        futile.logger::flog.info(paste("Genes: ",G))
    }
    # quiet=FALSE
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
                               quiet=MCMC_inferCNV_obj@args$quietly)
    update(model, 200, progress.bar=ifelse(MCMC_inferCNV_obj@args$quietly,"none","text"))
    # run the rjags model 
    ## set the parameters to return from sampling 
    parameters <- c('theta', 'epsilon')#, 'gamma')
    samples <- rjags::coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(MCMC_inferCNV_obj@args$quietly,"none","text"))
    return(samples)
}
##########################################################################################################################################################

# Function to plot the probability for each cell line of being in a particular state 
plot_cell_prob <- function(df, title){
    df$mag = c(1:6)
    long_data <- reshape::melt(df, id = "mag")
    ggplot2::ggplot(long_data, ggplot2::aes_string(x = 'variable', y = value, fill = as.factor('mag')))+
        ggplot2::geom_bar(stat="identity", width = 1) +
        ggplot2::coord_flip() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),panel.border = ggplot2::element_blank(),
            axis.text=ggplot2::element_text(size=20),
            plot.title = ggplot2::element_text(hjust = 0.5,size = 22),
            #legend.position = "none",
            legend.position="bottom",
            axis.text.x = ggplot2::element_text(size = 16),
            axis.text.y = ggplot2::element_text(size = 16),
            axis.title.x = ggplot2::element_text(size = 18),
            axis.title.y = ggplot2::element_text(size = 18))+
        ggplot2::labs(title = title) +
        #fill = "CNV States") + 
        ggplot2::xlab("Cell") +
        ggplot2::ylab("Probability")+
        ggplot2::scale_x_discrete(breaks =seq(1, ncol(df), 9))
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
    ggplot2::ggplot(data = means, ggplot2::aes_string(y = 'Probability', x= 'State', fill = as.factor('State'))) +
        ggplot2::geom_bar(stat = "identity")
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
#' @return Returns a MCMC_inferCNV_obj and posterior probability of being in one of six Copy Number Variation states 
#' (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by inferCNV's HMM. 
#' @export

inferCNVBayesNet <- function( 
                              file_dir,
                              infercnv_obj,
                              HMM_obj,
                              BayesMaxPNormal,
                              model_file = system.file("BUGS_Mixture_Model",package = "infercnv"),
                              CORES = NULL,
                              out_dir,
                              postMcmcMethod = NULL,
                              plotingProbs = TRUE,
                              quietly = TRUE) {
    
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
        if (as.integer(CORES) > detectCores()){
            error_message <- paste("Too many cores previded. The following system has ",detectCores(), " cores.",
                                   "Please select an appropriate amount.")
            futile.logger::flog.error(error_message)
            stop(error_message)
        } 
    }
    
    if(out_dir != "." & !file.exists(out_dir)){
        # create the output directory 
        dir.create(file.path(out_dir))
        futile.logger::flog.info(paste("Creating the following Directory: ", out_dir))
    }
    args_parsed <- list("file_dir" = file_dir,
                        "model_file" = model_file,
                        "CORES" = CORES,
                        "out_dir"= out_dir,
                        "plotingProbs" = TRUE,
                        "postMcmcMethod"=postMcmcMethod,
                        "BayesMaxPNormal" = BayesMaxPNormal,
                        "quietly" = quietly)
    #################################
    # LOAD DATA & INITIALIZE OBJECT #
    #################################
    
    ## create the S4 object
    MCMC_inferCNV_obj <- new("MCMC_inferCNV")
    MCMC_inferCNV_obj <- initializeObject(MCMC_inferCNV_obj, args_parsed, infercnv_obj)
    MCMC_inferCNV_obj <- getStates(MCMC_inferCNV_obj, HMM_obj)
    
    #############
    # MEAN & SD #
    #############
    ## Get mean and sd of expression in the predicted cnv areas 
    MCMC_inferCNV_obj <- MeanSD(MCMC_inferCNV_obj)
    
    
    # check and print the number of genes in each cnv 
    if (args_parsed$quietly == FALSE) {
        number_of_genes <- sapply(cellGene(MCMC_inferCNV_obj), function(i) length(i$Genes))
        print(paste("Number of genes for each CNV: ", paste(number_of_genes, sep = " ",collapse = " ")))
        # check the lengths of each cell group 
        sapply(cellGene(MCMC_inferCNV_obj), function(x) length(x$Cells))
    }
    
    ################################
    # Run MCMC Sampling            #
    ################################
    start_time <- Sys.time()
    MCMC_inferCNV_obj <- runMCMC(MCMC_inferCNV_obj)
    end_time <- Sys.time()
    futile.logger::flog.info("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes")
    
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
    
    ##############################
    # Return new infercnv object #
    ##############################
    infercnv_obj <- returningInferCNV(MCMC_inferCNV_obj, infercnv_obj)
    
    return(infercnv_obj)
}

##########################
# Command Line Arguments #
##########################
## Uncomment to use the command line arguments 
# if (!is.null(args)){
#     # Set Constants 
#     args_parsed <- optparse::parse_args(pargs)
#     file_dir <- args_parsed$file_dir
#     model_file <- args_parsed$model_file
#     CORES <- args_parsed$CORES
#     out_dir <- args_parsed$out_dir
#     inferCNVBayesNet(file_dir,
#                      model_file,
#                      CORES,
#                      out_dir)
# } 