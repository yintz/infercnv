#!/usr/bin/env Rscript

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
#' @slot bugs_model BUGS model.
#' @slot sig fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line
#' @slot mu Mean values to be used for determining the distribution of each cell line
#' @slot group_id ID's given to the cell clusters.
#' @slot cell_gene List containing the Cells and Genes that make up each CNV.
#' @slot mcmc Simulation output from sampling.
#' @slot combined_mcmc Combined chains for simulation output from sampling.
#' @slot cnv_probabilities Probabilities of each CNV belonging to a particular state from 0 (least likely)to 1 (most likely).
#' @slot cell_probabilities Probabilities of each cell being in a particular state, from 0 (least likely)to 1 (most likely).
#' @slot args Input arguments given by the user
#' @slot cnv_regions ID for each CNV found by the HMM
#' @slot States States that are identified and (depending on posterior MCMC input methods) modified.
#'
#' @exportClass MCMC_inferCNV
#' @name MCMC_inferCNV-class
#' @rdname MCMC_inferCNV-class
#' @keywords classes
#'
# Requires:
# infercnv, rjags, ggplot2, parallel, futile.logger, reshape
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
                                                     States = "ANY",
                                                     combined_mcmc = "list"),
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
#' @rdname cellGene-method
#' @keywords internal
#' @noRd
setGeneric(name = "cellGene",
           def = function(obj) standardGeneric("cellGene"))
#' @rdname cellGene-method
#' @aliases cellGene
#' @noRd
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
#' @rdname MeanSD-method
#' @keywords internal
#' @noRd
setGeneric(name="MeanSD",
           def=function(obj)
               { standardGeneric("MeanSD") }
)

#' @rdname MeanSD-method
#' @aliases MeanSD
#' @noRd
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

#' Add the probability threshold for the arguments in the MCMC infercnv object.
#'
#' This function adds the variable BayesMaxPNormal to the arguments slot of the the MCMC infercnv object.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param BayesMaxPNormal probability to be used as a threshold for CNV or cell removal.
#'
#' @return MCMC_inferCNV_obj S4 object.
#'
#' @rdname setBayesMaxPNormal-method
#' @keywords internal
#' @noRd
setGeneric(name = "setBayesMaxPNormal",
           def = function(obj, BayesMaxPNormal) standardGeneric("setBayesMaxPNormal"))

#' @rdname setBayesMaxPNormal-method
#' @aliases setBayesMaxPNormal
#' @noRd
setMethod(f = "setBayesMaxPNormal",
          signature = "MCMC_inferCNV",
          definition=function(obj, BayesMaxPNormal) {
              obj@args$BayesMaxPNormal <- BayesMaxPNormal
              return(obj)
              })

#' Create a list that holds Genes and Cells for each separate identified CNV
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param pred_cnv_genes_df Data for genes in each predicted CNV.
#' @param cell_groups_df Data for each cell in the predicted CNV's.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname getGenesCells-method
#' @keywords internal
#' @noRd
setGeneric(name="getGenesCells",
           def=function(obj, pred_cnv_genes_df, cell_groups_df)
               { standardGeneric("getGenesCells") }
)

#' @rdname getGenesCells-method
#' @aliases getGenesCells
#' @noRd
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
#' @param args_parsed The arguments given to the function.
#' @param infercnv_obj InferCNV object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname initializeObject-method
#' @keywords internal
#' @noRd
setGeneric(name="initializeObject",
           def=function(obj, args_parsed, infercnv_obj)
               { standardGeneric("initializeObject") }
)

#' @rdname initializeObject-method
#' @aliases initializeObject
#' @noRd
setMethod(f="initializeObject",
          signature="MCMC_inferCNV",
          definition=function(obj, args_parsed, infercnv_obj)
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
              cell_groups_df <- read.table(cell_groups_PATH, header = T, check.names = FALSE, sep="\t")
              pred_cnv_genes_df <- read.table(pred_cnv_genes_PATH, header = T, check.names = FALSE, sep="\t")

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
#' @param HMM_obj The HMM inferCNV object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname getStates-method
#' @keywords internal
#' @noRd
setGeneric(name="getStates",
           def=function(obj, HMM_obj)
               { standardGeneric("getStates") }
)

#' @rdname getStates-method
#' @aliases getStates
#' @noRd
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
#' @rdname getProbabilities-method
#' @keywords internal
#' @noRd
setGeneric(name="getProbabilities",
           def=function(obj)
               { standardGeneric("getProbabilities") }
)

#' @rdname getProbabilities-method
#' @aliases getProbabilities
#' @noRd
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

              combinedMCMC <-
              for(j in 1:length(obj@mcmc)){
                  # combine the chains
                  obj@combined_mcmc[[j]] <- do.call(rbind, obj@mcmc[[j]])
                  # run function to get probabilities
                  ## Thetas
                  cnv_probabilities[[j]] <- cnv_prob(obj@combined_mcmc[[j]])
                  ## Epsilons
                  cell_probabilities[[j]] <- cell_prob(obj@combined_mcmc[[j]])
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
#' @rdname withParallel-method
#' @keywords internal
#' @noRd
setGeneric(name="withParallel",
           def=function(obj)
               { standardGeneric("withParallel") }
)

#' @rdname withParallel-method
#' @aliases withParallel
#' @noRd
setMethod(f="withParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              par_func <- function(i){
                  if (obj@args$quietly == FALSE) {
                      futile.logger::flog.info(paste("Sampleing Number: ", i))
                  }
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- obj@group_id[ obj@cell_gene[[i]]$Cells ] # subset the tumor ids for the cells wanted
                      gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes, obj@cell_gene[[i]]$Cells]
                      return(run_gibb_sampling(gene_exp, obj))
                  } else {
                      return(list(NULL))
                  }
              }
              mc.cores = ifelse(.Platform$OS.type == 'unix', as.integer(obj@args$CORES), 1) # if windows, can only use 1 here
              futile.logger::flog.info(paste("Running Sampling Using Parallel with ", obj@args$CORES, "Cores"))
              obj@mcmc <- parallel::mclapply(1:length(obj@cell_gene),
                                             FUN = par_func,
                                             mc.cores = mc.cores)
              return(obj)
          }
)

#' Run simulations in Non-Parallel mode
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname nonParallel-method
#' @keywords internal
#' @noRd
setGeneric(name="nonParallel",
           def=function(obj)
               { standardGeneric("nonParallel") }
)

#' @rdname nonParallel-method
#' @aliases nonParallel
#' @noRd
setMethod(f="nonParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Running Gibbs sampling in Non-Parallel Mode."))
              # Iterate over the CNV's and run the Gibbs sampling.
              obj@mcmc <- lapply(1:length(obj@cell_gene), function(i){
                  if (obj@args$quietly == FALSE) {
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
#' @rdname removeCNV-method
#' @keywords internal
#' @noRd
setGeneric(name="removeCNV",
           def=function(obj)
               { standardGeneric("removeCNV") }
)

#' @rdname removeCNV-method
#' @aliases removeCNV
#' @noRd
setMethod(f="removeCNV",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Mean values of the probability distribution of the CNV states p(CNV == {states 1:6})
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              futile.logger::flog.info(paste("Attempting to removing CNV(s) with a probability of being normal above ", obj@args$BayesMaxPNormal))
              futile.logger::flog.info(paste("Removing ",length(which(cnv_means[3,] > obj@args$BayesMaxPNormal)), " CNV(s) identified by the HMM."))
              if (any(cnv_means[3,] > obj@args$BayesMaxPNormal)){
                  remove_cnv <- which(cnv_means[3,] > obj@args$BayesMaxPNormal)

                  if (obj@args$quietly == FALSE) { print("CNV's being removed have the following posterior probabilities of being a normal state: ") }

                  lapply(remove_cnv, function(i) {
                      if (obj@args$quietly == FALSE) {
                          print( paste(obj@cell_gene[[i]]$cnv_regions, ", Genes: ", length(obj@cell_gene[[i]]$Genes), " Cells: ", length(obj@cell_gene[[i]]$Cells)) )
                          # print(paste(paste( "Probabilities: "), cnv_means[,i]))
                        }
                      ## Change the states to normal states
                      obj@States[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- 3
                  })
                  ## Remove the CNV's from the following matrices
                  obj@cell_gene <- obj@cell_gene[-remove_cnv]
                  obj@cell_probabilities <- obj@cell_probabilities[-remove_cnv]
                  obj@cnv_probabilities <- obj@cnv_probabilities[-remove_cnv]
                  cnv_means <- cnv_means[,-remove_cnv]
                  # obj@mcmc <- obj@mcmc[-remove_cnv]
                  # obj@combined_mcmc <- obj@combined_mcmc[-remove_cnv]
                  futile.logger::flog.info(paste("Total CNV's after removing: ", length(obj@cell_gene)))
              }

              # Write the state probabilities for each CNV to a table.
              ## set column names to the CNV ID
              cnv_regions <- sapply(obj@cell_gene, function(i) { as.character(i$cnv_regions) })
              colnames(cnv_means) <- cnv_regions
              ## set row names to the states 1:6
              row.names(cnv_means) <- c(sprintf("State:%s",1:6))
              write.table(cnv_means,file = file.path(obj@args$out_dir, "CNV_State_Probabilities.dat"), col.names = TRUE, row.names=TRUE, quote=FALSE, sep="\t")
              return(obj)
          }
)

#' Run simulations and remove cells from cnv's that are predicted to be normal
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname removeCells-method
#' @keywords internal
#' @noRd
setGeneric(name="removeCells",
           def=function(obj)
               { standardGeneric("removeCells") }
)

#' @rdname removeCells-method
#' @aliases removeCells
#' @noRd
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
#' Run MCMC simulations using rjags. Also returns a plot the probability of each CNV being
#' normal before running any kind of post MCMC modification.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname runMCMC-method
#' @keywords internal
#' @noRd
setGeneric(name="runMCMC",
           def=function(obj)
               { standardGeneric("runMCMC") }
)

#' @rdname runMCMC-method
#' @aliases runMCMC
#' @noRd
setMethod(f="runMCMC",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Run MCMC
              if(obj@args$CORES == 1){
                  obj <- nonParallel(obj)
              } else {
                  obj <- withParallel(obj)
              }

              # Get the probability of of each cell line and complete CNV belonging to a specific state
              obj <- getProbabilities(obj)

              return(obj)
          }
)


##########################
# Plotting functions #
##########################
#' Get the probability of each cnv being a normal state and plot these probabilities.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param PNormal Option to add specific title to plot.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname postProbNormal-method
#' @keywords internal
#' @noRd
setGeneric(name="postProbNormal",
           def=function(obj, PNormal)
               { standardGeneric("postProbNormal") }
)

#' @rdname postProbNormal-method
#' @aliases postProbNormal
#' @noRd
setMethod(f="postProbNormal",
          signature="MCMC_inferCNV",
          definition=function(obj, PNormal)
          {
              if (obj@args$plotingProbs == TRUE){
                  # get probability of the cnv's belonging to each state
                  cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
                  # Adjust the probabilities so greater probability corresponds to less likely to be normal
                  normal_prob <- 1 - cnv_means[3,]
                  obj@expr.data[,] <- 0
                  lapply(1:length(normal_prob), function(i) {
                      ## change the states to normal states
                      obj@expr.data[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- normal_prob[i]
                  })
                  if (!is.null(PNormal)){
                      title <- sprintf(" (1 - Probabilities of Normal) With Threshold %s",obj@args$BayesMaxPNormal)
                  }else{
                      title <- sprintf(" (1 - Probabilities of Normal) Before Filtering")
                  }
                  infercnv::plot_cnv(infercnv_obj          = obj,
                                     #k_obs_groups         = 4,
                                     #cluster_by_groups    = cluster_by_groups,
                                     title                 = title,
                                     output_filename       = file.path(file.path(obj@args$out_dir),"infercnv.NormalProbabilities"),
                                     write_expr_matrix     = FALSE,
                                     x.center              = 0,
                                     x.range               = c(0,1)
                  )
              }
          }
)

#' Plots the probability for each cnv belonging to a specific state and the probability of
#' each cell line belonging to a specific states.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname plotProbabilities-method
#' @keywords internal
#' @noRd
setGeneric(name="plotProbabilities",
           def=function(obj)
               { standardGeneric("plotProbabilities") }
)

#' @rdname plotProbabilities-method
#' @aliases plotProbabilities
#' @noRd
setMethod(f="plotProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              if (obj@args$plotingProbs == TRUE){
                  futile.logger::flog.info(paste("Creating Plots for CNV and cell Probabilities."))
                  # Plotting
                  ## plots the probability of each cell line being a particular state
                  ## plots the probability of a cnv being a particular state

                  # Plot the probabilities of epsilons
                  ep <- function(df){df[,grepl('epsilon', colnames(df))]}
                  epsilons <- lapply(obj@combined_mcmc, function(x) ep(x))
                  ## add threshold to the plot title if given
                  if (!is.null(obj@args$BayesMaxPNormal)) {
                      file_CELLplot <- sprintf("cellProbs.%s.pdf",obj@args$BayesMaxPNormal)
                  } else{
                      file_CELLplot <- "cellProbs.pdf"
                  }
                  pdf(file = file.path(file.path(obj@args$out_dir),file_CELLplot), onefile = TRUE)
                  lapply(1:length(obj@cell_probabilities), function(i){
                      print(plot_cell_prob(as.data.frame(obj@cell_probabilities[[i]]), as.character(obj@cell_gene[[i]]$cnv_regions)))
                  })
                  dev.off()

                  ## Plot the probability of each state for a CNV
                  ## add threshold to the plot title if given
                  if (!is.null(obj@args$BayesMaxPNormal)) {
                      file_CNVplot <- sprintf("cnvProbs.%s.pdf",obj@args$BayesMaxPNormal)
                  } else{
                      file_CNVplot <- "cnvProbs.pdf"
                  }
                  pdf(file = file.path(file.path(obj@args$out_dir), file_CNVplot), onefile = TRUE)
                  lapply(1:length(obj@cell_probabilities), function(i){
                      print(plot_cnv_prob(obj@cnv_probabilities[[i]], as.character(obj@cell_gene[[i]]$cnv_regions)))
                  })
                  dev.off()
              }
          }
)


#' Return the InferCNV Object with the new adjucted CNV's
#'
#' Returns Infercnv Object
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param infercnv_obj Current inferCNV object that will be adjusted based on the results of the Bayesian Network Model.
#'
#' @return An inferCNV object
#'
#' @rdname returningInferCNV-method
#' @keywords internal

setGeneric(name = "returningInferCNV",
           def = function(obj, infercnv_obj)
               { standardGeneric("returningInferCNV") }
)
#' @rdname returningInferCNV-method
#' @aliases returningInferCNV
#' @export
#' 
#' @examples
#' # data(data)
#' # data(annots)
#' # data(genes)
#' #
#' # infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data, 
#' #                                                gene_order_file=genes,
#' #                                                annotations_file=annots,
#' #                                                ref_group_names=c("normal"))
#' # infercnv_obj <- infercnv::run(infercnv_obj,
#' #                               cutoff=1,
#' #                               out_dir="../example_output", 
#' #                               cluster_by_groups=TRUE, 
#' #                               denoise=TRUE,
#' #                               HMM=TRUE,
#' #                               num_threads=2,
#' #                               no_plot=TRUE)
#' # files <- list.files("../example_output", full.names = TRUE)
#' # HMM_obj <- readRDS(files[grep("hmm_mode-samples.infercnv_obj",files)])
#' # mcmc_obj <- infercnv::inferCNVBayesNet( infercnv_obj   = infercnv_obj,
#' #                               HMM_obj         = HMM_obj,
#' #                               file_dir        = "../example_output",
#' #                               postMcmcMethod  = "removeCNV",
#' #                               out_dir         = "../example_output",
#' #                               quietly         = TRUE,
#' #                               CORES           = 2,
#' #                               plotingProbs    = FALSE,
#' #                               diagnostics     = FALSE)
#'
#' load(HMM_obj)
#' load(mcmc_obj)
#'
#' hmm.infercnv_obj <- infercnv::returningInferCNV(mcmc_obj, HMM_obj)
#'
#'
setMethod(f = "returningInferCNV",
          signature = "MCMC_inferCNV",
          definition=function(obj, infercnv_obj) {
              NewStates <- obj@States
              infercnv_obj@expr.data <- NewStates
              return(infercnv_obj)
          }
)


#' Create Diagnostic Plots And Summaries.
#'
#' Create Diagnostic Plots And Summaries in order to determine if convergence has occured.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname mcmcDiagnosticPlots-method
#' @keywords internal
#' @noRd
setGeneric(name="mcmcDiagnosticPlots",
           def=function(obj)
           { standardGeneric("mcmcDiagnosticPlots") }
)

#' @rdname mcmcDiagnosticPlots-method
#' @aliases mcmcDiagnosticPlots
#' @noRd
setMethod(f="mcmcDiagnosticPlots",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Creating Diagnostic Plots."))
              ###########################
              # trace and denisty plots
              ###########################
              #--------------------------------------
              # trace and denisty plots for each cnv
              #--------------------------------------
              ## get the theta values
              if (obj@args$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Trace and Density Plots.")) }
              cnvProb <- function(combined_samples) {
                  thetas <- combined_samples[,grepl('theta', colnames(combined_samples))]
              }
              cnvMCMCList <- lapply(1:length(obj@mcmc), function(i){
                  lapply(obj@mcmc[[i]], cnvProb)
              })
              # trace and denisty plots
              pdf(file = file.path(file.path(obj@args$out_dir),"CNVDiagnosticPlots.pdf"), onefile = TRUE)
              lapply(1:length(cnvMCMCList), function(i){
                  plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              #---------------------------------------
              # trace and denisty plots for each cell
              #---------------------------------------
              ## get the theta values
              if (obj@args$quietly == FALSE) { futile.logger::flog.info(paste("Plotting Cell Trace and Density Plots.")) }
              cellProb <- function(samples) {
                  epsilons <- samples[,grepl('epsilon', colnames(samples))]
              }

              cellMCMCList <- lapply(1:length(obj@mcmc), function(i){
                  lapply(obj@mcmc[[i]], cellProb)
              })
              # trace and denisty plots
              pdf(file = file.path(file.path(obj@args$out_dir),"CellDiagnosticPlots.pdf"), onefile = TRUE)
              lapply(1:length(cellMCMCList), function(i){
                  plot(coda::mcmc.list(cellMCMCList[[i]]))
              })
              dev.off()


              ###########################
              # Auto Correlation Plots
              ###########################
              #---------------------------------------
              # Auto Correlation for each CNV
              #---------------------------------------
              if (obj@args$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Autocorrelation Plots.")) }
              pdf(file = file.path(file.path(obj@args$out_dir),"CNVautocorrelationPlots.pdf"), onefile = TRUE)
              lapply(1:length(cnvMCMCList), function(i){
                  autocorr.plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              ###########################
              # Gelman Plots
              ###########################
              #---------------------------------------
              # Gelman for each CNV
              #---------------------------------------
              if (obj@args$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Gelman Plots.")) }
              pdf(file = file.path(file.path(obj@args$out_dir),"CNVGelmanPlots.pdf"), onefile = TRUE)
              lapply(1:length(cellMCMCList), function(i){
                  gelman.plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              ###########################
              # Summary Tables
              ###########################
              if (obj@args$quietly == FALSE) { futile.logger::flog.info(paste("Creating CNV Statistical Summary Tables.")) }
              # Function to initialize the summary tables
              theta_table <- function(x,y,w){
                  mu<- unlist(summary(x[[1]][,w])[[1]][,1])
                  stdev<- unlist(summary(x[[1]][,w])[[1]][,2])
                  q2.5<- unlist(summary(x[[1]][,w])[[2]][,1])
                  q50<- unlist(summary(x[[1]][,w])[[2]][,3])
                  q97.5<- unlist(summary(x[[1]][,w])[[2]][,5])
                  gewek = unlist(geweke.diag(x[[1]][,w], frac1=0.1, frac2=0.5))[1:length(w)]
                  df = data.frame(mu,stdev,q2.5,q50,q97.5,gewek)
                  colnames(df) <- c('Mean','St.Dev','2.5%','50%','97.5%', "Geweke")
                  rownames(df) <- c(w)
                  #return(knitr::kable(df, caption = y))
                  return(df)
              }
              # Function to get the theta (state CNV probabilities) values
              getThetas <- function(df){ df[,grepl('theta', colnames(df))] }
              # List of statistical summary tables
              summary_table <- lapply(1:length(obj@mcmc), function(i) {
                  title <- sprintf("CNV %s Summary Table", obj@cell_gene[[i]]$cnv_regions)
                  thetas <- lapply(obj@mcmc[[i]], function(x) getThetas(x))
                  w = row.names(summary(as.mcmc(thetas))[[1]])
                  return(theta_table(coda::as.mcmc(thetas), title, w))
              })
              # Theme for the grob tables
              theme.1 <- gridExtra::ttheme_default(core = list(fg_params = list(parse=TRUE, cex = 0.5)),
                                                   colhead = list(fg_params=list(parse=TRUE, cex = 0.5)),
                                                   rowhead = list(fg_params=list(parse=TRUE, cex = 0.5)))
              # List of tables, table for each CNV
              plot_list <- lapply(1:length(summary_table), function(i) {
                  ## Create table grob object
                  table <- gridExtra::tableGrob(summary_table[[i]],rows = c("State 1","State 2","State 3","State 4","State 5","State 6"), theme = theme.1)
                  ## Create the title for the table as a seperate grob object
                  title <- sprintf("%s CNV Summary Table", obj@cell_gene[[i]]$cnv_regions)
                  title <- gridExtra::tableGrob(summary_table[[i]][1,1],rows=NULL, cols=c(title))
                  ## Combine the summary table grob and the title grob
                  tab <- gridExtra::gtable_combine(title[1,], table, along=2)
                  # Adjust the position of the title
                  tab$layout[1, c("l", "r")] <- c(7, 2)
                  tab$layout[2, c("l", "r")] <- c(7, 2)
                  return(tab)
              })
              # Combine all the tablles together as one column
              test <- gridExtra::gtable_combine(plot_list, along = 2)
              # Save the tables to a PDF document
              pdf(file = file.path(file.path(obj@args$out_dir),"CNVSummaryTablels.pdf") , paper = "a4", onefile = TRUE, height = 0, width = 0)
              print(gridExtra::marrangeGrob(grobs = test, nrow = 5, ncol = 1))
              dev.off()
          }
)



##########################
# Command line arguments #
##########################
# pargs <- optparse::OptionParser()
pargs <- argparse::ArgumentParser()
pargs$add_argument(c("-f", "--infercnv_dir"),
                   type="character",
                   action='store_true',
                   dest="file_dir",
                   metavar="File_Directory",
                   help=paste("Path to files created by inferCNV.",
                              "[Default %default][REQUIRED]"))
pargs$add_argument(c("-m", "--model"),
                              type="character",
                              action='store_true',
                              dest="model_file",
                              metavar="Model_File_Path",
                              help=paste("Path to the BUGS Model file.",
                                         "[Default %default][REQUIRED]"))
pargs$add_argument(c("-p","--parallel"),
                              type="character",
                              action='store_true',
                              dest="CORES",
                              default = NULL,
                              metavar="Number_of_Cores",
                              help=paste("Option to run parallel by specifying the number of cores to be used.",
                                         "[Default %default]"))
pargs$add_argument(c("-o","--out_dir"),
                              type="character",
                              action='store_true',
                              dest="out_dir",
                              default = NULL,
                              metavar="Output_Directory",
                              help=paste("Option to set the output directory to save the outputs.",
                                         "[Default %default]"))
pargs$add_argument(c("-M","--method"),
                              type="character",
                              action='store_true',
                              dest="postMcmcMethod",
                              default = NULL,
                              metavar="Posterior_MCMC_Method",
                              help=paste("What actions to take after finishing the MCMC.",
                                         "[Default %default]"))
pargs$add_argument(c("-x","--plot"),
                              type="logical",
                              action='store_true',
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
run_gibb_sampling <- function(gene_exp,
                              MCMC_inferCNV_obj
                              ){
    if (is.null(ncol(gene_exp))){
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
    stats::update(model, 200, progress.bar=ifelse(MCMC_inferCNV_obj@args$quietly,"none","text"))
    # run the rjags model
    ## set the parameters to return from sampling
    parameters <- c('theta', 'epsilon')
    samples <- rjags::coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(MCMC_inferCNV_obj@args$quietly,"none","text"))
    return(samples)
}

# Function to plot the probability for each cell line of being in a particular state
plot_cell_prob <- function(df, title){
    df$mag = c(1:6)
    long_data <- reshape::melt(df, id = "mag")
    long_data$mag <- as.factor(long_data$mag)
    ggplot2::ggplot(long_data, ggplot2::aes_string(x = 'variable', y = 'value', fill = 'mag'))+
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
        ggplot2::labs(fill = "States")+
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
plot_cnv_prob <- function(df,title){
    colnames(df) <- c(1:6)
    df <- melt(df)
    colnames(df) <- c("row", "State", "Probability")
    states <- as.factor(df$State)
    ggplot2::ggplot(data = df, ggplot2::aes_string(y = 'Probability', x= 'State', fill = 'states')) +
        ggplot2::geom_boxplot()+
        ggplot2::labs(title = title) +
        ggplot2::theme(plot.title = element_text(hjust = 0.5))
}


###############################################
# Main Function to run Bayesion Network Model #
###############################################
#' @title inferCNVBayesNet: Run Bayesian Network Mixture Model To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by
#' inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster and each individual
#' cell line in the CNV.
#'
#' @param file_dir Location of the directory of the inferCNV outputs.
#' @param infercnv_obj InferCNV object.
#' @param HMM_obj InferCNV object with HMM states in expression data.
#' @param model_file Path to the BUGS Model file.
#' @param CORES Option to run parallel by specifying the number of cores to be used. (Default: 1)
#' @param out_dir (string) Path to where the output file should be saved to.
#' @param postMcmcMethod What actions to take after finishing the MCMC.
#' @param plotingProbs Option for adding plots of Cell and CNV probabilities. (Default: TRUE)
#' @param quietly Option to print descriptions along each step. (Default: TRUE)
#' @param diagnostics Option to plot Diagnostic plots and tables. (Default: FALSE)
#'
#' @return Returns a MCMC_inferCNV_obj and posterior probability of being in one of six Copy Number Variation states
#' (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by inferCNV's HMM.
#'
#' @export
#'
#' @examples
#' data(data)
#' data(annots)
#' data(genes)
#'
#' infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data, 
#'                                                gene_order_file=genes,
#'                                                annotations_file=annots,
#'                                                ref_group_names=c("normal"))
#' infercnv_obj <- infercnv::run(infercnv_obj,
#'                               cutoff=1,
#'                               out_dir="../example_output", 
#'                               cluster_by_groups=TRUE, 
#'                               denoise=TRUE,
#'                               HMM=TRUE,
#'                               num_threads=2,
#'                               no_plot=TRUE)
#'
#' files <- list.files("../example_output", full.names = TRUE)
#' HMM_obj <- readRDS(files[grep("hmm_mode-samples.infercnv_obj",files)])
#' mcmc_obj <- infercnv::inferCNVBayesNet( infercnv_obj   = infercnv_obj,
#'                               HMM_obj         = HMM_obj,
#'                               file_dir        = "../example_output",
#'                               postMcmcMethod  = "removeCNV",
#'                               out_dir         = "../example_output",
#'                               quietly         = TRUE,
#'                               CORES           = 2,
#'                               plotingProbs    = FALSE,
#'                               diagnostics     = FALSE)

inferCNVBayesNet <- function(
                              file_dir,
                              infercnv_obj,
                              HMM_obj,
                              out_dir,
                              model_file      = system.file("BUGS_Mixture_Model",package = "infercnv"),
                              CORES           = 1,
                              postMcmcMethod  = NULL,
                              plotingProbs    = TRUE,
                              quietly         = TRUE,
                              diagnostics     = FALSE) {

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
    if (!(CORES == 1)){
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
                        "quietly" = quietly,
                        "BayesMaxPNormal" = 0)
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
    ## Run Gibbs sampling and time the process
    start_time <- Sys.time()
    MCMC_inferCNV_obj <- runMCMC(MCMC_inferCNV_obj)
    end_time <- Sys.time()
    futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))

    ## Save the MCMC.infercnv_object as an RDS
    saveRDS(MCMC_inferCNV_obj, file = file.path(MCMC_inferCNV_obj@args$out_dir, "MCMC_inferCNV_obj.rds"))

    ########
    # Plot #
    ########
    if (diagnostics == TRUE){
        mcmcDiagnosticPlots(MCMC_inferCNV_obj)
    }
    postProbNormal(MCMC_inferCNV_obj,
               PNormal = NULL)

    return(MCMC_inferCNV_obj)
}

#############################################################
# Function to modify CNV's identified base on probabilities #
#############################################################
#' @title filterHighPNormals: Filter the HMM identified CNV's by the CNV's posterior probability
#' of belonging to a normal state.
#'
#' @description The following function will filter the HMM identified CNV's by the CNV's posterior
#' probability of belonging to a normal state identified by the function inferCNVBayesNet(). Will filter
#' CNV's based on a user desired threshold probability. Any CNV with a probability of being normal above
#' the threshold will be removed.
#'
#' @param MCMC_inferCNV_obj MCMC infernCNV object.
#' @param BayesMaxPNormal Option to filter CNV or cell lines by some probability threshold.
#'
#' @return Returns a MCMC_inferCNV_obj With removed CNV's.
#'
#' @export
#' 
#' @examples
#' # data(data)
#' # data(annots)
#' # data(genes)
#'
#' # infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=data, 
#' #                                                gene_order_file=genes,
#' #                                                annotations_file=annots,
#' #                                                ref_group_names=c("normal"))
#' # infercnv_obj <- infercnv::run(infercnv_obj,
#' #                               cutoff=1,
#' #                               out_dir="../example_output", 
#' #                               cluster_by_groups=TRUE, 
#' #                               denoise=TRUE,
#' #                               HMM=TRUE,
#' #                               num_threads=2,
#' #                               no_plot=TRUE)
#' # files <- list.files("../example_output", full.names = TRUE)
#' # HMM_obj <- readRDS(files[grep("hmm_mode-samples.infercnv_obj",files)])
#' # mcmc_obj <- infercnv::inferCNVBayesNet( infercnv_obj   = infercnv_obj,
#' #                               HMM_obj         = HMM_obj,
#' #                               file_dir        = "../example_output",
#' #                               postMcmcMethod  = "removeCNV",
#' #                               out_dir         = "../example_output",
#' #                               quietly         = TRUE,
#' #                               CORES           = 2,
#' #                               plotingProbs    = FALSE,
#' #                               diagnostics     = FALSE)
#' # hmm.infercnv_obj <- infercnv::returningInferCNV(mcmc_obj, HMM_obj)
#'
#' load(mcmc_obj)
#'
#' mcmc_obj <- infercnv::filterHighPNormals( MCMC_inferCNV_obj = mcmc_obj, 
#'                               BayesMaxPNormal   = 0.5)
#'
#'
#'
#'
#'

filterHighPNormals <- function( MCMC_inferCNV_obj,
                                BayesMaxPNormal) {
    MCMC_inferCNV_obj <- setBayesMaxPNormal( obj             = MCMC_inferCNV_obj,
                                             BayesMaxPNormal = BayesMaxPNormal )
    ## Either Remove CNV's based on CNV posterier probabilities ("removeCNV")
    ## or remove cell lines based on cell line posterior probabilities ("removeCells")
    if(!(is.null(MCMC_inferCNV_obj@args$postMcmcMethod))){
        if(MCMC_inferCNV_obj@args$postMcmcMethod == "removeCNV"){
            MCMC_inferCNV_obj <- removeCNV(MCMC_inferCNV_obj)
        } else {
            MCMC_inferCNV_obj <- removeCells(MCMC_inferCNV_obj)
        }
    }

    ## Plot the resuls

    plotProbabilities(MCMC_inferCNV_obj)
    postProbNormal(MCMC_inferCNV_obj,
                   PNormal = TRUE)

    return(MCMC_inferCNV_obj)
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
