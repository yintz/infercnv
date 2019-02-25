

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
           {
               standardGeneric("MeanSD")
           }
)

#' @rdname MeanSD-method
#' @aliases MeanSD
#' 
setMethod(f="MeanSD",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              gene_expr_by_cnv = infercnv:::.get_gene_expr_by_cnv(obj@.hspike)
              cnv_mean_sd = infercnv:::.get_gene_expr_mean_sd_by_cnv(gene_expr_by_cnv)
              cnv_sd <- cbind(lapply(cnv_mean_sd,function(x){x$sd}))
              cnv_mean <- cbind(lapply(cnv_mean_sd,function(x){x$mean}))
              ## Sort so in order of {x0,...,x1,..,x3} and get into a vector format
              obj@mu <- unlist(cbind(cnv_mean[sort(row.names(cnv_mean)),]))
              obj@sig <-  unlist(cbind(cnv_sd[sort(row.names(cnv_sd)),]))
              obj@sig <- 1/(obj@sig^2)
              print(paste("Means: ", obj@mu, collapse = ""))
              print(paste("Sig: ", obj@sig, collapse = ""))
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
           {
               standardGeneric("getGenesCells")
           }
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
           def=function(obj, args_parsed)
           {
               standardGeneric("initializeObject")
           }
)

#' @rdname initializeObject-method
#' @aliases initializeObject
#' 
setMethod(f="initializeObject",
          signature="MCMC_inferCNV",
          definition=function(obj, args_parsed)
          {
              futile.logger::flog.info(paste("Loading InferCNV Object."))
              files <- list.files(args_parsed$file_dir, full.names = TRUE)
              infercnv_obj_PATH <- files[grep(files, pattern = "preliminary.infercnv_obj")] #"preliminary.infercnv_obj")]
              infercnv_obj = readRDS(infercnv_obj_PATH)
              # Validate the inferCNV Object 
              infercnv:::validate_infercnv_obj(infercnv_obj)
              
              ## create the S4 object
              obj <- MCMC_inferCNV(infercnv_obj)
              ## add the command line arguments 
              obj@args <- args_parsed
              
              ## Load the files for cnv predictions 
              cell_groups_PATH <- files[grep(files, pattern = "12_HMM_preds.cell_groupings")]
              pred_cnv_genes_PATH <- files[grep(files, pattern = "12_HMM_preds.pred_cnv_genes.dat")]
              # cell_groups_df <- read.table(cell_groups_PATH, header = T)
              cell_groups_df <- read.csv(cell_groups_PATH, sep = "\t", header = T, check.names = FALSE)
              pred_cnv_genes_df <- read.table(pred_cnv_genes_PATH, header = T, check.names = FALSE)
              
              ## Load Mixture Model File
              futile.logger::flog.info(paste("Loading BUGS Model."))
              obj@bugs_model <- readChar(obj@args$model_file,file.info(obj@args$model_file)$size)
              
              # cnv region id's 
              obj@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
              futile.logger::flog.info(paste("Total CNV's: ", length(obj@cnv_regions)))
              
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
           def=function(obj)
           {
               standardGeneric("getStates")
           }
)

#' @rdname getStates-method
#' @aliases getStates
#' 
setMethod(f="getStates",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Loading InferCNV Object."))
              files <- list.files(obj@args$file_dir, full.names = TRUE)
              # Add the HMM defined states
              HMM_obj_PATH <- files[grep(files, pattern = "12_HMM_pred.infercnv_obj")]
              HMM_obj = readRDS(HMM_obj_PATH)
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
           {
               standardGeneric("getProbabilities")
           }
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
           {
               standardGeneric("withParallel")
           }
)

#' @rdname withParallel-method
#' @aliases withParallel
#' 
setMethod(f="withParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Running Sampling Using Parallel with ",CORES,"Cores"))
              obj@mcmc <- parallel::mclapply(1:length(obj@cell_gene), function(i){ 
                  futile.logger::flog.info(paste("Sampleing Number: ", i))
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- group_id[obj@cell_gene[[i]]$Cells] # subset the tumor ids for the cells wanted
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

#' Run simulations NOT in Parallel
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
           {
               standardGeneric("nonParallel")
           }
)

#' @rdname nonParallel-method
#' @aliases nonParallel
#' 
setMethod(f="nonParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
            futile.logger::flog.info(paste("Not Running MCMC in Parallel Mode."))
            obj@mcmc <- lapply(1:length(obj@cell_gene), function(i){
                futile.logger::flog.info(paste("Sample Number: ", i))
                if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                    tumor_grouping <- obj@group_id[obj@cell_gene[[i]]$Cells] # subset the tumor ids for the cells wanted
                    gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes,obj@cell_gene[[i]]$Cells]
                    return(run_gibb_sampling(gene_exp, obj))
                } else {
                    return(list(NULL))
                }
            })
            return(obj)
          }
)

#' Run simulations and remove cnv's that are predicted to be normal 
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
           {
               standardGeneric("removeCNV")
           }
)

#' @rdname removeCNV-method
#' @aliases removeCNV
#' 
setMethod(f="removeCNV",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              if (any(cnv_means[3,] >.5)){
                  remove_cnv <- which(cnv_means[3,] >.5)
                  futile.logger::flog.info(paste("Removing ",length(remove_cnv), " CNV's identified by the HMM."))
                  print("CNV's being removed have the following posterior probabilities of being a normal state: ")
                  
                  lapply(remove_cnv, function(i) {
                      print(paste("Genes: ", length(obj@cell_gene[[i]]$Genes), " Cells: ", length(obj@cell_gene[[i]]$Cells) ))
                      ## change the states to normal states
                      obj@States[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- 3
                  })
                  obj@cell_gene <- obj@cell_gene[-remove_cnv]
              }
              saveRDS(obj, file = "MCMC_inferCNV_obj_removed_cnvs")
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
           {
               standardGeneric("removeCells")
           }
)

#' @rdname removeCells-method
#' @aliases removeCells
#' 
setMethod(f="removeCells",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              if (any(do.call(cbind, obj@cell_probabilities)[3,] >.5)){
                  lapply(1:length(obj@cell_probabilities), function(i) {
                      idx <- which(obj@cell_probabilities[[i]][3,] >.5)
                      if(length(idx) > 0){
                          ## change the states to normal states
                          obj@States[ obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells[ idx ] ] <<- 3
                          ## remove these cells from the cnv
                          obj@cell_gene[[i]]$Cells <<- obj@cell_gene[[i]]$Cells[- idx]
                      }
                  })
                  saveRDS(obj, file = "MCMC_inferCNV_obj_removed_cells")
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
           {
               standardGeneric("runMCMC")
           }
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
              saveRDS(obj, file = "MCMC_inferCNV_obj")
              
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
           {
               standardGeneric("postProbNormal")
           }
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
              infercnv::plot_cnv(infercnv_obj=obj,
                                 #k_obs_groups=4,
                                 #cluster_by_groups=cluster_by_groups,
                                 out_dir=obj@args$output_dir,
                                 title=sprintf("NormalProbabilities"),
                                 output_filename=sprintf("infercnv.NormalProbabilities"),
                                 write_expr_matrix=FALSE,
                                 x.center= 0,
                                 x.range=c(0,1)
              )
          }
)

#' Main plotting function, plots the probability for each cnv belonging to a specific state and the probability of 
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
           {
               standardGeneric("plotProbabilities")
           }
)

#' @rdname plotProbabilities-method
#' @aliases plotProbabilities
#' 
setMethod(f="plotProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              # Main plotting function 
              ## plots the probability of each cell line being a particular state
              ## plots the probability of a cnv being a particular state 
              
              # Plot the probabilities of epsilons 
              ep <- function(df){df[,grepl('epsilon', colnames(df))]}
              epsilons <- lapply(obj@mcmc, function(x) ep(x))
              
              pdf(file = './cellProbs.pdf', onefile = TRUE)
              lapply(1:length(obj@cell_probabilities), function(i){ 
                  print(plot_cell_prob(as.data.frame(obj@cell_probabilities[[i]]), as.character(obj@cell_gene[[i]]$cnv_regions)))
              })
              dev.off()
              
              ## Plot the probability of each state for a CNV 
              pdf(file = './cnvProbs.pdf', onefile = TRUE)
              lapply(obj@cnv_probabilities,function(x){
                  print(plot_cnv_prob(x))
              })
              dev.off()
          }
)



