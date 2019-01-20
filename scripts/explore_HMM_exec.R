#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)
library(futile.logger)
library(HiddenMarkov)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

cnv_mean_sd=infercnv:::get_spike_dists(infercnv_obj@.hspike)
cnv_level_to_mean_sd_fit=infercnv:::get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)
t=1/6
p_val=0.05
hclust_method='ward.D2'


flog.info(sprintf("predict_CNV_via_HMM_on_tumor_subclusters(p_val=%g)", p_val))
HMM_info  <- infercnv:::.get_HMM(cnv_mean_sd, t)
chrs = unique(infercnv_obj@gene_order$chr)
expr.data = infercnv_obj@expr.data
gene_order = infercnv_obj@gene_order
hmm.data = expr.data
hmm.data[,] = -1 #init to invalid state

tumor_subclusters <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=F)


##########################################
#chrs = c('chr1')
##########################################


##############################################
## From HiddenMarkovPackage
getj <- function (x, j)  {
    if (is.null(x)) 
        return(NULL)
    n <- length(x)
    for (i in 1:n) x[[i]] <- x[[i]][j]
    return(x)
}


Viterbi.dthmm <- function (object, ...){
    x <- object$x
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    n <- length(x)
    m <- nrow(object$Pi) # transition matrix
    nu <- matrix(NA, nrow = n, ncol = m)  # scoring matrix
    y <- rep(NA, n) # final trace
    pseudocount = 1e-20
    
    emissions <- matrix(NA, nrow = n, ncol = m) 
    
    ## init first row
    emission <- pnorm(abs(x[1]-object$pm$mean)/object$pm$sd, log=T, lower.tail=F)
    emission <- 1 / (-1 * emission)
    emission <- emission / sum(emission)
    
    emissions[1,] <- log(emission)
    
    nu[1, ] <- log(object$delta) + # start probabilities
        emissions[1,]
    
    logPi <- log(object$Pi) # convert transition matrix to log(p)
    
    for (i in 2:n) {
        
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        #nu[i, ] <- apply(matrixnu + logPi, 2, max) +
        #              dfunc(x=x[i], object$pm, getj(object$pn, i),
        #                    log=TRUE)

        
        #emission <- dfunc(x=x[i], object$pm, getj(object$pn, i), log=T)
        ## normalize emission p-values
        ## first add pseudcounts
        #missions[i, ] <- emissions[i, ] + pseudocount
        #emissions[i, ] <- emissions[i, ] / sum(emissions[i, ]) 
 
        #emissions[i, ] <- log(emissions[i, ])
                

        emission <- pnorm(abs(x[i]-object$pm$mean)/object$pm$sd, log=T, lower.tail=F)
        emission <- 1 / (-1 * emission)
        emission <- emission / sum(emission)
        
        emissions[i, ] <- log(emission)
        
        nu[i, ] <- apply(matrixnu + logPi, 1, max) + emissions[i, ] 
                
    }
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")


    ## traceback
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1))
        y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])

    return(y)
}


##########################################


for (chr in chrs) {
    print(chr)
    chr_gene_idx = which(gene_order$chr == chr)
    
    ## run through each cell for this chromosome:
    for (tumor_subcluster_name in names(tumor_subclusters)) {
        print(tumor_subcluster_name)
        tumor_subcluster_cells_idx <- tumor_subclusters[[tumor_subcluster_name]]
                
        gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx,drop=F])
        
        num_cells = length(tumor_subcluster_cells_idx)
        
        state_emission_params <- infercnv:::.get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
        
        hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                   HMM_info[['state_transitions']],
                                   HMM_info[['delta']],
                                   "norm",
                                   state_emission_params)

        hmm_trace <- Viterbi.dthmm(hmm)
        
        print(hmm_trace)
        
        hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <- hmm_trace
    }
}




