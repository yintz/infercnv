

## Based on number of cells in a clade
.ZHMM_get_normal_sd_trend_by_num_cells_fit <- function(infercnv_obj, z_p_val=0.05, plot=TRUE) {
    
    if (!has_reference_cells(infercnv_obj)) {
        stop("Error, cannot tune parameters without reference 'normal' cells defined")
    }


    normal_samples = infercnv_obj@reference_grouped_cell_indices
    
    normal_expr_vals <- infercnv_obj@expr.data[,unlist(normal_samples)]
    
    mu = mean(normal_expr_vals)
    sigma = sd(normal_expr_vals)
    nrounds = 100
    sds = c()
    ngenes = nrow(normal_expr_vals)

    num_normal_samples = length(normal_samples)
    message("-got ", num_normal_samples, " samples") 
    for (ncells in seq_len(100)) {
        means = c()
        
        for(i in 1:nrounds) {
            ## pick a random gene
            rand.gene = sample(1:ngenes, size=1)
            
            ## pick a random normal cell type
            rand.sample = sample(1:num_normal_samples, size=1)
            #message("rand.sample: " , rand.sample)
            
            vals = sample(infercnv_obj@expr.data[rand.gene, normal_samples[[rand.sample]] ], size=ncells, replace=T)
            m_val = mean(vals)
            means = c(means,  m_val)
        }
        sds = c(sds, sd(means))
    }
    
    ## fit linear model
    num_cells = 1:length(sds)
    fit = lm(log(sds) ~ log(num_cells)) #note, hbadger does something similar, but not for the hmm cnv state levels
    
    if (plot) {
        plot(log(num_cells), log(sds), main='log(sd) vs. log(num_cells)')
    }
    
    ## store mean_delta for the single gene for convenience sake
    mean_delta = determine_mean_delta_via_Z(sigma, p=z_p_val)    

    message("mean_delta: ", mean_delta, ", at sigma: ", sigma, ", and pval: ", z_p_val)
    
    normal_sd_trend = list(mu=mu,
                           sigma=sigma,
                           fit=fit,
                           mean_delta=mean_delta)
    
    return(normal_sd_trend)
    
}


.ZHMM_get_HMM <- function(normal_sd_trend, num_cells, t, z_p_val=0.05) {

    ## Here we do something very similar to HoneyBadger
    ## which is to estimate the mean/var for the CNV states
    ## based on the KS statistic and the expression values
    ## for the normal cells.
    
    
    ##                      states: 0.5,    1,     1 .5
    state_transitions = matrix( c(1-5*t,  t,     t,
                                  t,     1-5*t,  t,
                                  t,      t,     1-5*t),
                               byrow=T,
                               nrow=3)
    
    delta=c(t,     1-5*t,    t) # more likely normal,
        
    mu = normal_sd_trend$mu # normal cell mean
    
    if (num_cells == 1) {
        sigma = normal_sd_trend$sigma
        mean_delta = normal_sd_trend$mean_delta
    } else {
        ## use the var vs. num cells trend
        sigma <- exp(predict(normal_sd_trend$fit,
                             newdata=data.frame(num_cells=num_cells))[[1]])  

        ## use t-test to determine adjusted mean value
        
     
        #mean_delta <- determine_mean_delta_via_Z(sigma, p=z_p_val)

        mean_delta = normal_sd_trend$mean_delta
    }
    
    
    state_emission_params = list(mean=c(
                                     mu - mean_delta, # state 0.5 = deletion
                                     mu, # state 1 = neutral
                                     mu + mean_delta), # state 1.5 = amplification
                                 
                                 sd=c(sigma,
                                      sigma,
                                      sigma)  # shared variance as in HB
                                 )
    
    HMM_info = list(state_transitions=state_transitions,
                    delta=delta,
                    state_emission_params=state_emission_params)
    
    return(HMM_info)
}


ZHMM_predict_CNV_via_HMM_on_indiv_cells  <- function(infercnv_obj,
                                                     z_p_val=0.05,
                                                     normal_sd_trend=.ZHMM_get_normal_sd_trend_by_num_cells_fit(infercnv_obj, z_p_val),
                                                     t=1e-6) {
    
    flog.info("predict_CNV_via_HMM_on_indiv_cells()")
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    HMM_info  <- .ZHMM_get_HMM(normal_sd_trend, num_cells = 1, t=t, z_p_val=z_p_val)

    message(HMM_info)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(seq_len(ncol(expr.data)), function(cell_idx) {
            
            gene_expr_vals = as.vector(expr.data[chr_gene_idx,cell_idx])

            hmm <- HiddenMarkov::dthmm(x=gene_expr_vals,
                                       Pi=HMM_info[['state_transitions']],
                                       delta=HMM_info[['delta']],
                                       distn="norm",
                                       pm=HMM_info[['state_emission_params']])
            
            
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}
            
ZHMM_predict_CNV_via_HMM_on_tumor_subclusters  <- function(infercnv_obj,
                                                           z_p_val=0.05,
                                                           normal_sd_trend=.ZHMM_get_normal_sd_trend_by_num_cells_fit(infercnv_obj, z_p_val),
                                                           t=1e-6
                                                           ) {
    
    
    flog.info(sprintf("ZHMM_predict_CNV_via_HMM_on_tumor_subclusters(z_p_val=%g)", z_p_val))

    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("No subclusters defined, so instead running on whole samples")
        return(ZHMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, z_p_val, normal_sd_trend, t));
    }
    
    
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    tumor_subclusters <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=F)
    
    ## add the normals, so they get predictions too:
    tumor_subclusters <- c(tumor_subclusters, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_subclusters, function(tumor_subcluster_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx,drop=F])
            
            num_cells = length(tumor_subcluster_cells_idx)

            HMM_info  <- .ZHMM_get_HMM(normal_sd_trend, num_cells=num_cells, t=t, z_p_val=z_p_val)
                        
                                    
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       HMM_info[['state_emission_params']])
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}



ZHMM_predict_CNV_via_HMM_on_whole_tumor_samples  <- function(infercnv_obj,
                                                        z_p_val=0.05,
                                                        normal_sd_trend=.ZHMM_get_normal_sd_trend_by_num_cells_fit(infercnv_obj, z_p_val),
                                                        t=1e-6
                                                        ) {
    
    
    flog.info("predict_CNV_via_HMM_on_whole_tumor_samples")
    
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    ## add the normals, so they get predictions too:
    tumor_samples <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_samples, function(tumor_sample_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_sample_cells_idx,drop=F])
            
            num_cells = length(tumor_sample_cells_idx)
            
            HMM_info  <- .ZHMM_get_HMM(normal_sd_trend, num_cells=num_cells, t=t, z_p_val=z_p_val)
            
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       HMM_info[['state_emission_params']])
            
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_sample_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}

ZHMM_assign_HMM_states_to_proxy_expr_vals <- function(infercnv_obj) {
    
    expr.data = infercnv_obj@expr.data
    
    expr.data[expr.data == 1] <- 0.5
    expr.data[expr.data == 2] <- 1
    expr.data[expr.data == 3] <- 1.5
    
    infercnv_obj@expr.data <- expr.data

    return(infercnv_obj)

}


determine_mean_delta_via_Z <- function(sigma, p) {
    
    ## want tails of the distribution to minimially overlap at the p-value
    
    delta = abs(qnorm(p=p, mean=0, sd=sigma)) 

    flog.info(sprintf("determine mean delta (sigma: %g, p=%g) -> %g", sigma, p, delta))
    
    return(2 * delta)
    
}

