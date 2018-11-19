
get_spike_dists <- function(hspike_obj) {

    if (is.null(hspike_obj)) {
        flog.error("get_spike_dists(hspike_obj): Error, hspike obj is null")
        stop("Error")
    }
    
    gene_expr_by_cnv = .get_gene_expr_by_cnv(hspike_obj)
        
    cnv_mean_sd = .get_gene_expr_mean_sd_by_cnv(gene_expr_by_cnv)

    return(cnv_mean_sd)
    
}



.get_gene_expr_by_cnv <- function(hspike_obj) {
    
    chr_info = .get_hspike_chr_info()

    spike_cell_idx = unlist(hspike_obj@observation_grouped_cell_indices)

    spike.expr.data = hspike_obj@expr.data[, spike_cell_idx]

    gene_expr_by_cnv = list()
    for (info in chr_info) {
        chr_name = info$name
        cnv = sprintf("cnv:%g", info$cnv)
        chr_gene_idx = which(hspike_obj@gene_order$chr == chr_name)
        gene.expr = c(spike.expr.data[chr_gene_idx, ])
        
        if (cnv %in% names(gene_expr_by_cnv)) {
            gene_expr_by_cnv[[cnv]] = c(gene_expr_by_cnv[[cnv]], gene.expr)
        }  else {
            gene_expr_by_cnv[[cnv]] = gene.expr
        }
    }
    
    return(gene_expr_by_cnv)
}


.get_gene_expr_mean_sd_by_cnv <- function(gene_expr_by_cnv) {

    cnv_mean_sd = list()
    
    for (cnv_level in names(gene_expr_by_cnv) ) {
        gene_expr = gene_expr_by_cnv[[ cnv_level ]]

        gene_expr_mean = mean(gene_expr)
        gene_expr_sd = sd(gene_expr)

        cnv_mean_sd[[ cnv_level ]] = list(mean=gene_expr_mean, sd=gene_expr_sd)
    }
    
    return(cnv_mean_sd)
    
}

.plot_gene_expr_by_cnv <- function(gene_expr_by_cnv, cnv_mean_sd) {
    
    df  = do.call(rbind, lapply(names(gene_expr_by_cnv), function(x) { data.frame(cnv=x, expr=gene_expr_by_cnv[[x]]) }))

    p = df %>% ggplot(aes(expr,  fill=cnv, colour=cnv))  +  geom_density(alpha=0.1)

    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd)) 

    return(p)
}



predict_CNV_via_HMM  <- function(infercnv_obj, cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike), t=1e-6) {
    
    flog.info("predict_CNV_via_HMM()")
    
    ##                      states: 0,   0.5,    1,     1.5,     2,    3
    state_transitions = matrix( c(1-5*t,  t,     t,      t,      t,    t,
                                  t,     1-5*t,  t,      t,      t,    t,
                                  t,      t,     1-5*t,  t,      t,    t,
                                  t,      t,     t,    1-5*t,    t,    t,
                                  t,      t,     t,      t,    1-5*t,  t,
                                  t,      t,     t,      t,      t,  1-5*t),
                               byrow=T,
                               nrow=6)

    delta=c(t,      t,     t,    1-5*t,    t,    t) # more likely normal,

    state_emission_params = list(mean=c(cnv_mean_sd[["cnv:0.01"]]$mean,
                                        cnv_mean_sd[["cnv:0.5"]]$mean,
                                        cnv_mean_sd[["cnv:1"]]$mean,
                                        cnv_mean_sd[["cnv:1.5"]]$mean,
                                        cnv_mean_sd[["cnv:2"]]$mean,
                                        cnv_mean_sd[["cnv:3"]]$mean),

                                 sd=c(cnv_mean_sd[["cnv:0.01"]]$sd,
                                      cnv_mean_sd[["cnv:0.5"]]$sd,
                                      cnv_mean_sd[["cnv:1"]]$sd,
                                      cnv_mean_sd[["cnv:1.5"]]$sd,
                                      cnv_mean_sd[["cnv:2"]]$sd,
                                      cnv_mean_sd[["cnv:3"]]$sd) )
    
    chrs = unique(infercnv_obj@gene_order$chr)


    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(seq_len(ncol(expr.data)), function(cell_idx) {
            
            gene_expr_vals = as.vector(expr.data[chr_gene_idx,cell_idx])
            
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       state_transitions,
                                       delta,
                                       "norm",
                                       state_emission_params)
            
            hmm_trace <- HiddenMarkov::Viterbi(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}
            
