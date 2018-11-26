
get_spike_dists <- function(hspike_obj) {
    
    if (is.null(hspike_obj)) {
        flog.error("get_spike_dists(hspike_obj): Error, hspike obj is null")
        stop("Error")
    }
    if (! is.null(hspike_obj@.hspike)) {
        flog.error("get_spike_dists() should have an hspike obj as param, and this doesn't look like an hspike")
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


## Based on number of cells in a clade
get_hspike_cnv_mean_sd_trend_by_num_cells_fit <- function(hspike_obj, plot=F) {
    
    gene_expr_by_cnv <- .get_gene_expr_by_cnv(hspike_obj)
    cnv_level_to_mean_sd = list()
    
    for (cnv_level in names(gene_expr_by_cnv) ) {
        expr_vals = gene_expr_by_cnv[[ cnv_level ]]
        nrounds = 100
        
        sds = c()
        for (ncells in seq_len(100)) {
            means = c()
            
            for(i in 1:nrounds) {
                vals = sample(expr_vals, size=ncells, replace=T)
                m_val = mean(vals)
                means = c(means,  m_val)
            }
            sds = c(sds, sd(means))
        }
        cnv_level_to_mean_sd[[ cnv_level ]] <- sds
    }
    
    if (plot) {
        df = do.call(rbind, lapply(names(cnv_level_to_mean_sd), function(cnv_level) {
            sd=cnv_level_to_mean_sd[[ cnv_level ]]
            data.frame(cnv=cnv_level, num_cells=1:length(sd), sd=sd)
        }))

        p  = df %>% ggplot(aes(x=log(num_cells),y=log(sd), color=cnv)) + geom_point()
        pdf("hspike_obs_mean_var_trend.pdf")
        flog.info("plotting: hspike_obs_mean_var_trend.pdf")
        plot(p)
        dev.off()
    }

    ## fit linear model
    cnv_level_to_mean_sd_fit = list()
    for (cnv_level in names(cnv_level_to_mean_sd) ) {
        sd_vals=cnv_level_to_mean_sd[[ cnv_level ]]
        num_cells = 1:length(sd_vals)

        flog.info(sprintf("fitting num cells vs. variance for cnv level: %s", cnv_level))
        fit = lm(log(sd_vals) ~ log(num_cells)) #note, hbadger does something similar, but not for the hmm cnv state levels
        cnv_level_to_mean_sd_fit[[ cnv_level ]] = fit
    }
    
    return(cnv_level_to_mean_sd_fit)
    
}







.get_HMM <- function(cnv_mean_sd, t) {

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



    HMM_info = list(state_transitions=state_transitions,
                    delta=delta,
                    state_emission_params=state_emission_params)
    
    return(HMM_info)
}

    

predict_CNV_via_HMM_on_indiv_cells  <- function(infercnv_obj, cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike), t=1e-6) {
    
    flog.info("predict_CNV_via_HMM_on_indiv_cells()")
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
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
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       HMM_info[['state_emission_params']])
            
            hmm_trace <- HiddenMarkov::Viterbi(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}
            
predict_CNV_via_HMM_on_tumor_subclusters  <- function(infercnv_obj,
                                                      cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                      cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                      t=1e-6) {


    flog.info("predict_CNV_via_HMM_on_tumor_subclusters()")
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
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
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx])
            
            num_cells = length(tumor_subcluster_cells_idx)

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
                        
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            hmm_trace <- HiddenMarkov::Viterbi(hmm)
            
            hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}


.get_state_emission_params <- function(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit, plot=FALSE) {
        
    for (cnv_level in names(cnv_mean_sd)) {
        fit <- cnv_level_to_mean_sd_fit[[cnv_level]]
        sd <- exp(predict(fit, newdata=data.frame(num_cells=num_cells))[[1]])
        cnv_mean_sd[[cnv_level]]$sd <- sd
    }
    
    if (plot) {
        .plot_cnv_mean_sd_for_num_cells(num_cells, cnv_mean_sd)
    }

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
    
    
    return(state_emission_params)
}

.plot_cnv_mean_sd_for_num_cells <- function(num_cells, cnv_mean_sd) {

    pdf(sprintf("state_emissions.%d-cells.pdf", num_cells))
    
    p = ggplot(data.frame(x=c(0,3)), aes(x)) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd),
                      aes(colour='cnv:0')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd),
                      aes(colour='cnv:0.5')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd),
                      aes(colour='cnv:1')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd),
                      aes(colour='cnv:1.5')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd),
                      aes(colour='cnv:2')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd),
                      aes(colour='cnv:3')) +
        scale_colour_manual("Function",
                            values=c("blue", "red", "green", "magenta", "orange", "cyan"),
                            breaks=c("cnv:0","cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3"))
    
    plot(p)

    dev.off()

}


.compare_obs_vs_fit_sc_cnv_sd <- function(infercnv_obj,
                                          cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                          cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)) {
    
    cnv_levels = names(cnv_mean_sd)

    df = do.call(rbind, lapply(cnv_levels, function(cnv_level) {
        sd_obs=cnv_mean_sd[[cnv_level]]$sd
        sd_pred=exp(predict(cnv_level_to_mean_sd_fit[[cnv_level]], newdata=data.frame(num_cells=1)))
        return(data.frame(cnv=cnv_level, sd_obs=sd_obs, sd_pred=sd_pred))

    }))
    
    return(df)
}




    
