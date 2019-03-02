
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
    
    chr_info = .get_hspike_chr_info(1,1) # dummy values for now

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

    delta=c(t,      t,     1-5*t,    t,    t,   t) # more likely normal,
    
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
            
            hmm <- HiddenMarkov::dthmm(x=gene_expr_vals,
                                       Pi=HMM_info[['state_transitions']],
                                       delta=HMM_info[['delta']],
                                       distn="norm",
                                       pm=HMM_info[['state_emission_params']])
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}
            
predict_CNV_via_HMM_on_tumor_subclusters  <- function(infercnv_obj,
                                                      cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                      cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                      t=1e-6
                                                      ) {
    
    
    flog.info(sprintf("predict_CNV_via_HMM_on_tumor_subclusters"))
    
    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("No subclusters defined, so instead running on whole samples")
        return(predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, cnv_mean_sd, cnv_level_to_mean_sd_fit, t));
    }
    
    
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
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx,drop=F])
            
            num_cells = length(tumor_subcluster_cells_idx)

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
                        
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}



predict_CNV_via_HMM_on_whole_tumor_samples  <- function(infercnv_obj,
                                                        cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                        cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                        t=1e-6
                                                        ) {
    
    
    flog.info("predict_CNV_via_HMM_on_whole_tumor_samples")
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
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

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
                        
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_sample_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
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


get_predicted_CNV_regions <- function(infercnv_obj, by=c("consensus", "subcluster", "cell")) {
    by = match.arg(by)

    flog.info(sprintf("get_predicted_CNV_regions(%s)", by)) 
    
    cell_groups = NULL

    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("get_predicted_CNV_regions() - no subclusters defined, resetting reporting mode to consensus")
        by <- "consensus"
    }
    
    if (by == "consensus") {
        cell_groups = infercnv_obj@observation_grouped_cell_indices
    } else if (by == "subcluster") {
        cell_groups = unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=F)
    } else if (by == "cell") {
        cell_groups = lapply(unlist(infercnv_obj@observation_grouped_cell_indices), function(x) x) 
    }
    else {
        stop("Error, shouldn't get here ... bug")
    }
    
    cnv_regions = list()

    cnv_counter_start = 0
    for (cell_group_name in names(cell_groups)) {
        
        cell_group = cell_groups[[cell_group_name]]
        #flog.info(sprintf("cell group %s -> %s", cell_group_name, cell_group))

        flog.info(sprintf("-processing cell_group_name: %s, size: %d", cell_group_name, length(cell_group)))
                
        cell_group_mtx = infercnv_obj@expr.data[,cell_group,drop=F]
        cell_group_names = colnames(cell_group_mtx)

        state_consensus <- .get_state_consensus(cell_group_mtx)
        
        names(state_consensus) <- rownames(cell_group_mtx)
        cnv_gene_regions <- .define_cnv_gene_regions(state_consensus, infercnv_obj@gene_order, cnv_counter_start)
        cnv_ranges <- .get_cnv_gene_region_bounds(cnv_gene_regions)

        consensus_state_list = list(cell_group_name=cell_group_name,
                                    cells=cell_group_names,
                                    gene_regions=cnv_gene_regions,
                                    cnv_ranges=cnv_ranges)
        
        cnv_regions[[length(cnv_regions)+1]] = consensus_state_list

        cnv_counter_start = cnv_counter_start + length(cnv_gene_regions)
        
    }
    
    return(cnv_regions)
    
}


generate_cnv_region_reports <- function(infercnv_obj,
                                        output_filename_prefix,
                                        out_dir,
                                        by=c("consensus", "subcluster", "cell") ) {

    
    cnv_regions <- get_predicted_CNV_regions(infercnv_obj, by)

    ## cell clusters defined.
    cell_clusters_outfile = paste(out_dir, paste0(output_filename_prefix, ".cell_groupings"), sep="/")

    cell_clusters_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        cells = x$cells
        ret_df = data.frame(cell_group_name=cell_group_name, cell=cells)
        return(ret_df)
    })
    
    cell_clusters_df = do.call(rbind, cell_clusters_df)
    flog.info(sprintf("-writing cell clusters file: %s", cell_clusters_outfile))
    write.table(cell_clusters_df, file=cell_clusters_outfile, row.names=F, quote=F, sep="\t")
    
    ## regions DF:
    regions_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_regions.dat"), sep="/")
    
    regions_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        cnv_ranges = x$cnv_ranges
        ret_df = cbind(cell_group_name=cell_group_name, cnv_ranges)
        return(ret_df)
    })

    regions_df = do.call(rbind, regions_df)

    ## remove the neutral calls:
    regions_df = regions_df[regions_df$state != 3, ]
    flog.info(sprintf("-writing cnv regions file: %s", regions_outfile)) 
    write.table(regions_df, regions_outfile, row.names=F, sep="\t", quote=F)


    ## write the per-gene reports of cnv:
    if (by == "cell") {
        flog.warn("Note, HMM reporting is being done by 'cell', so this may use more memory, write more info to disk, take more time, ...")
    }
    gene_cnv_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        gene_region_list = x$gene_regions
        gene_region_names = names(gene_region_list)
        gene_region_df = lapply(gene_region_names, function(gene_region_name) {
            df = gene_region_list[[gene_region_name]]
            df = cbind(cell_group_name, gene_region_name, df)
            rownames(df) <- NULL
            return(df)
        })
        gene_region_df = do.call(rbind, gene_region_df)
        ## remove state 3 = neutral
        gene_region_df = gene_region_df[gene_region_df$state != 3,]
        return(gene_region_df)
    })
    gene_cnv_df = do.call(rbind, gene_cnv_df)
    
    ## write output file:
    gene_cnv_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_genes.dat"), sep="/") 
    flog.info(sprintf("-writing per-gene cnv report: %s", gene_cnv_outfile))
    write.table(gene_cnv_df, gene_cnv_outfile, row.names=F, sep="\t", quote=F)

    ## write file containing all genes that were leveraged in the predictions:
    gene_order_outfile = paste(out_dir, paste0(output_filename_prefix, ".genes_used.dat"), sep="/")
    flog.info(sprintf("-writing gene ordering info: %s", gene_order_outfile))
    write.table(infercnv_obj@gene_order, file=gene_order_outfile, quote=F, sep="\t")
    
    
    return
    
}
    
.get_state_consensus <- function(cell_group_matrix) {

    consensus  = apply(cell_group_matrix, 1, function(x) {
        t = table(x)
        names(t)[order(t, decreasing=T)[1]]
    })

    consensus <- as.numeric(consensus)
    
    return(consensus)
}


.define_cnv_gene_regions <- function(state_consensus, gene_order, cnv_region_counter) {

    regions = list()

    gene_names = rownames(gene_order)
    
    chrs = unique(gene_order$chr)
    for (chr in chrs) {
        gene_idx = which(gene_order$chr==chr)
        if (length(gene_idx) < 2) { break; }
        
        chr_states = state_consensus[gene_idx]
        prev_state = chr_states[1]
        ## pos_begin = paste(gene_order[gene_idx[1],,drop=T], collapse=",")
        pos_begin = gene_order[gene_idx[1],,drop=T]
        
        cnv_region_counter = cnv_region_counter + 1

        cnv_region_name = sprintf("%s-region_%d", chr, cnv_region_counter)
        current_cnv_region = data.frame(state=prev_state,
                                        gene=gene_names[gene_idx[1]],
                                        chr=pos_begin$chr,
                                        start=pos_begin$start,
                                        end=pos_begin$stop) 
        regions[[cnv_region_name]] = current_cnv_region
        
        for (i in seq(2,length(gene_idx))) {
            state = chr_states[i]
            pos_end = gene_order[gene_idx[i-1],,drop=T]
            next_gene_entry = data.frame(state=state,
                                         gene=gene_names[gene_idx[i]],
                                         chr=pos_end$chr,
                                         start=pos_end$start,
                                         end=pos_end$stop)
            
            if (state != prev_state) {
                ## state transition
                ## start new cnv region
                cnv_region_counter = cnv_region_counter + 1
                cnv_region_name = sprintf("%s-region_%d", chr, cnv_region_counter)
                regions[[cnv_region_name]] = next_gene_entry
            } else {
                ## append gene to current cnv region
                regions[[cnv_region_name]] = rbind(regions[[cnv_region_name]], next_gene_entry)
            }

            prev_state = state
        }
        
    }

    return(regions)
}

# cnv_gene_regions = .define_cnv_bounds(hmm_path, infercnv_obj@gene_order)

.get_cnv_gene_region_bounds <- function(cnv_gene_regions) {

    bounds = do.call(rbind, lapply(names(cnv_gene_regions), function(x) {
        cnv_name = x
        df = cnv_gene_regions[[cnv_name]]
        state = df$state[1]
        chr = df$chr[1]
        start = min(df$start)
        end = max(df$end)
        
        return(data.frame(cnv_name, state, chr, start, end))
    }) )

    rownames(bounds) <- NULL
    
    return(bounds)
}


    


## Adapted from the HiddenMarkov package:
Viterbi.dthmm.adj <- function (object, ...){
    x <- object$x

    if (length(x) < 2) {
        ## not enough run a trace on
        return(3); # neutral state
    }
    
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    n <- length(x)
    m <- nrow(object$Pi) # transition matrix
    nu <- matrix(NA, nrow = n, ncol = m)  # scoring matrix
    y <- rep(NA, n) # final trace
    pseudocount = 1e-20
    
    emissions <- matrix(NA, nrow = n, ncol = m) 
    
    ## ###############################
    ## restrict to constant variance to avoid nonsensical results:

    ## object$pm$sd = max(object$pm$sd)
    object$pm$sd = median(object$pm$sd) # max is too high
    
    ## ###############################    


    ## init first row

    emission <- pnorm(abs(x[1]-object$pm$mean)/object$pm$sd, log=T, lower.tail=F)
    emission <- 1 / (-1 * emission)
    emission <- emission / sum(emission)
    
    emissions[1,] <- log(emission)
    
    nu[1, ] <- log(object$delta) + # start probabilities
        emissions[1,]
    
    
    #nu[1, ] <- log(object$delta) + # start probabilities
    #    dfunc(x=x[1], # mean expr val for gene_1
    #          object$pm, # parameters (mean, sd) for norm dist
    #          HiddenMarkov:::getj(object$pn, 1), # NULL value
    #          log=TRUE) # returns p-values as log(p)
    
    logPi <- log(object$Pi) # convert transition matrix to log(p)
    
    for (i in 2:n) {
        
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        #nu[i, ] <- apply(matrixnu + logPi, 2, max) +
        #              dfunc(x=x[i], object$pm, getj(object$pn, i),
        #                    log=TRUE)

        
        emission <- pnorm(abs(x[i]-object$pm$mean)/object$pm$sd, log=T, lower.tail=F)
        emission <- 1 / (-1 * emission)
        emission <- emission / sum(emission)
        
        emissions[i, ] <- log(emission)
        
        nu[i, ] <- apply(matrixnu + logPi, 2, max) + emissions[i, ] 
                
    }
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")


    ## traceback
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1))
        y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])

    return(y)
}


assign_HMM_states_to_proxy_expr_vals <- function(infercnv_obj) {

    expr.data = infercnv_obj@expr.data

    expr.data[expr.data == 1] <- 0
    expr.data[expr.data == 2] <- 0.5
    expr.data[expr.data == 3] <- 1
    expr.data[expr.data == 4] <- 1.5
    expr.data[expr.data == 5] <- 2
    expr.data[expr.data == 6] <- 3

    infercnv_obj@expr.data <- expr.data

    return(infercnv_obj)

}

