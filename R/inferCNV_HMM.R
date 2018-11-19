
get_spike_dists <- function(hspike_obj) {

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

