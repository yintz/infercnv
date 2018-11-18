
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
        cnv = info$cnv
        chr_gene_idx = which(hspike_obj@gene_order$chr == chr_name)
        gene.expr = spike.expr.data[chr_gene_idx, ]

        gene_expr_by_cnv[[cnv]] = c(gene_expr_by_cnv[[cnv]], gene.expr)
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
        
