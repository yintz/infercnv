
define_signif_tumor_subclusters <- function(infercnv_obj, p_val, hclust_method, cluster_by_groups, partition_method, restrict_to_DE_genes=FALSE) {
    
    flog.info(sprintf("define_signif_tumor_subclusters(p_val=%g", p_val))
    
    # tumor_groups <- infercnv_obj@observation_grouped_cell_indices

    res = list()

    if (restrict_to_DE_genes) {
        normal_expr_data = infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices) ]
    }
    
    tumor_groups = list()

    if (cluster_by_groups) {
        tumor_groups <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    }
    else {
        if(length(infercnv_obj@reference_grouped_cell_indices) > 0) {
            tumor_groups <- list(all_observations=unlist(infercnv_obj@observation_grouped_cell_indices, use.names=FALSE), all_references=unlist(infercnv_obj@reference_grouped_cell_indices, use.names=FALSE))
        }
        else {
            tumor_groups <- list(all_observations=unlist(infercnv_obj@observation_grouped_cell_indices, use.names=FALSE))
        }
    }

    for (tumor_group in names(tumor_groups)) {

        flog.info(sprintf("define_signif_tumor_subclusters(), tumor: %s", tumor_group))
        
        tumor_group_idx <- tumor_groups[[ tumor_group ]]
        tumor_group_cell_names <- colnames(infercnv_obj@expr.data[,tumor_group_idx])
        tumor_expr_data <- infercnv_obj@expr.data[,tumor_group_idx, drop=FALSE]

        if (restrict_to_DE_genes) {
            p_vals <- .find_DE_stat_significance(normal_expr_data, tumor_expr_data)
            
            DE_gene_idx = which(p_vals < p_val)
            tumor_expr_data = tumor_expr_data[DE_gene_idx, , drop=FALSE]
            
        }
        
        tumor_subcluster_info <- .single_tumor_subclustering(tumor_group, tumor_group_idx, tumor_group_cell_names, tumor_expr_data, p_val, hclust_method, partition_method)
        
        res$hc[[tumor_group]] <- tumor_subcluster_info$hc
        res$subclusters[[tumor_group]] <- tumor_subcluster_info$subclusters

    }
         
    infercnv_obj@tumor_subclusters <- res


    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- define_signif_tumor_subclusters(infercnv_obj@.hspike, p_val, hclust_method, cluster_by_groups, partition_method, restrict_to_DE_genes)
    }
        
    
    return(infercnv_obj)
}



.single_tumor_subclustering <- function(tumor_name, tumor_group_idx, tumor_group_cell_names, tumor_expr_data, p_val, hclust_method,
                                        partition_method=c('qnorm', 'pheight', 'qgamma', 'shc', 'none')
                                        ) {
    
    partition_method = match.arg(partition_method)
    
    tumor_subcluster_info = list()
    
    if (ncol(tumor_expr_data) > 2) {

        hc <- hclust(dist(t(tumor_expr_data)), method=hclust_method)
        
        tumor_subcluster_info$hc = hc
        
        heights = hc$height

        grps <- NULL
        
        if (partition_method == 'pheight') {

            cut_height = p_val * max(heights)
            flog.info(sprintf("cut height based on p_val(%g) = %g and partition_method: %s", p_val, cut_height, partition_method))
            grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height

            
        } else if (partition_method == 'qnorm') {

            mu = mean(heights)
            sigma = sd(heights)
            
            cut_height = qnorm(p=1-p_val, mean=mu, sd=sigma)
            flog.info(sprintf("cut height based on p_val(%g) = %g and partition_method: %s", p_val, cut_height, partition_method))
            grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height
            
        } else if (partition_method == 'qgamma') {

            # library(fitdistrplus)
            gamma_fit = fitdist(heights, 'gamma')
            shape = gamma_fit$estimate[1]
            rate = gamma_fit$estimate[2]
            cut_height=qgamma(p=1-p_val, shape=shape, rate=rate)
            flog.info(sprintf("cut height based on p_val(%g) = %g and partition_method: %s", p_val, cut_height, partition_method))
            grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height
            
        #} else if (partition_method == 'shc') {
        #    
        #    grps <- .get_shc_clusters(tumor_expr_data, hclust_method, p_val)
            
        } else if (partition_method == 'none') {
            
            grps <- cutree(hc, k=1)
            
        } else {
            stop("Error, not recognizing parition_method")
        }
        
        # cluster_ids = unique(grps)
        # flog.info(sprintf("cut tree into: %g groups", length(cluster_ids)))
        
        tumor_subcluster_info$subclusters = list()
        
        ordered_idx = tumor_group_idx[hc$order]
        s = split(grps,grps)
        flog.info(sprintf("cut tree into: %g groups", length(s)))

        start_idx = 1

        # for (g in cluster_ids) {
        for (g in names(s)) {
            
            split_subcluster = paste0(tumor_name, "_s", g)
            flog.info(sprintf("-processing %s,%s", tumor_name, split_subcluster))
            
            # subcluster_indices = tumor_group_idx[which(grps == g)]
            end_idx = start_idx + length(s[[g]]) - 1
            subcluster_indices = tumor_group_idx[hc$order[start_idx:end_idx]]
            subcluster_names = tumor_group_cell_names[hc$order[start_idx:end_idx]]
            start_idx = end_idx + 1
            
            tumor_subcluster_info$subclusters[[ split_subcluster ]] = subcluster_indices
            names(tumor_subcluster_info$subclusters[[ split_subcluster ]]) = subcluster_names

        }
    }
    else {
        tumor_subcluster_info$hc = NULL # can't make hc with a single element, even manually, need to have workaround in plotting step
        tumor_subcluster_info$subclusters[[paste0(tumor_name, "_s1") ]] = tumor_group_idx
    }
    
    return(tumor_subcluster_info)
}


#.get_shc_clusters <- function(tumor_expr_data, hclust_method, p_val) { 
#
# library(sigclust2)
#    
#    flog.info(sprintf("defining groups using shc, hclust_method: %s, p_val: %g", hclust_method, p_val))
#    
#    shc_result = sigclust2::shc(t(tumor_expr_data), metric='euclidean', linkage=hclust_method, alpha=p_val)
#
#    cluster_idx = which(shc_result$p_norm <= p_val)
#        
#    grps = rep(1, ncol(tumor_expr_data))
#    names(grps) <- colnames(tumor_expr_data)
#    
#    counter = 1
#    for (cluster_id in cluster_idx) {
#        labelsA = unlist(shc_result$idx_hc[cluster_id,1])
#
#        labelsB = unlist(shc_result$idx_hc[cluster_id,2])
#
#        counter = counter + 1
#        grps[labelsB] <- counter
#    }
#    
#    return(grps)
#}
        



.find_DE_stat_significance <- function(normal_matrix, tumor_matrix) {
    
    run_t_test<- function(idx) {
        vals1 = unlist(normal_matrix[idx,,drop=TRUE])
        vals2 = unlist(tumor_matrix[idx,,drop=TRUE])
        
        ## useful way of handling tests that may fail:
        ## https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html

        res = try(t.test(vals1, vals2), silent=TRUE)
        
        if (is(res, "try-error")) return(NA) else return(res$p.value)
        
    }

    pvals = sapply(seq(nrow(normal_matrix)), run_t_test)

    return(pvals)
}




##### Below is deprecated.... use inferCNV_tumor_subclusters.random_smoothed_trees
## Random Trees

.partition_by_random_trees <- function(tumor_name, tumor_expr_data, hclust_method, p_val) {

    grps <- rep(sprintf("%s.%d", tumor_name, 1), ncol(tumor_expr_data))
    names(grps) <- colnames(tumor_expr_data)

    grps <- .single_tumor_subclustering_recursive_random_trees(tumor_expr_data, hclust_method, p_val, grps)

    
    return(grps)

}


.single_tumor_subclustering_recursive_random_trees <- function(tumor_expr_data, hclust_method, p_val, grps.adj, min_cluster_size_recurse=10) {

    tumor_clade_name = unique(grps.adj[names(grps.adj) %in% colnames(tumor_expr_data)])
    message("unique tumor clade name: ", tumor_clade_name)
    if (length(tumor_clade_name) > 1) {
        stop("Error, found too many names in current clade")
    }
    
    hc <- hclust(dist(t(tumor_expr_data)), method=hclust_method)

    rand_params_info = .parameterize_random_cluster_heights(tumor_expr_data, hclust_method)

    h_obs = rand_params_info$h_obs
    h = h_obs$height
    max_height = rand_params_info$max_h
    
    max_height_pval = 1
    if (max_height > 0) {
        ## important... as some clades can be fully collapsed (all identical entries) with zero heights for all
        e = rand_params_info$ecdf
        max_height_pval = 1- e(max_height)
    }

    #message(sprintf("Lengths(h): %s", paste(h, sep=",", collapse=",")))
    #message(sprintf("max_height_pval: %g", max_height_pval))
    
    if (max_height_pval <= p_val) {
        ## keep on cutting.
        cut_height = mean(c(h[length(h)], h[length(h)-1]))
        message(sprintf("cutting at height: %g",  cut_height))
        grps = cutree(h_obs, h=cut_height)
        print(grps)
        uniqgrps = unique(grps)
        
        message("unique grps: ", paste0(uniqgrps, sep=",", collapse=","))
        for (grp in uniqgrps) {
            grp_idx = which(grps==grp)
            
            message(sprintf("grp: %s  contains idx: %s", grp, paste(grp_idx,sep=",", collapse=","))) 
            df = tumor_expr_data[,grp_idx,drop=FALSE]
            ## define subset.
            subset_cell_names = colnames(df)
            
            subset_clade_name = sprintf("%s.%d", tumor_clade_name, grp)
            grps.adj[names(grps.adj) %in% subset_cell_names] <- subset_clade_name

            if (length(grp_idx) > min_cluster_size_recurse) {
                ## recurse
                grps.adj <- .single_tumor_subclustering_recursive_random_trees(tumor_expr_data=df,
                                                                               hclust_method=hclust_method,
                                                                               p_val=p_val,
                                                                               grps.adj)
            } else {
                message("paritioned cluster size too small to recurse further")
            }
        }
    } else {
        message("No cluster pruning: ", tumor_clade_name)
    }
    
    return(grps.adj)
}


.parameterize_random_cluster_heights <- function(expr_matrix, hclust_method, plot=TRUE) {
    
    ## inspired by: https://www.frontiersin.org/articles/10.3389/fgene.2016.00144/full

    t_tumor.expr.data = t(expr_matrix) # cells as rows, genes as cols
    d = dist(t_tumor.expr.data)

    h_obs = hclust(d, method=hclust_method)

        
    # permute by chromosomes
    permute_col_vals <- function(df) {

        num_cells = nrow(df)

        for (i in seq(ncol(df) ) ) {
            
            df[, i] = df[sample(x=seq_len(num_cells), size=num_cells, replace=FALSE), i]
        }
        
        df
    }
    
    h_rand_ex = NULL
    max_rand_heights = c()
    num_rand_iters=100
    for (i in seq_len(num_rand_iters)) {
        #message(sprintf("iter i:%d", i))
        rand.tumor.expr.data = permute_col_vals(t_tumor.expr.data)
        
        rand.dist = dist(rand.tumor.expr.data)
        h_rand <- hclust(rand.dist, method=hclust_method)
        h_rand_ex = h_rand
        max_rand_heights = c(max_rand_heights, max(h_rand$height))
    }
    
    h = h_obs$height

    max_height = max(h)
    
    message(sprintf("Lengths for original tree branches (h): %s", paste(h, sep=",", collapse=",")))
    message(sprintf("Max height: %g", max_height))

    message(sprintf("Lengths for max heights: %s", paste(max_rand_heights, sep=",", collapse=",")))
    
    e = ecdf(max_rand_heights)
    
    pval = 1- e(max_height)
    message(sprintf("pval: %g", pval))
    
    params_list <- list(h_obs=h_obs,
                        max_h=max_height,
                        rand_max_height_dist=max_rand_heights,
                        ecdf=e,
                        h_rand_ex = h_rand_ex
                        )
    
    if (plot) {
        .plot_tree_height_dist(params_list)
    }
    
    
    return(params_list)
    
}


.plot_tree_height_dist <- function(params_list, plot_title='tree_heights') {

    mf = par(mfrow=(c(3,1)))

    ## density plot
    rand_height_density = density(params_list$rand_max_height_dist)
    
    xlim=range(params_list$max_h, rand_height_density$x)
    ylim=range(rand_height_density$y)
    plot(rand_height_density, xlim=xlim, ylim=ylim, main=paste(plot_title, "density"))
    abline(v=params_list$max_h, col='red')

        
    ## plot the clustering
    h_obs = params_list$h_obs
    h_obs$labels <- NULL #because they're too long to display
    plot(h_obs)
    
    ## plot a random example:
    h_rand_ex = params_list$h_rand_ex
    h_rand_ex$labels <- NULL
    plot(h_rand_ex)
            
    par(mf)
        
}

.get_tree_height_via_ecdf <- function(p_val, params_list) {
    
    h = quantile(params_list$ecdf, probs=1-p_val)

    return(h)
}


