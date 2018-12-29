
define_signif_tumor_subclusters <- function(infercnv_obj, p_val, hclust_method, partition_method) {
    
    flog.info(sprintf("define_signif_tumor_subclusters(p_val=%g", p_val))
    
    tumor_groups <- infercnv_obj@observation_grouped_cell_indices

    res = list()
        
    for (tumor_group in names(tumor_groups)) {

        flog.info(sprintf("define_signif_tumor_subclusters(), tumor: %s", tumor_group))
        
        tumor_group_idx <- tumor_groups[[ tumor_group ]]
        tumor_expr_data <- infercnv_obj@expr.data[,tumor_group_idx]

        tumor_subcluster_info <- .single_tumor_subclustering(tumor_group, tumor_group_idx, tumor_expr_data, p_val, hclust_method, partition_method)
                
        res$hc[[tumor_group]] <- tumor_subcluster_info$hc
        res$subclusters[[tumor_group]] <- tumor_subcluster_info$subclusters

    }
         
    infercnv_obj@tumor_subclusters <- res
    
    return(infercnv_obj)
}


.single_tumor_subclustering <- function(tumor_name, tumor_group_idx, tumor_expr_data, p_val, hclust_method,
                                        partition_method=c('shc', 'qnorm', 'pheight', 'qgamma') ) {
    
    partition_method = match.arg(partition_method)
    
    tumor_subcluster_info = list()
    
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

        library(fitdistrplus)
        gamma_fit = fitdist(heights, 'gamma')
        shape = gamma_fit$estimate[1]
        rate = gamma_fit$estimate[2]
        cut_height=qgamma(p=0.9, shape=shape, rate=rate)
        flog.info(sprintf("cut height based on p_val(%g) = %g and partition_method: %s", p_val, cut_height, partition_method))
        grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height
        
    } else if (partition_method == 'shc') {
        
        grps <- .get_shc_clusters(tumor_expr_data, hclust_method, p_val)
    }
    else {
        stop("Error, not recognizing parition_method")
    }
        
    cluster_ids = unique(grps)
    flog.info(sprintf("cut tree into: %g groups", length(cluster_ids)))
    
    tumor_subcluster_info$subclusters = list()
    
    for (g in cluster_ids) {
        split_subcluster = paste0(tumor_name, "_s", g)
        flog.info(sprintf("-processing %s,%s", tumor_name, split_subcluster))
        
        subcluster_indices = tumor_group_idx[which(grps == g)]
        
        tumor_subcluster_info$subclusters[[ split_subcluster ]] = subcluster_indices
        
    }
    
    return(tumor_subcluster_info)
}


.get_shc_clusters <- function(tumor_expr_data, hclust_method, p_val) { 

    library(sigclust2)
    
    flog.info(sprintf("defining groups using shc, hclust_method: %s, p_val: %g", hclust_method, p_val))
    
    shc_result = shc(t(tumor_expr_data), metric='euclidean', linkage=hclust_method, alpha=p_val)

    cluster_idx = which(shc_result$p_norm <= p_val)
        
    grps = rep(1, ncol(tumor_expr_data))
    names(grps) <- colnames(tumor_expr_data)
    
    counter = 1
    for (cluster_id in cluster_idx) {
        labelsA = unlist(shc_result$idx_hc[cluster_id,1])

        labelsB = unlist(shc_result$idx_hc[cluster_id,2])

        counter = counter + 1
        grps[labelsB] <- counter
    }
    
    return(grps)
}
        
