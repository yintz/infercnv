
define_signif_tumor_subclusters <- function(infercnv_obj, p_val, hclust_method, use_random_trees=FALSE) {
    
    flog.info(sprintf("define_signif_tumor_subclusters(p_val=%g", p_val))
    
    tumor_groups <- infercnv_obj@observation_grouped_cell_indices

    res = list()
        
    for (tumor_group in names(tumor_groups)) {

        flog.info(sprintf("define_signif_tumor_subclusters(), tumor: %s", tumor_group))
        
        tumor_group_idx <- tumor_groups[[ tumor_group ]]
        tumor_expr_data <- infercnv_obj@expr.data[,tumor_group_idx]

                
        if (use_random_trees) {
            tumor_subcluster_info <- .single_tumor_subclustering_random_trees(tumor_group, tumor_group_idx, tumor_expr_data, p_val, hclust_method)
        } else {
            tumor_subcluster_info <- .single_tumor_subclustering(tumor_group, tumor_group_idx, tumor_expr_data, p_val, hclust_method)
        }
        
        res$hc[[tumor_group]] <- tumor_subcluster_info$hc
        res$subclusters[[tumor_group]] <- tumor_subcluster_info$subclusters

    }
    
     
    infercnv_obj@tumor_subclusters <- res
    
    return(infercnv_obj)
}


.single_tumor_subclustering <- function(tumor_name, tumor_group_idx, tumor_expr_data, p_val, hclust_method) {
    
    tumor_subcluster_info = list()
    
    hc <- hclust(dist(t(tumor_expr_data)), method=hclust_method)
    
    tumor_subcluster_info$hc = hc
    
    heights = hc$height

    mu = mean(heights)
    sigma = sd(heights)

    cut_height = qnorm(p=1-p_val, mean=mu, sd=sigma)
    
    flog.info(sprintf("cut height based on p_val(%g) = %g", p_val, cut_height))
    
    grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height
    
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








#### Experimental stuff below that might get deleted altogether


.single_tumor_subclustering_random_trees <- function(tumor_name, tumor_group_idx, tumor_expr_data, p_val, hclust_method) {

    tumor_subcluster_info = list()
    rand_params_info = .parameterize_random_cluster_heights(tumor_expr_data, hclust_method)

    tumor_subcluster_info$rand_params_info <- rand_params_info
    
    hc = rand_params_info$h_obs
    tumor_subcluster_info$hc = hc
        
    cut_height = .get_tree_height_via_ecdf(p_val, rand_params_info)
    
    flog.info(sprintf("cut height based on p_val(%g) = %g", p_val, cut_height))
            
    grps <- cutree(hc, h=cut_height) # will just be one cluster if height > max_height
    
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


.parameterize_random_cluster_heights <- function(expr_matrix, hclust_method, plot=F) {
    
    ## inspired by: https://www.frontiersin.org/articles/10.3389/fgene.2016.00144/full
    
    h_obs <- hclust(dist(t(expr_matrix)), method=hclust_method)
    
    num_samples = ncol(expr_matrix)

    norm_heights = h_obs$height
    
    ## run 10 randomizations
    nrand_iter=10
    rand_hc_example = NULL
    
    df = do.call(cbind, lapply(1:nrand_iter, function(x) {
        flog.info(sprintf("randomization(%d)", x))

        ## shuffle genes for each cell 
        
        expr_rand = apply(expr_matrix, 1, function(row) {
            row <- sample(row,size=length(row), replace=F)
            return(row)
        })
        flog.info("clustering randomized matrix")
        h_rand = hclust(dist(expr_rand), method='complete')
        rand_hc_example <<- h_rand
                
        data.frame(h_rand$height)
    }))
    
    rand_height_dist = as.numeric(as.matrix(df))

    
    e = ecdf(rand_height_dist)
    #pvals = sapply(h_obs_heights, function(x) 1-e(x))
    
    mu = mean(rand_height_dist)
    sigma = sd(rand_height_dist)
        
    params_list <- list(h_obs=h_obs,
                        norm_heights=norm_heights,
                        rand_mu=mu,
                        rand_sigma=sigma,
                        rand_height_dist=rand_height_dist,
                        ecdf=e,
                        h_rand_ex = rand_hc_example,
                        num_samples=num_samples
                        )
    
    if (plot) {
        .plot_tree_height_dist(params_list)
    }
    
    
    return(params_list)
    
}

.plot_tree_height_dist <- function(params_list, plot_title='tree_heights') {

    mf = par(mfrow=(c(3,2)))

    ## density plot
    rand_height_density = density(params_list$rand_height_dist)
    h_obs_height_density = density(params_list$norm_heights)
    
    xlim=range(h_obs_height_density$x, rand_height_density$x)
    ylim=range(h_obs_height_density$y, rand_height_density$y)
    plot(rand_height_density, xlim=xlim, ylim=ylim, main=paste(plot_title, "density"))
    points(h_obs_height_density)

    ## p-val for observed heights based on ecdf 
    h = params_list$h_obs$height
    pvals = 1 - sapply(h, params_list$ecdf)
    plot(h, pvals, main=paste(plot_title, "height pvals"), xlab='height', ylab='pval')
    
    ## plot the clustering
    h_obs = params_list$h_obs
    h_obs$labels <- NULL #because they're too long to display
    plot(h_obs)
    
    ## plot a random example:
    plot(params_list$h_rand_ex)

    r = range(h, params_list$rand_height_dist)
    qqplot(h, params_list$rand_height_dist, xlim=r, ylim=r)
    abline(a=0,b=1, col='red')
    
    plot.new() # last one.
    
    par(mf)
        
}

.get_tree_height_via_ecdf <- function(p_val, params_list) {
    
    h = quantile(params_list$ecdf, probs=1-p_val)

    return(h)
}


        

    
