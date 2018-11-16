#' @title apply_median_filtering
#' 
#' @description Apply a median filtering to the expression matrix within each tumor bounds
#'
#' @param infercnv_obj infercnv_object
#' 
#' @param window_size Size of the window side centered on the data point to filter.
#' 
#' @return infercnv_obj with median filtering applied to observations
#'
#' @export
#'

apply_median_filtering <- function(infercnv_obj,
                                   window_size=11) {

    if (window_size%%2 != 1 | window_size < 2) {
        flog.error("::apply_median_filtering: Error, window_size is an even or < 2. Please specify an odd number >= 3.")
    }
  
    half_window = (window_size - 1) / 2
    tumor_groupings = infercnv_obj@observation_grouped_cell_indices
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (tumor_type in names(tumor_groupings)) {
        tumor_indices = tumor_groupings[[ tumor_type ]] 
      
        for (chr in chrs) {
            chr_genes_indices = which(gene_chr_listing == chr)
          
            working_data = infercnv_obj@expr.data[chr_genes_indices, tumor_indices]
            xdim = dim(working_data)[1]
            ydim = dim(working_data)[2]
            results = working_data
            
            if (xdim >= window_size & ydim >= window_size) {
                # for (posx in ((half_window + 1):(xdim - (half_window + 1)))) {
                for (posx in 1:xdim) {
                    posxa <- ifelse(posx <= (half_window + 1), 1, (posx - (half_window + 1)))
                    posxb <- ifelse(posx >= (xdim - (half_window + 1)), xdim, (posx + (half_window + 1)))
                    for ( posy in 1:ydim) {
                        posya <- ifelse(posy <= (half_window + 1), 1, (posy - (half_window + 1)))
                        posyb <- ifelse(posy >= (ydim - (half_window + 1)), ydim, (posy + (half_window + 1)))
                        results[posx, posy] = median(working_data[posxa:posxb, posya:posyb])
                    }
                }
            }
            infercnv_obj@expr.data[chr_genes_indices, tumor_indices] = results
        }
    }
    
    return(infercnv_obj)
}


.apply_heatmap_median_filtering <- function(infercnv_obj,
                                            window_size=11) {

    if (window_size%%2 != 1 | window_size < 2) {
      flog.error("::apply_median_filtering: Error, window_size is an even or < 2. Please specify an odd number >= 3.")
    }
    
    half_window = (window_size - 1) / 2
    tumor_groupings = infercnv_obj@observation_grouped_cell_indices
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (tumor_type in names(tumor_groupings)) {
        
        indices = infercnv_obj@tumor_subclusters[["subclusters"]][[ tumor_type ]]
        if(is.list(indices)) {
            tumor_indices_list = indices
        }
        else { # is.vector(indices)
            tumor_indices_list = list(indices)
        }
        
        #tumor_indices = tumor_groupings[[ tumor_type ]] 
        for (tumor_indices in tumor_indices_list) {
        
            for (chr in chrs) {
                chr_genes_indices = which(gene_chr_listing == chr)
                
                working_data = infercnv_obj@expr.data[chr_genes_indices, tumor_indices]
                xdim = dim(working_data)[1]
                ydim = dim(working_data)[2]
                results = working_data
                
                if (xdim >= window_size & ydim >= window_size) {
                    # for (posx in ((half_window + 1):(xdim - (half_window + 1)))) {
                      for (posx in 1:xdim) {
                          posxa <- ifelse(posx <= (half_window + 1), 1, (posx - (half_window + 1)))
                          posxb <- ifelse(posx >= (xdim - (half_window + 1)), xdim, (posx + (half_window + 1)))
                          for ( posy in 1:ydim) {
                            posya <- ifelse(posy <= (half_window + 1), 1, (posy - (half_window + 1)))
                            posyb <- ifelse(posy >= (ydim - (half_window + 1)), ydim, (posy + (half_window + 1)))
                            results[posx, posy] = median(working_data[posxa:posxb, posya:posyb])
                        }
                    }
                }
              infercnv_obj@expr.data[chr_genes_indices, tumor_indices] = results
            }
        }
    }
    
    return(infercnv_obj)
}
