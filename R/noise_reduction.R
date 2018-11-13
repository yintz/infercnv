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
    
    for (tumor_type in names(tumor_groupings)) {
        tumor_indices = tumor_groupings[[ tumor_type ]] 
        
        working_data = infercnv_obj@expr.data[, tumor_indices]
        xdim = dim(working_data)[1]
        ydim = dim(working_data)[2]
        results = working_data
        if (xdim >= window_size & ydim >= window_size) {
            for (posx in ((half_window + 1):(xdim - (half_window + 1)))) {
                for ( posy in ((half_window + 1):(ydim - (half_window + 1)))) {
                    results[posx, posy] = median(working_data[(posx - half_window):(posx + half_window), (posy - half_window):(posy + half_window)])
                }
            }
        }
        infercnv_obj@expr.data[, tumor_indices] = results
    }
    
    return(infercnv_obj)
}
