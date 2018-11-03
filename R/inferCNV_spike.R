
#' @title spike_in_variation_chrs()
#'
#' Adds a 'SPIKE'-in to the observations set at different thresholds of loss/gain to
#' aid in tracking the effect of infercnv operations and for defining the final scaling.
#'
#' 
#' @param infercnv_obj An infercnv object populated with raw count data
#'
#' @param spike_in_chrs : define the chromsomes that will serve as signal for gain/loss
#'                        default: picks chrosomes in order of size
#'
#' @param min_genes_per_chr  : default 100
#' 
#' @param spike_in_multiplier_vec : factors that define relative expression for gain/loss
#'                                  and must match ordering of spike_in_chrs above
#'                                  default: c(0.01, 2.0)
#' 
#' @param max_cells max number of cells to incorporate in the spike-in
#'
#' @export

spike_in_variation_chrs <- function(infercnv_obj,
                                    spike_in_chrs=NULL,
                                    spike_in_multiplier_vec=c(0.01, 2.0),
                                    max_cells=100,
                                    min_genes_per_chr=100) {


    if (is.null(spike_in_chrs)) {
        num_chrs_want = length(spike_in_multiplier_vec)
        spike_in_chrs = .select_longest_chrs(infercnv_obj, num_chrs_want) 
        flog.info(paste("Selecting longest chrs for adding spike:", paste(spike_in_chrs, collapse=",")))
    } else {
        flog.info(paste("Using specified chrs for adding spike:", paste(spike_in_chrs, collapse=",")))
    }
    
    min_genes_selected = min_genes_per_chr
    
    ## get the gene ordering:
    gene_selection_listing = list()
    for (chr in spike_in_chrs) {
        chr_genes = which(infercnv_obj@gene_order$chr == chr)
        if (length(chr_genes) < min_genes_selected) {
            flog.error(sprintf("Error, have %d genes found for chr %s, < min %d required",
                               length(chr_genes),
                               chr,
                               min_genes_selected))
            stop("error")
        }
        gene_selection_listing[[chr]] = chr_genes
    }

    infercnv_obj <- .spike_in_variation_genes_via_modeling(infercnv_obj, gene_selection_listing, spike_in_multiplier_vec, max_cells=max_cells)
    
    return(infercnv_obj)
}


##' @title .spike_in_variation_genes_via_modeling()
##'
##' Creates the spike-in based on a list of genes. 
##'
##' @param infercnv_obj An infercnv object populated with raw count data
##'
##' @param gene_selection_listing : list of [[chr]] = [1,2,3,4,...] corresponding to indices (rows) of genes,
##'                                 and should match order of spike_in_multiplier_vec below.
##'
##' @param spike_in_multiplier_vec : vector of factors corresponding to gain/loss multipliers matching order in the gene list above.
##'
##' @param max_cells : max number of cells to include in the spike-in
##' 
##' @keywords internal
##' @noRd
##'
.spike_in_variation_genes_via_modeling <- function(infercnv_obj, gene_selection_listing, spike_in_multiplier_vec, max_cells=max_cells) {

    mvtable = .get_mean_var_table(infercnv_obj)
    normal_cells_idx = infercnv::get_reference_grouped_cell_indices(infercnv_obj)
    normal_cells_expr = infercnv_obj@expr.data[,normal_cells_idx]

    ## apply spike-in multiplier vec
    for (i in 1:length(spike_in_multiplier_vec)) {
        
        gene_indices = gene_selection_listing[[i]]
        multiplier = spike_in_multiplier_vec[i]
        
        normal_cells_expr[gene_indices, ] = normal_cells_expr[gene_indices, ] * multiplier
        
    }

    ## get simulated matrix
    sim_matrix = .get_simulated_cell_matrix(mvtable, normal_cells_expr, max_cells)

    ## integrate into expr data and count data matrices
    ncol_begin = ncol(infercnv_obj@expr.data) + 1
    ncol_end = ncol_begin + max_cells - 1
    
    infercnv_obj@expr.data = cbind( infercnv_obj@expr.data, sim_matrix )
    infercnv_obj@count.data = cbind( infercnv_obj@count.data, sim_matrix ) # just so it validates... not useful, otherwise
    
    infercnv_obj@observation_grouped_cell_indices[['SPIKE']] = ncol_begin:ncol_end

    validate_infercnv_obj(infercnv_obj)

    
    return(infercnv_obj)
}


##' .get_simulated_cell_matrix()
##'
##' generates a simulated grouping of cells vs. genes based on a cell expression matrix and
##' the mean/variance relationship for all genes in all cell groupings.
##'
##' Cells are simulated as so:
##'    A random cell is selected from the normal cell expression matrix.
##'    The expression of each gene is treated as a targeted mean expression value.
##'    The variance is chosen based on a spline fit to the mean/variance relationship provided.
##'    A random expression value is generated from a normal distribution with corresponding (mean, variance)
##'
##' Genes are named according to the input expression matrix, and cells are named 'spike_{number}'.
##' 
##' @param mean_var_table : a data.frame containing three columns: group_name, mean, variance of expression per gene per grouping.
##'
##' @param normal_cell_expr : expression matrix of normal cells to guide the simulation. Should be total sum normalized.
##'
##' @param num_cells : number of cells to simulate
##'
##' @return matrix containing simulated expression values.
##'
##' @keywords internal
##' @noRd
##' 

.get_simulated_cell_matrix <- function(mean_var_table, normal_cell_expr, num_cells) {
    
    # should be working on the total sum count normalized data.
    # model the mean variance relationship

    s = smooth.spline(log2(mean_var_table$m+1), log2(mean_var_table$v+1))
        
    spike_cell_names = paste0("spike_", 1:num_cells)
    
    ngenes = nrow(normal_cell_expr)
    
    sim_expr_val <- function(gene_idx, rand_cell_idx) {
        m = normal_cell_expr[gene_idx, rand_cell_idx]
        #v = predict(s, log2(m+1))$y
        #v = max(0, 2^v-1)
        #val = max(0, rnorm(n=1, mean=m, sd=sqrt(v)))
        val = max(0, rnbinom(n=1, mu=m, size=0.1)) 
        return(val)
    }
    
    sim_cell_matrix = matrix(rep(0,ngenes*num_cells), nrow=ngenes)
    rownames(sim_cell_matrix) = rownames(normal_cell_expr)
    colnames(sim_cell_matrix) = spike_cell_names
    
    for (i in 1:num_cells) {
        rand_cell_idx = floor(runif(1) * ncol(normal_cell_expr)+1)
        newvals = sapply(1:ngenes, FUN=sim_expr_val, rand_cell_idx)
        sim_cell_matrix[,i] = newvals
    }

    return(sim_cell_matrix)
}
    
##' .get_mean_var_table()
##'
##' Computes the gene mean/variance table based on all defined cell groupings (reference and observations)
##'
##' @param infercnv_obj An infercnv object populated with raw count data
##'
##' @return data.frame with 3 columns: group_name, mean, variance
##' 
##'
##' @keywords internal
##' @noRd
##'

.get_mean_var_table <- function(infercnv_obj) {

    group_indices = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)

    mean_var_table = NULL
    
    for (group_name in names(group_indices)) {
        flog.info(sprintf("processing group: %s", group_name))
        expr.data = infercnv_obj@expr.data[, group_indices[[ group_name ]] ]
        m = rowMeans(expr.data)
        v = apply(expr.data, 1, var)
        if (is.null(mean_var_table)) {
            mean_var_table = data.frame(g=group_name, m=m, v=v)
        } else {
            mean_var_table = rbind(mean_var_table, data.frame(g=group_name, m=m, v=v))
        }
    }
        
    return(mean_var_table)
}

##' get_spike_in_average_bounds()
##'
##' return mean bounds for expression of all cells in the spike-in
##'
##' @param infercnv_obj An infercnv object populated with raw count data
##'
##' @return c(left_bound, right_bound)
##' 
##' @keywords internal
##' @noRd
##'


.get_spike_in_average_bounds <- function(infercnv_obj) {

    spike_in_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ 'SPIKE' ]]
    spike.expr.data = infercnv_obj@expr.data[,spike_in_cell_idx]

    bounds = .get_average_bounds(spike.expr.data)

    return(bounds)
}


#' remove_spike()
#'
#' Removes the spiked-in group named 'SPIKE' from the infercnv_obj
#'
#' @param infercnv_obj An infercnv object populated with raw count data
#'
#' @return infercnv_obj 
#'
#' @export

remove_spike <- function(infercnv_obj) {

    flog.info("Removing spike")
    
    spike_in_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ 'SPIKE' ]]

    infercnv_obj@expr.data = infercnv_obj@expr.data[, -spike_in_cell_idx]

    infercnv_obj@observation_grouped_cell_indices[[ 'SPIKE' ]] <- NULL # deletes it.

    return(infercnv_obj)

}



#' scale_cnv_by_spike()
#'
#' Scales expression data according to the expression value bounds in the SPIKE group.
#'
#' Assumes data is centered at 1
#' Expression below 1 is scaled according to the left spike bound set to zero.
#' Expression above 1 is scaled according to the right spike bound set to two.
#'
#' @param infercnv_obj An infercnv object populated with raw count data
#'
#' @return infercnv_obj
#'
#' @export


scale_cnv_by_spike <- function(infercnv_obj) {

    # everything here should be centered at 1 (no change).
    
    spike_bounds = .get_spike_in_average_bounds(infercnv_obj)

    left_bound = spike_bounds[1]
    right_bound = spike_bounds[2]

                                        # zero gets set to left bound
                                        # right bound gets set to 2x

    scale_by_spike <- function(x) {
        if (x < 1) {
            x = 1 - ( (1-x)/(1-left_bound) )
            if (x < 0) { x = 0 }
        } else if (x > 1) {
            x = 1 + ( (x-1) / (right_bound - 1) )
        }
        return(x)
    }

    infercnv_obj@expr.data <- apply(infercnv_obj@expr.data, 1:2, scale_by_spike)

    return(infercnv_obj)
}


# selects the specified number of chrs having the largest number of (expressed) genes
.select_longest_chrs <- function(infercnv_obj, num_chrs_want) {

    # get count of chrs
    counts = infercnv_obj@gene_order %>% count(.data$chr, sort=TRUE)

    return(counts$chr[1:num_chrs_want])
        
}
