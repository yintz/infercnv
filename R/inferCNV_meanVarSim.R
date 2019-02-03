.get_simulated_cell_matrix_using_meanvar_trend <- function(infercnv_obj, gene_means, num_cells, include.dropout=FALSE) {

    # should be working on the total sum count normalized data.
    # model the mean variance relationship


    mean_var_table = .get_mean_var_table(infercnv_obj)

    logm = log(mean_var_table$m + 1)
    logv = log(mean_var_table$v + 1)

    mean_var_spline = smooth.spline(logv ~ logm)

    ngenes = length(gene_means)


    dropout_logistic_params <- NULL

    if (include.dropout) {

        mean_p0_table <- .get_mean_vs_p0_table(infercnv_obj)

        dropout_logistic_params <- infercnv:::.get_logistic_params(mean_p0_table)
    }

    spike_cell_names = paste0('sim_cell_', 1:num_cells)

    sim_cell_matrix = matrix(rep(0,ngenes*num_cells), nrow=ngenes)
    rownames(sim_cell_matrix) = names(gene_means)
    colnames(sim_cell_matrix) = spike_cell_names

    sim_expr_vals <- function(gene_idx) {
        m = gene_means[gene_idx]
        return(.sim_expr_val_mean_var(m, mean_var_spline, dropout_logistic_params))
    }

    for (i in 1:num_cells) {
        newvals = sapply(1:ngenes, FUN=sim_expr_vals)
        sim_cell_matrix[,i] = newvals
    }

    return(sim_cell_matrix)
}


##' @keywords internal
##' @noRd
##'

.sim_expr_val_mean_var <- function(m,  mean_var_spline, dropout_logistic_params) {

    # include drop-out prediction

    val = 0
    if (m > 0) {
        logm = log(m+1)
        pred_log_var = predict(mean_var_spline, logm)$y

        var = max(exp(pred_log_var)-1, 0)

        val = round(max(rnorm(n=1, mean=m, sd=sqrt(var)), 0))

        if (! is.null(dropout_logistic_params)) {

            dropout_prob <- predict(dropout_logistic_params$spline, log(m))$y[1]

            if (runif(1) <= dropout_prob) {
                ## a drop-out
                val = 0
            }
        }
    }

    return(val)
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

##' .get_spike_in_average_bounds()
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

