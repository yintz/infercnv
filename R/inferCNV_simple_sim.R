
##' .get_simulated_cell_matrix()
##'
##' generates a simulated grouping of cells vs. genes based on a cell expression matrix and
##' the mean/variance relationship for all genes in all cell groupings.
##'
##' Cells are simulated as so:
##'    The mean for genes in the normal cells are computed
##'    A random expression value is chosen for each gene using a negative binomial distribution with specified common dispersion
##'
##' Genes are named according to the input expression matrix, and cells are named 'spike_{number}'.
##'
##' @param mean_var_table : a data.frame containing three columns: group_name, mean, variance of expression per gene per grouping.
##'
##' @param normal_cell_expr : expression matrix of normal cells to guide the simulation. Should be total sum normalized.
##'
##' @param num_cells : number of cells to simulate
##'
##' @param common_dispersion: in rnbinom(size=1/common_dispersion).  Set to very small number to achieve Poisson.
##'
##' @return matrix containing simulated expression values.
##'
##' @keywords internal
##' @noRd
##'

.get_simulated_cell_matrix <- function(gene_means, mean_p0_table, num_cells, common_dispersion) {

    # should be working on the total sum count normalized data.
    # model the mean variance relationship

    ngenes = length(gene_means)

    dropout_logistic_params <- NULL
    if (! is.null(mean_p0_table)) {
        tryCatch (
            dropout_logistic_params <- .get_logistic_params(mean_p0_table),
            error=function(x) { cat(sprintf("(%s), zero inflation couldn't be estimated from data. Using just neg binom now\n", x)) }
        )
    }

    spike_cell_names = paste0('spike_cell_', seq_len(num_cells))

    sim_cell_matrix = matrix(rep(0,ngenes*num_cells), nrow=ngenes)
    rownames(sim_cell_matrix) = names(gene_means)
    colnames(sim_cell_matrix) = spike_cell_names

    sim_expr_vals <- function(gene_idx) {
        m = gene_means[gene_idx]
        return(.sim_expr_val(m, dropout_logistic_params, common_dispersion=common_dispersion))
    }

    for (i in seq_len(num_cells)) {
        newvals = sapply(seq_len(ngenes), FUN=sim_expr_vals)
        sim_cell_matrix[,i] = newvals
    }

    return(sim_cell_matrix)
}

##' @keywords internal
##' @noRd
##'

.sim_expr_val <- function(m,  dropout_logistic_params, common_dispersion=0.1, use_spline=TRUE) {
    
    # include drop-out prediction

    val = 0
    if (m > 0) {

        val = rnbinom(n=1, mu=m, size=1/common_dispersion)

        if ( (! is.null(dropout_logistic_params)) & val > 0) {

            if (use_spline) {
                dropout_prob <- predict(dropout_logistic_params$spline, log(val))$y[1]
            } else {
                dropout_prob <- .logistic_midpt_slope(x=log(val), midpt=dropout_logistic_params$midpt, slope=dropout_logistic_params$slope)
            }
            if (runif(1) <= dropout_prob) {
                ## a drop-out
                val = 0
            }
        }
    }

    return(val)
}



#' Computes probability of seeing a zero expr val as a function of the mean gene expression
#' The p(0 | mean_expr) is computed separately for each sample grouping.
#'
#' @keywords internal
#' @noRd
#'

.get_mean_vs_p0_table <- function(infercnv_obj) {

    group_indices = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)


    mean_vs_p0_table <- .get_mean_vs_p0_table_from_matrix(infercnv_obj@expr.data, group_indices)

    return(mean_vs_p0_table)

}

.get_mean_vs_p0_table_from_matrix <- function(expr.matrix, cell_groupings=NULL) {

    if (is.null(cell_groupings)) {
        ## use all cells as single group
        cell_groupings = list(allcells=seq(ncol(expr.matrix)))
    }

    mean_p0_table = NULL

    for (group_name in names(cell_groupings)) {
        #flog.info(sprintf("processing group: %s", group_name))
        expr.data = expr.matrix[, cell_groupings[[ group_name ]] ]

        group_mean_p0_table <- .get_mean_vs_p0_from_matrix(expr.data)
        group_mean_p0_table[[ 'group_name' ]] <- group_name

        if (is.null(mean_p0_table)) {
            mean_p0_table = group_mean_p0_table
        } else {
            mean_p0_table = rbind(mean_p0_table, group_mean_p0_table)
        }
    }

    return(mean_p0_table)
}

#' Computes probability of seeing a zero expr val as a function of the mean gene expression
#' based on the input expression  matrix.
#'
#' @keywords internal
#' @noRd
#'

.get_mean_vs_p0_from_matrix <- function(expr.data) {
    ncells = ncol(expr.data)
    m = rowMeans(expr.data)
    numZeros = apply(expr.data, 1, function(x) { sum(x==0) })

    pZero = numZeros/ncells

    mean_p0_table = data.frame(m=m, p0=pZero)

    return(mean_p0_table)
}


#'
#' Logistic function
#'
#' InferCNV note: Standard function here, but lifted from
#' Splatter (Zappia, Phipson, and Oshlack, 2017)
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0
#'
#' Implementation of the logistic function
#'
#' @param x value to apply the function to.
#' @param x0 midpoint parameter. Gives the centre of the function.
#' @param k shape parameter. Gives the slope of the function.
#'
#' @return Value of logistic function with given parameters
#'
#' @keywords internal
#' @noRd
#'
.logistic_midpt_slope <- function(x, midpt, slope) {
    1 / (1 + exp(-slope * (x - midpt)))
}



#' Given the mean, p0 table, fits the data to a logistic function to compute
#' the shape of the logistic distribution.
#'
#' @keywords internal
#' @noRd
#'

.get_logistic_params <- function(mean_p0_table) {

    mean_p0_table <- mean_p0_table[mean_p0_table$m > 0, ] # remove zeros, can't take log.

    x = log(mean_p0_table$m)
    y = mean_p0_table$p0

    df = data.frame(x,y)

    # write.table(df, "_logistic_params", quote=FALSE, sep="\t")  # debugging...

    logistic_params = list()

    tryCatch ( {
                
        fit <- nls(y ~ .logistic_midpt_slope(x, midpt = x0, slope = k),
                   data = df,
                   start = list(x0 = mean(x), k = -1)) # borrowed/updated from splatter
                   
                   
                   
                   logistic_params[[ 'midpt' ]] <- summary(fit)$coefficients["x0", "Estimate"]
                   logistic_params[[ 'slope' ]] <- summary(fit)$coefficients["k", "Estimate"]
        },
    
        error=function(x) { cat(sprintf("(%s), couldn't fit logistic, but no worries, going to use a spline\n", x)) }
    )
    
    
    ## also fit a spline
    s = smooth.spline(x, mean_p0_table$p0)
    logistic_params[[ 'spline' ]] = s


    return(logistic_params)
}



.estimate_common_dispersion <- function(expr.data) {

    ## estimate common disp from these data:
    ## creds to splatter
    design <- matrix(1, ncol(expr.data), 1)

    disps <- edgeR::estimateDisp(expr.data, design = design)

    common_dispersion = disps$common.dispersion

    flog.info(sprintf("-edgeR::estimateDisp() -> %g", common_dispersion))

    return(common_dispersion)

}



KS_plot <- function(title, tumor_expr, hspike_expr, names=NULL) {

    tumor_ecdf = ecdf(tumor_expr)
    hspike_ecdf = ecdf(hspike_expr)
    val_range = range(tumor_expr, hspike_expr)
    step = (val_range[2] - val_range[1])/100
    vals = seq(val_range[1], val_range[2], step)

    tumor_cdf = tumor_ecdf(vals)
    hspike_cdf = hspike_ecdf(vals)

    cdfs = data.frame(vals,
                      tumor_cdf,
                      hspike_cdf)

    if  ( (! is.null(names)) & length(names) == 2) {

        colnames(cdfs)[2] <- names[1]
        colnames(cdfs)[3] <- names[2]

        name1 <- names[1]
        name2 <- names[2]
    } else {
        name1 = 'tumor_cdf'
        name2 = 'hspike_cdf'
    }

    ks_point = which.max(abs(cdfs[,2] - cdfs[,3]))
    ks_point_info = cdfs[ks_point,]
    ##message("KS point info: ", paste(ks_point_info, collapse=', '))

    cdfs = cdfs %>% gather(name1, name2, key='type', value='cdf')

    p = ggplot(cdfs, aes_string(x=vals, y='cdf')) +
        geom_line(aes_string(color='type', linetype='type')) +
        geom_segment(aes(x=ks_point_info$vals,
                         y=ks_point_info[[name1]],
                         xend=ks_point_info$vals,
                         yend=ks_point_info[[name2]]), color='magenta', size=2) +
        ggtitle(title) + xlab("expr.val") + ylab("cdf")

    plot(p)

}



.mean_vs_p0_to_stats <- function(mean_vs_p0_table) {

    logm <- log(mean_vs_p0_table$m + 1)
    p0 <- mean_vs_p0_table$p0

    x_approx_mid <- median(logm[which(p0>0.2 & p0 < 0.8)])

    x <- logm
    y <- p0
    df <- data.frame(x,y)

    fit <- nls(y ~ .logistic(x, x0 = x0, k = k), data = df,
               start = list(x0 = x_approx_mid, k = -1))

    logistic_x <- x
    logistic_y <- predict(fit, newdata=x)
 
    ## also try fitting a spline
    spline.fit <- smooth.spline(x,y)
    spline.pts = predict(spline.fit, newdata=x)
 
    ret = list(logistic_x = logistic_x,
               logistic_y = logistic_y,
               spline_x <- spline.pts$x,
               spline_y <- spline.pts$y,
               spline.fit <- spline.fit,
               logistic.fit <- fit)
    

    return(ret)
}


