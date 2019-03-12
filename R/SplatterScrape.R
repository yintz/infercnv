################################################################################################################
## Code here is primarily scraped from the Splatter package, leveraging just the pieces needed and further
## customized to our needs.
##
## Be sure to explore the original Splatter code as the source for these functions
## https://github.com/Oshlack/splatter
## and paper:
## Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA
## sequencing data. Genome Biology (2017).
##
## All code was 'scraped', 'lifted', whatever you want to call it, after discussions with the author Luke Zappia of the Splatter package
## as this form was easiest for integration of splatter sim methods customized to our needs.
## All attribution for single cell simulation methods is given to Zappia et al. and we're hugely thankful for being able to utilize it here.
#################################################################################################################


.estimateSingleCellParamsSplatterScrape <- function(counts,
                                                    include.dropout=FALSE,
                                                    use.spline.dropout.fit=FALSE # logistic is default.
                                                    ) {

    # scraped from splatter
    params = list()

    params[['include.dropout']] <- include.dropout
    params[['use.spline.dropout.fit']] <- use.spline.dropout.fit

    ## Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    params <- .splatEstMean(norm.counts, params)

    params <- .splatEstLib(counts, params)

    params <- .splatEstOutlier(norm.counts, params)

    params <- .splatEstBCV(counts, params)

    params <- .splatEstDropout(norm.counts, params)

    params[['nGenes']] <- nrow(counts)
    params[['nCells']] <- ncol(counts)

    print(params)

    return(params)
}


.splatEstMean <- function(norm.counts, params) {

    # library(fitdistrplus)

    means <- rowMeans(norm.counts)
    means <- means[means != 0]

    means <- .winsorize(means, q = 0.1)

    fit <- fitdistrplus::fitdist(means, "gamma", method = "mge",
                                 gof = "CvM")
    if (fit$convergence > 0) {
        warning("Fitting means using the Goodness of Fit method failed, ",
                "using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
    }

    params[['mean.shape']] <- unname(fit$estimate["shape"])
    params[['mean.rate']] <- unname(fit$estimate["rate"])

    return(params)
}

.winsorize <- function(x, q) {

    lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)

    if (diff(lohi) < 0) { lohi <- rev(lohi) }

    x[!is.na(x) & x < lohi[1]] <- lohi[1]
    x[!is.na(x) & x > lohi[2]] <- lohi[2]

    return(x)
}



.splatEstLib <- function(counts, params) {

    lib.sizes <- colSums(counts)

    if (length(lib.sizes) > 5000) {
        message("NOTE: More than 5000 cells provided. ",
                "5000 sampled library sizes will be used to test normality.")
        lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
    } else {
        lib.sizes.sampled <- lib.sizes
    }

    norm.test <- shapiro.test(lib.sizes.sampled)
    lib.norm <- norm.test$p.value > 0.2

    if (lib.norm) {
        fit <- fitdistrplus::fitdist(lib.sizes, "norm")
        lib.loc <- unname(fit$estimate["mean"])
        lib.scale <- unname(fit$estimate["sd"])
        message("NOTE: Library sizes have been found to be normally ",
                "distributed instead of log-normal. You may want to check ",
                "this is correct.")
    } else {
        fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
        lib.loc <- unname(fit$estimate["meanlog"])
        lib.scale <- unname(fit$estimate["sdlog"])
    }

    params[['lib.loc']] <- lib.loc
    params[['lib.scale']] <- lib.scale
    params[['lib.norm']] <- lib.norm

    return(params)
}


.splatEstOutlier <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)
    lmeans <- log(means)

    med <- median(lmeans)
    mad <- mad(lmeans)

    bound <- med + 2 * mad

    outs <- which(lmeans > bound)

    prob <- length(outs) / nrow(norm.counts)

    params[['out.prob']] <- prob

    if (length(outs) > 1) {
        facs <- means[outs] / median(means)
        fit <- fitdistrplus::fitdist(facs, "lnorm")

        params[['out.facLoc']] <- unname(fit$estimate["meanlog"])
        params[['out.facScale']] <- unname(fit$estimate["sdlog"])
    }

    return(params)
}


.splatEstBCV <- function(counts, params) {

    # Add dummy design matrix to avoid print statement
    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)


    ## linear adjustment to bcv is based on somulations as per splatter code documentation.
    params[['bcv.common']] <- 0.1 + 0.25 * disps$common.dispersion
    params[['bcv.df']] <- disps$prior.df

    return(params)
}


.splatEstDropout <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)

    x <- log(means)

    obs.zeros <- rowSums(norm.counts == 0)

    y <- obs.zeros / ncol(norm.counts)

    df <- data.frame(x, y)

    colnames(df) <- c('log_means', 'pct_zeros')
    #write.table(df, file="dropout.dat", quote=FALSE, sep="\t")
    #plot(df$log_means, df$pct_zeros)

    x_approx_mid <- median(x[which(y>0.2 & y < 0.8)]) # bhaas-added to avoid error: Error in nls(y ~ .logistic(x, x0 = x0, k = k), data = df, start = list(x0 = 0,  : singular gradient

    fit <- nls(y ~ .logistic(x, x0 = x0, k = k), data = df,
               start = list(x0 = x_approx_mid, k = -1))

    mid <- summary(fit)$coefficients["x0", "Estimate"]
    shape <- summary(fit)$coefficients["k", "Estimate"]

    #points(x, predict(fit, newdata=x), col='green')

    params[['dropout.mid']] <- mid
    params[['dropout.shape']] <- shape


    ## also try fitting a spline
    spline.fit <- smooth.spline(x,y)
    params[['dropout.spline.fit']] <- spline.fit
    spline.pts = predict(spline.fit, newdata=x)
    #points(spline.pts$x, spline.pts$y, col='magenta')
    #legend('topright', c('logistic', 'spline'), col=c('green', 'magenta'), pch=1)

    
    return(params)
}

.logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}


#####################################
### End of Splat Estimation routines
#####################################
## Beginning of Splat Simulation routines
#########################################

.simulateSingleCellCountsMatrixSplatterScrape <- function(params,
                                                          use.genes.means=NULL
                                                          ) {

    if ( (! is.null(use.genes.means)) && length(use.genes.means) != params[['nGenes']]) {
        stop("Error, use.genes.means provided but not matching the params nGenes count")
    }

    # library(SingleCellExperiment)

    ## Get the parameters we are going to use
    nCells <- params[["nCells"]]
    nGenes <- params[["nGenes"]]

    # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    ## Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features,
                                colData = cells,
                                metadata = list(Params = params))

    message("Simulating library sizes...")
    sim <- .splatSimLibSizes(sim, params)

    message("Simulating gene means...")
    sim <- .splatSimGeneMeans(sim, params, use.genes.means)

    sim <- .splatSimBatchCellMeans(sim, params)

    sim <- .splatSimSingleCellMeans(sim, params)

    message("Simulating BCV...")
    sim <- .splatSimBCVMeans(sim, params)

    message("Simulating counts...")
    sim <- .splatSimTrueCounts(sim, params)

    message("Simulating dropout (if needed)...")
    sim <- .splatSimDropout(sim, params)

    return(sim)
}

.splatSimLibSizes <- function(sim, params) {

    nCells <- params[["nCells"]]
    lib.loc <- params[["lib.loc"]]
    lib.scale <- params[["lib.scale"]]
    lib.norm <- params[["lib.norm"]]

    if (lib.norm) {
        exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
        min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
        exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
    } else {
        exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    }

    colData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}


.splatSimGeneMeans <- function(sim, params, use.genes.means) {

    nGenes <- params[["nGenes"]]
    mean.shape <- params[["mean.shape"]]
    mean.rate <- params[["mean.rate"]]
    out.prob <- params[["out.prob"]]
    out.facLoc <- params[["out.facLoc"]]
    out.facScale <- params[["out.facScale"]]

    if (! is.null(use.genes.means)) {
        base.means.gene <- use.genes.means
    } else {
        ## Simulate base gene means
        base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
    }
    ## Add expression outliers
    outlier.facs <- .getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                     out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene

    return(sim)
}

.getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {

    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
    facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected

    return(factors)
}

.splatSimBatchCellMeans <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    gene.means <- rowData(sim)$GeneMean

    nCells <- params[["nCells"]]
    nGenes <- params[["nGenes"]]

    batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)

    batch.means.cell <- batch.facs.cell * gene.means

    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- gene.names
    assays(sim)$BatchCellMeans <- batch.means.cell

    return(sim)
}




.splatSimSingleCellMeans <- function(sim, params) {

    nCells <- params[["nCells"]]
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    cell.means.gene <- batch.means.cell
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    assays(sim)$CellMeans <- base.means.cell # default, updated under .splatSimBCVMeans()

    return(sim)
}


.splatSimBCVMeans <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- params[["nGenes"]]
    nCells <- params[["nCells"]]
    bcv.common <- params[["bcv.common"]]
    bcv.df <- params[["bcv.df"]]
    base.means.cell <- assays(sim)$BaseCellMeans

    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    } else {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(base.means.cell)))
    }

    means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                                scale = base.means.cell * (bcv ^ 2)),
                         nrow = nGenes, ncol = nCells)

    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    assays(sim)$BCV <- bcv
    assays(sim)$CellMeans <- means.cell

    return(sim)

}

.splatSimTrueCounts <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- params[["nGenes"]]
    nCells <- params[["nCells"]]
    cell.means <- assays(sim)$CellMeans

    true.counts <- matrix(rpois(nGenes * nCells, lambda = cell.means),
                          nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    assays(sim)$TrueCounts <- true.counts

    return(sim)
}

.splatSimDropout <- function(sim, params) {

    include.dropout <- params[["include.dropout"]]
    true.counts <- assays(sim)$TrueCounts
    dropout.mid <- params[["dropout.mid"]]
    dropout.shape <- params[["dropout.shape"]]
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nCells <- params[["nCells"]]
    nGenes <- params[["nGenes"]]
    nBatches <- params[["nBatches"]]
    nGroups <- params[["nGroups"]]
    cell.means <- assays(sim)$CellMeans
    dropout.spline.fit <- params[['dropout.spline.fit']]

    if (include.dropout) {

        if ( params[['use.spline.dropout.fit']] ) {
            ## Generate probabilites based on expression
            drop.prob <- sapply(seq_len(nCells), function(idx) {
                eta <- log(cell.means[, idx])
                pvals <- predict(dropout.spline.fit, eta)$y
                pvals[is.na(pvals)] <- 0
                pvals[pvals<0] <- 0
                pvals[pvals>1] <- 1
                return(pvals)
            })


        } else {
            # using logistic
            dropout.mid <- rep(dropout.mid, nCells)
            dropout.shape <- rep(dropout.shape, nCells)

            ## Generate probabilites based on expression
            drop.prob <- sapply(seq_len(nCells), function(idx) {
                eta <- log(cell.means[, idx])
                return(.logistic(eta, x0 = dropout.mid[idx], k = dropout.shape[idx]))
            })
        }

        print(drop.prob)

        # Decide which counts to keep
        keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob),
                       nrow = nGenes, ncol = nCells)

        counts <- true.counts * keep

        colnames(drop.prob) <- cell.names
        rownames(drop.prob) <- gene.names
        colnames(keep) <- cell.names
        rownames(keep) <- gene.names

        assays(sim)$DropProb <- drop.prob
        assays(sim)$Dropout <- !keep
    } else {
        counts <- true.counts
    }

    BiocGenerics::counts(sim) <- counts

    return(sim)
}
