#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
library(infercnv)
library(tidyverse)


parser = ArgumentParser()
parser$add_argument("--counts_matrix", help="raw counts matrix file", required=TRUE, nargs=1)
parser$add_argument("--sim_method", help="simulation method: splatter, simple, meanvar", required=TRUE)
parser$add_argument("--include_dropout", default=FALSE, action='store_true', help='include dropout modeling')
args = parser$parse_args()


include.dropout = args$include_dropout


data = read.table(args$counts_matrix)
data = as.matrix(data)

orig.counts = data

if (! any(args$sim_method %in% c('splatter', 'simple', 'meanvar'))) {
    stop(sprintf("Error, not recognizing sim method: %s", args$sim_method))
}


#' normalize first:
cs = colSums(data)
median_cs = median(cs)
data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
data <- data * median_cs

gene_means <- rowMeans(data)

num_cells = ncol(data)

## sim the tumor matrix
sim_method = args$sim_method
if (sim_method == 'simple') {
    message('-using simple sim')

    mean_p0_table <- NULL
    if (include.dropout) {
        mean_p0_table <- infercnv:::.get_mean_vs_p0_from_matrix(data)
    }

    sim_matrix <- infercnv:::.get_simulated_cell_matrix(gene_means,
                                                        mean_p0_table=mean_p0_table,
                                                        num_cells=num_cells,
                                                        common_dispersion=0.1)
} else if (sim_method == 'splatter') {
    message('-using splatter sim')

    params <- infercnv:::.estimateSingleCellParamsSplatterScrape(orig.counts)

    params[['nCells']] <- num_cells
    params[['include.dropout']] <- include.dropout

    gene_means[gene_means == 0] <- 1e-3
    sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
    sim_matrix <- counts(sim_matrix)

} else if (sim_method == 'meanvar') {
    message('-using meanvar sim')
    ##tumor_sim_matrix <- infercnv:::.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(gene_means, data, args$num_tumor_cells)
    sim_matrix <- infercnv:::.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(gene_means, data, num_cells, include.dropout=include.dropout)

} else {
    stop(sprintf("not recognizing --sim_method: %s", args$sim_method))
}


## Plotting
if (include.dropout) {
    sim_method <- sprintf("%s-With_Dropout", sim_method)
} else {
    sim_method <- sprintf("%s-NO_Dropout", sim_method)
}

rownames(sim_matrix) <- names(gene_means)
colnames(sim_matrix) <- colnames(data)
sim_matrix_filename <- sprintf("sim.%s.counts.matrix", sim_method)
message("-writing matrix")
write.table(sim_matrix, sim_matrix_filename, quote=F, sep="\t")


message("-plotting QQ plot")
png(sprintf("sim_vs_orig_counts.%s.qqplots.png", sim_method))
qqplot(log(as.numeric(data)+1), log(as.numeric(sim_matrix)+1), main='orig vs. full sim')
abline(a=0,b=1,col='red')

message("-plotting KS plot")
png(sprintf("sim_vs_orig_counts.%s.KS.png", sim_method))
infercnv:::KS_plot(sprintf("KS, %s", sim_method), log(as.numeric(data)+1), log(as.numeric(sim_matrix)+1), names=c('orig',  sim_method))


