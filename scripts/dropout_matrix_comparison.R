#!/usr/bin/env Rscript

set.seed(1234)

suppressPackageStartupMessages(library("argparse"))

library(tidyverse)
library(infercnv)

parser = ArgumentParser()
parser$add_argument("--matrix1", required=T, nargs=1)
parser$add_argument("--matrix2", required=T, nargs=1)
parser$add_argument("--output", required=T, nargs=1, help="output filename pdf")

args = parser$parse_args()


#' learn distribution parameters:
data1 = as.matrix(read.table(args$matrix1, header=T, row.names=1))
data2 = as.matrix(read.table(args$matrix2, header=T, row.names=1))


pdf(args$output)


data1.mean_vs_p0 <- infercnv:::.get_mean_vs_p0_from_matrix(data1)
data2.mean_vs_p0 <- infercnv:::.get_mean_vs_p0_from_matrix(data2)

plot_mean_vs_p0_with_data <- function(title='title', mean_vs_p0_table) {

    logm <- log(mean_vs_p0_table$m + 1)
    p0 <- mean_vs_p0_table$p0

    plot(logm, p0, pch='.', main=title)

    x_approx_mid <- median(logm[which(p0>0.2 & p0 < 0.8)])

    x <- logm
    y <- p0
    df <- data.frame(x,y)

    fit <- nls(y ~ infercnv:::.logistic(x, x0 = x0, k = k), data = df,
               start = list(x0 = x_approx_mid, k = -1))

    logistic_x <- x
    logistic_y <- predict(fit, newdata=x)
    points(x, logistic_y, col='green')

    ## also try fitting a spline
    spline.fit <- smooth.spline(x,y)
    spline.pts = predict(spline.fit, newdata=x)
    points(spline.pts$x, spline.pts$y, col='magenta')
    legend('topright', c('logistic', 'spline'), col=c('green', 'magenta'), pch=1)

    ret = list(logistic_x = logistic_x,
               logistic_y = logistic_y,
               spline_x <- spline.pts$x,
               spline_y <- spline.pts$y)


    return(ret)
}


p1 <- plot_mean_vs_p0_with_data(args$matrix1, data1.mean_vs_p0)
p2 <- plot_mean_vs_p0_with_data(args$matrix2, data2.mean_vs_p0)


## plot both logistics in a single plot
plot(p1$logistic_x, p1$logistic_y, col='blue')
points(p2$logistic_x, p2$logistic_y, col='magenta')
legend('topright', c(args$matrix1, args$matrix2), col=c('blue', 'magenta'), pch=1)



