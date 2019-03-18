#' Generated SmartSeq2 expression data with 10 normal cells and 10 tumor cells.
#' This is only to demonstrate how to use methods, not actual data to be used in an analysis.
#'
#' @format A data frame with 8252 rows (genes) and 20 columns (cells)
#'
#'
"data"

#' Generated classification for 10 normal cells and 10 tumor cells.
#'
#' @format A data frame with 20 rows (cells) and 1 columns (classification)
#'
#'
"annots"

#' Downsampled gene coordinates file from GrCh37
#'
#' @format A data frame with 10338 rows (genes) and 3 columns (chr, start, end)
#'
#'
"genes"

#' infercnv object result of the processing of run() in the example, to be used for other examples.
#'
#' @format An infercnv object
#'
#'
"infercnv_obj"

#' infercnv object result of the processing of run() in the HMM example, to be used for other examples.
#'
#' @format An infercnv object containing HMM predictions
#'
#'
"HMM_obj"

#' infercnv object result of the processing of inferCNVBayesNet in the example, to be used for other examples.
#'
#' @format An infercnv object containing posterior probability of CNV states
#'
#'
"mcmc_obj"
