#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
parser$add_argument("--no_scale_data", help="dont scale the data (ie. already scaled)", required=F, action='store_true', default=FALSE)
args = parser$parse_args()

library(infercnv)
library(ggplot2)
library(futile.logger)
library(HoneyBADGER)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

require(biomaRt) ## for gene coordinates
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = 'hsapiens_gene_ensembl',
                    host = "jul2015.archive.ensembl.org")

do_scale=TRUE
if (args$no_scale_data) {
    do_scale=FALSE
}


run_hbadger <- function(tumor_group_name, normal_matrix, tumor_matrix) {

    hb <- new('HoneyBADGER', name=tumor_group_name)

    ref_normal <- rowMeans(normal_matrix)
    
    hb$setGexpMats(tumor_matrix, ref_normal, mart.obj, filter=FALSE, scale=do_scale, verbose=TRUE)
    
    pdf(sprintf("%s-hb.pdf", tumor_group_name))

    hb$plotGexpProfile() ## initial visualization


    hb$setMvFit(verbose=TRUE)
    hb$setGexpDev(verbose=TRUE)
    hb$calcGexpCnvBoundaries(init=TRUE, verbose=FALSE)
    

    ## double check what CNVs were identified
    bgf <- hb$bound.genes.final
    genes <- hb$genes
    regions.genes <- range(genes[unlist(bgf)])
    
    print(regions.genes)

    if (length(regions.genes) == 0) {
        message("No cnv regions identified")
        return()
    }
    
    ## Indeed, our initial HMM has identified a number of candidate CNVs to test. We can now retest all identified CNVs on all cells to derive the final posterior probability of each CNV in each cell. We can cluster cells on these posterior probabilities and visualize them as a heatmap.
    
    hb$retestIdentifiedCnvs(retestBoundGenes = TRUE, retestBoundSnps = FALSE, verbose=FALSE)
    
    ## look at final results
    results <- hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
    print(head(results[,1:7]))
    write.table(results[,1:7], sprintf("%s-hb.cnvs.tsv", tumor_group_name), quote=F, sep="\t")
    
    
    ## visualize as heatmap 
    trees <- hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))
    
    ## order cells
    hc <- trees$hc
    order <- hc$labels[hc$order]
    ## plot all chromosomes
    hb$plotGexpProfile(cellOrder=order)
    
    
    ## plot just identified cnvs
    hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['amp']])
    
    hb$plotGexpProfile(cellOrder=order, region=hb$cnvs[['gene-based']][['del']])
    
    
}




normal_matrix = infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices), drop=F]

tumor_groups = infercnv_obj@observation_grouped_cell_indices

tumor_group_names = names(tumor_groups)
tumor_group_name = tumor_group_names[1] # for debugging
for (tumor_group_name in tumor_group_names) {
    tumor_grp_idx = tumor_groups[[tumor_group_name]]

    tumor_matrix = infercnv_obj@expr.data[,tumor_grp_idx]

    run_hbadger(tumor_group_name, normal_matrix, tumor_matrix)
}
