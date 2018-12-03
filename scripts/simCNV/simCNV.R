#!/usr/bin/env Rscript

set.seed(1234)

suppressPackageStartupMessages(library("argparse"))
library(edgeR)
library(infercnv)
library(dplyr)
library(Matrix)
library(stringr)

parser = ArgumentParser()

parser$add_argument("--normal_counts_matrix", required=T, nargs=1)
parser$add_argument("--CNV_spec", required=F, nargs='*', default=NULL)
parser$add_argument("--gene_order_file", required=T, nargs=1)
parser$add_argument("--num_normal_cells", type="integer", default=100, nargs=1)
parser$add_argument("--num_tumor_cells", type="integer", default=100, nargs='+')
parser$add_argument("--readRDS", action="store_true", default=FALSE, help="load matrix obj via readRDS")
parser$add_argument("--use_real_normals", action="store_true", default=FALSE, help="samples from normals instead of simulating them; note all tumor cells are simulated")
parser$add_argument("--no_run_infercnv", action="store_true", default=FALSE, help="do not run inferCNV automatically after simulating data")
parser$add_argument("--infercnv_cutoff", type='double', default=1.0, help="cutoff for infercnv filtering. (typically use 1 for smart-seq2, use 0.1 for 10x genomics)")
parser$add_argument("--copy_normal_to_simnormal", action="store_true")

args = parser$parse_args()


#' learn distribution parameters:
data = NULL
if (args$readRDS) {
    data = readRDS(args$normal_counts_matrix)
    ## ensure unique rownames
    rn = rownames(data)
    rn.uniq = unique(rn)
    if (length(rn) !=  length(rn.uniq)) {
        message("-removing non-unique rownames")
        m = match(rn.uniq, rn)
        data = data[m,]
    }
} else {
    data = read.table(args$normal_counts_matrix, header=T, row.names=1)
}

dim(data)

#' normalize first:
cs = colSums(data)
median_cs = median(cs)
data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
data <- data * median_cs

## get sim params
gene_means = rowMeans(data)

mean_p0_table <- infercnv:::.get_mean_vs_p0_from_matrix(data)
logistic_params <- infercnv:::.get_logistic_params(mean_p0_table)

## make simulated normals

normal_sim_matrix <- NULL
if (args$use_real_normals) {
    message("Sampling from normal cells")
    normal_sim_matrix <- data[, sample(x=1:ncol(data), size=args$num_normal_cells, replace=T)]
    normal_sim_matrix <- as.matrix(normal_sim_matrix)
    
} else {
    message("Simulating normal cells")
    normal_sim_matrix <- infercnv:::.get_simulated_cell_matrix(gene_means, mean_p0_table, args$num_normal_cells)
}
colnames(normal_sim_matrix) = paste0("normal_", 1:ncol(normal_sim_matrix))

cell_annots_df = data.frame("cells"=colnames(normal_sim_matrix), "class"="normal")

merged_matrix <- normal_sim_matrix

#' make simulated tumors:
gene_ordering = read.table(args$gene_order_file, header=F, row.names=NULL, stringsAsFactors=F)
colnames(gene_ordering) = c('gene', 'chr', 'lend', 'rend')


if (args$copy_normal_to_simnormal) {
    message("debugging - copying normal over to tumor")
    tumor_sim_matrix = normal_sim_matrix
    colnames(tumor_sim_matrix) = paste0("tumor_", 1:ncol(tumor_sim_matrix))

    cell_annots_df = rbind(cell_annots_df,
                           data.frame("cells"=colnames(tumor_sim_matrix), "class"="tumor") )

    merged_matrix = cbind(merged_matrix, tumor_sim_matrix)
    
    
} else {
    message("Simulating tumor cells")
 
    if (! is.null(args$CNV_spec)) {
        
        num_tumor_clusters = length(args$CNV_spec)
        if (length(args$num_tumor_cells) < length(args$CNV_spec)) {
            stop("Error, number of tumor cells vector is less than number of cnv_spec files")
        }
        
        for (cnv_file_idx in 1:length(args$CNV_spec)) {
            cnv_file = args$CNV_spec[cnv_file_idx]
            num_tumor_cells = args$num_tumor_cells[cnv_file_idx]
            
            cnv_info = read.table(cnv_file, header=F, row.names=NULL, stringsAsFactors=F)
            colnames(cnv_info) = c('chr', 'lend', 'rend', 'cnv')

            tumor_group_name = paste0("tumor_grp_", cnv_file_idx)
            message(sprintf("** Simulating tumor set: %s", tumor_group_name))

            gene_means_tumor = gene_means

            ## set amp/deletions for tumor group
            for (i in 1:nrow(cnv_info)) {
                r = unlist(cnv_info[i,,drop=T])
                message("simulating cnv for: ", paste(r, collapse="\t"))
                chrwant = r[1]
                chrlend = as.integer(r[2])
                chrrend = as.integer(r[3])
                cnv = as.double(r[4])
                
                chr_genes = gene_ordering %>% filter(chr==chrwant)
                
                if (chrlend > 0) {
                    ## if less than zero, using the whole chr
                    chr_genes = chr_genes %>% filter(lend >= chrlend & rend <= chrrend)
                }
                
                idx = which(names(gene_means) %in% chr_genes$gene)
                if (length(idx)==0) {
                    stop(sprintf("Error, found no genes on chr: %s", chrwant))
                }
                
                                
                gene_means_tumor[idx] = gene_means[idx] * cnv
            }

            ## sim the tumor matrix
            tumor_sim_matrix <- infercnv:::.get_simulated_cell_matrix(gene_means_tumor, mean_p0_table, num_tumor_cells)
            
            colnames(tumor_sim_matrix) = paste0(sprintf("%s_cell_", tumor_group_name), 1:ncol(tumor_sim_matrix))
            
            cell_annots_df = rbind(cell_annots_df,
                                   data.frame("cells"=colnames(tumor_sim_matrix), "class"="tumor") )
            
            merged_matrix = cbind(merged_matrix, tumor_sim_matrix)
                        
        }
    }
}



# write outputs
message("writing my.cell.annots")
write.table(cell_annots_df, file="my.cell.annots", quote=F, sep="\t", col.names=F, row.names=F)

message("writing merged.matrix")
## first, shuffle it

merged_matrix = merged_matrix[,sample(x=1:ncol(merged_matrix), size=ncol(merged_matrix),replace=F)]
write.table(merged_matrix, file="merged.matrix", quote=F, sep="\t")


if (! args$no_run_infercnv) {

    message("Running inferCNV")
    
    ## create the infercnv object
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix="merged.matrix",
                                        annotations_file="my.cell.annots",
                                        delim="\t",
                                        gene_order_file=args$gene_order_file,
                                        ref_group_names=c("normal"))
    
    out_dir=paste0("output_dir-", Sys.time())
    out_dir = str_replace(out_dir, " ", "_")
    ## perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=args$infercnv_cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir=out_dir, 
                                 cluster_by_groups=T, 
                                 plot_steps=T,
                                 include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
                                 )
                   
}

