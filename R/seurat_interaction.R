##' @keywords internal
##' @noRd
##'
#make_seurat_from_infercnv_obj <- function(infercnv_obj) {
#    return(CreateSeuratObject(counts = infercnv_obj@count.data, project="infercnv", min.cells = 3, min.features = 200))
#}
#
##' @keywords internal
##' @noRd
##'
#make_seurat_from_infercnv <- function(infercnv_output_path) {
#    if (file.exists(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))) {
#        return(make_seurat_from_infercnv_obj(readRDS(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))))
#    }
#    else {
#        stop()
#    }
#}


#' @title add_to_seurat()
#'
#' @description Add meta.data about CNAs to a Seurat object from an infercnv_obj
#'
#' @param seurat_obj Seurat object to add meta.data to
#'
#' @param infercnv_output_path Path to the output folder of the infercnv run to use
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @export
#'
add_to_seurat <- function(seurat_obj, infercnv_output_path, top_n = 10) {
    lfiles <- list.files(infercnv_output_path, full.names = FALSE)
    
    if (!file.exists(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))) {
        flog.warn(sprintf("::Could not find \"run.final.infercnv_obj\" file at: %s"), paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))
        stop()
    }
    infercnv_obj = readRDS(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))
    
    if(!(setequal(row.names(seurat_obj@meta.data), colnames(infercnv_obj@expr.data)) ||
         setequal(colnames(seurat_obj@assays$RNA), colnames(infercnv_obj@expr.data)))) {
        flog.warn("::Cell names in Seurat object and infercnv results do not match")
        stop()
    }
    
    ## add check that data row/col names match seurat obj
    
    if (any(grep(lfiles, pattern="HMM_CNV_predictions.HMM.*.Pnorm_0.[0-9]+"))) {
        ###### states are 0/0.5/1/1.5/2
        scaling_factor = 1
        if (any(grep(lfiles, pattern="HMM_CNV_predictions.HMMi6.*.Pnorm_0.[0-9]+"))) {
            center_state = 1
        }
        else if (any(grep(lfiles, pattern="HMM_CNV_predictions.HMMi3.*.Pnorm_0.[0-9]+"))) {
            center_state = 1
        }
        else {
            flog.warn("::Found filtered HMM predictions output, but they do not match any known model type.")
            stop()
        }
        # sort to take lowest BayesProb if there are multiple
        regions = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles, pattern="HMM_CNV_predictions.HMMi6.*.Pnorm_0.[0-9]+.pred_cnv_regions.dat")])[1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles, pattern="HMM_CNV_predictions.HMMi6.*.Pnorm_0.[0-9]+.pred_cnv_genes.dat")])[1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        # from_pbayes()
    }
    else if (any(grep(lfiles, pattern = "12_HMM_preds"))) {
        ###### states are 1/2/3/4/5/6
        scaling_factor = 2
        if (any(grep(lfiles, pattern = "12_HMM_predHMMi6"))) {
            center_state = 3
        }
        else if (any(grep(lfiles, pattern = "12_HMM_predHMMi3"))) {
            center_state = 2
        }
        else {
            flog.warn("::Found HMM predictions output, but they do not match any known model type")
            stop()
        }
        regions = read.table(paste(infercnv_output_path, "12_HMM_preds.pred_cnv_regions.dat", sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, "12_HMM_preds.pred_cnv_genes.dat", sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        # from_hmm()
    }
    else {
        flog.warn(sprintf("::Could not find any HMM predictions outputs at: %s", infercnv_output_path))
        stop()
    }
    
    features_to_add <- .get_features(infercnv_obj = infercnv_obj, 
                                     regions = regions, 
                                     hmm_genes = hmm_genes, 
                                     center_state = center_state, 
                                     scaling_factor = scaling_factor, 
                                     top_n = top_n)
    
    for (lv in levels(infercnv_obj@gene_order$chr)) {
        seurat_obj@meta.data[[paste0("has_cnv_", lv)]] = features_to_add$feature_vector_chrs_has_cnv[[lv]]
        seurat_obj@meta.data[[paste0("has_loss_", lv)]] = features_to_add$feature_vector_chrs_has_loss[[lv]]
        seurat_obj@meta.data[[paste0("has_dupli_", lv)]] = features_to_add$feature_vector_chrs_has_dupli[[lv]]
        seurat_obj@meta.data[[paste0("proportion_cnv_", lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion[[lv]]
        seurat_obj@meta.data[[paste0("proportion_loss_", lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion[[lv]]
        seurat_obj@meta.data[[paste0("proportion_dupli_", lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion[[lv]]
        seurat_obj@meta.data[[paste0("proportion_scaled_cnv_", lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]]
        seurat_obj@meta.data[[paste0("proportion_scaled_loss_", lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion_scaled[[lv]]
        seurat_obj@meta.data[[paste0("proportion_scaled_dupli_", lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]]
    }
    
    for (n in names(features_to_add)[grep(names(features_to_add), pattern = "top_")] ) {
        seurat_obj@meta.data[[n]] = features_to_add[[n]]
    }
    
    return(seurat_obj)
}

#' @title .get_features()
#'
#' @description Get data from infercnv objects to add to Seurat meta.data
#'
#' @param infercnv_obj infercnv hmm object
#'
#' @param regions Table with predicted CNAs regions from the HMM model
#'
#' @param hmm_genes Table with predicted CNAs genes from the HMM model
#'
#' @param center_state Value that represents the neutral state in the HMM results.
#'
#' @param scaling_factor Factor to multiply divergence from center_state to get CNA amplitude.
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @return all_features A list of all the calculated meta.data to add.
#'
#' @keywords internal
#' @noRd
#'
.get_features <- function(infercnv_obj, regions, hmm_genes, center_state, scaling_factor, top_n) {
    
    chr_gene_count = table(infercnv_obj@gene_order$chr)
    
    # features templates for initialization
    double_feature_vector = vector(mode="double", length=ncol(infercnv_obj@expr.data))
    names(double_feature_vector) = colnames(infercnv_obj@expr.data)
    logical_feature_vector = vector(mode="logical", length=ncol(infercnv_obj@expr.data))
    names(logical_feature_vector) = colnames(infercnv_obj@expr.data)
    
    # initialize features lists
    all_features = c()
    all_features$feature_vector_chrs_has_cnv = c()
    all_features$feature_vector_chrs_has_dupli = c()
    all_features$feature_vector_chrs_has_loss = c()
    all_features$feature_vector_chrs_gene_cnv_proportion = c()
    all_features$feature_vector_chrs_gene_dupli_proportion = c()
    all_features$feature_vector_chrs_gene_loss_proportion = c()
    all_features$feature_vector_chrs_gene_cnv_proportion_scaled = c()
    all_features$feature_vector_chrs_gene_dupli_proportion_scaled = c()
    all_features$feature_vector_chrs_gene_loss_proportion_scaled = c()
    for (lv in levels(infercnv_obj@gene_order$chr)) {
        all_features$feature_vector_chrs_has_cnv[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_has_dupli[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_has_loss[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_gene_cnv_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_dupli_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_loss_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_loss_proportion_scaled[[lv]] = double_feature_vector
    }
    
    # map for top_n mapping
    subclust_name_to_clust = list()

    for (clust in names(infercnv_obj@tumor_subclusters$subclusters)) {
        for (subclust in names(infercnv_obj@tumor_subclusters$subclusters[[clust]])) {
            subclust_name = paste(clust, subclust, sep=".")
            subclust_name_to_clust[[subclust_name]] = c(clust, subclust)
            res = regions[regions$cell_group_name == subclust_name, , drop=FALSE]
            gres = hmm_genes[hmm_genes$cell_group_name == subclust_name, , drop=FALSE]
            if (nrow(res) > 0) {
                for (c in unique(res$chr)) {
                    all_features$feature_vector_chrs_has_cnv[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                    all_features$feature_vector_chrs_gene_cnv_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_cnv_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (sum(abs(gres[(which(gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                }
                
                sub_gres = gres[gres$state < center_state, ]
                for (c in unique(sub_gres$chr)) {
                    all_features$feature_vector_chrs_has_loss[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                    all_features$feature_vector_chrs_gene_loss_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_loss_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (abs(sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                }
                
                sub_gres = gres[gres$state > center_state, ]
                for (c in unique(sub_gres$chr)) {
                    all_features$feature_vector_chrs_has_dupli[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                    all_features$feature_vector_chrs_gene_dupli_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_dupli_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state) / (chr_gene_count[[c]] * scaling_factor))
                }
            }
        }
    }
    
    sorted_regions = sort(table(hmm_genes$gene_region_name), decreasing=TRUE)
    sorted_regions_loss = sort(table(hmm_genes$gene_region_name[hmm_genes$state < center_state]), decreasing=TRUE)
    sorted_regions_dupli = sort(table(hmm_genes$gene_region_name[hmm_genes$state > center_state]), decreasing=TRUE)
    
    top_n_cnv = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions, top_n = top_n)
    top_n_loss = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions = sorted_regions_loss, top_n = top_n)
    top_n_dupli = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions = sorted_regions_dupli, top_n = top_n)
    
    for (i in seq_along(top_n_cnv)) {
        feature_name = paste0("top_cnv_", i)
        all_features[[feature_name]] = logical_feature_vector
        
        for (subclust_name in top_n_cnv[[i]]$subclust_name) {
            clust = subclust_name_to_clust[[subclust_name]][1]
            subclust = subclust_name_to_clust[[subclust_name]][2]
            all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
        }
    }
    for (i in seq_along(top_n_loss)) {
        feature_name = paste0("top_loss_", i)
        all_features[[feature_name]] = logical_feature_vector
        
        for (subclust_name in top_n_loss[[i]]$subclust_name) {
            clust = subclust_name_to_clust[[subclust_name]][1]
            subclust = subclust_name_to_clust[[subclust_name]][2]
            all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
        }
    }
    for (i in seq_along(top_n_dupli)) {
        feature_name = paste0("top_dupli_", i)
        all_features[[feature_name]] = logical_feature_vector
        
        for (subclust_name in top_n_dupli[[i]]$subclust_name) {
            clust = subclust_name_to_clust[[subclust_name]][1]
            subclust = subclust_name_to_clust[[subclust_name]][2]
            all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
        }
    }
    
    return(all_features)
}


#' @title .get_top_n_regions()
#'
#' @description Get top n largest CNA regions in number of genes
#'
#' @param hmm_genes Table with predicted CNAs genes from the HMM model
#'
#' @param sorted_region List of regions sorted by size in number of genes for the CNA type desired (gain/loss/both)
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @return all_features A list of all the calculated meta.data to add.
#'
#' @keywords internal
#' @noRd
#'
.get_top_n_regions <- function(hmm_genes, sorted_regions, top_n) {
    j = 1
    previous_region_chr = -1
    previous_region_start = -1
    previous_region_end = -1
    top_regions = vector("list", top_n)
    
    for (i in seq_len(nrow(sorted_regions))) {
        genes_in_region = hmm_genes[which(hmm_genes$gene_region_name %in% names(sorted_regions[i])), ]
        region_chr = genes_in_region$chr[1]
        region_start = min(genes_in_region$start)
        region_end = max(genes_in_region$end)
        # check if the current region is the same as the previous one for a different subcluster or not
        # if it is, extend the previous assignment without increasing the count of found top hits
        if (region_chr == previous_region_chr && region_start == previous_region_start && region_end == previous_region_end) {
            top_regions[[j]]$subclust_names = c(top_regions[[j]]$subclust_names, genes_in_region$cell_group_name[1])
            top_regions[[j]]$regions_names = c(top_regions[[j]]$regions_names, genes_in_region$gene_region_name[1])
        }
        else {
            top_regions[[j]]$subclust_names = genes_in_region$cell_group_name[1]
            top_regions[[j]]$regions_names = genes_in_region$gene_region_name[1]
            previous_region_chr = region_chr
            previous_region_start = region_start
            previous_region_end = region_end
            j = j + 1
        }
        if (j == top_n + 1) {
            break
        }
    }
    
    if (j < top_n + 1) { # if less non unique regions than top_n
        top_regions = top_regions[1:j]
    }
    
    return(top_regions)
}



