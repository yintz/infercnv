
#' @title mask_non_DE_genes_via_Seurat
#'
#' Masks expression data values for those genes defined as *not* differentially expressed.
#' Seurat FindMarkers() is used to define differentially expressed genes.
#' All tumor types are each compared pairwise to the normal cell type groups.
#' Any genes not found differentially expressed according to the defined parameters in any
#' of the tumor/normal comparisons are set to the center_val, effectively masking them.
#'
#' @param infercnv_obj
#'
#' @param thresh.use  from Seurat::FindMarkers : "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing thresh.use speeds up the function, but can miss weaker signals."
#'
#' @param min.pct   from Seurat::FindMarkers : "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed."
#'
#' @param p_val_threshold maximum allowed p-value
#'
#' @param test.use from Seurat::FindMarkers : "Denotes which test to use. Seurat currently implements "bimod" (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'. The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution"
#'
#' @param center_val value to set the non-DE gene expression values to in the matrix.  By default, sets to the mean of all input expression values.
#'
#' @param normalize_factor normalization factor used by Seurat in total sum count scaling. Typically set to 1e4 or 1e5. Default: estimates from the input data.
#' 


mask_non_DE_genes_via_Seurat <- function(infercnv_obj,
                                         thresh.use=0.25,
                                         min.pct=0.50,
                                         p_val_thresh = 0.05,
                                         test.use='bimod',
                                         center_val = mean(infercnv_obj@expr.data),
                                         normalize_factor=compute_normalization_factor(infercnv_obj) ) {
    
    
  

    all_DE_results = get_DE_genes_via_Seurat(infercnv_obj,
                                             thresh.use=thresh.use,
                                             min.pct=min.pct,
                                             p_val_thresh=p_val_thresh,
                                             test.use=test.use,
                                             normalize_factor=normalize_factor)


    infercnv_obj <- .mask_DE_genes(infercnv_obj, all_DE_results, mask_val=center_val)

    return(infercnv_obj)

}


#' @title .mask_DE_genes()
#' private function that does the actual masking of expression values in the matrix
#' according to the specified mask value
#'
#' @param infercnv_obj
#'
#' @param all_DE_results  DE results list with structure:
#'                           all_DE_results[[ tumor_type ]] = list(tumor=tumor_type, de_genes=c(geneA, geneB, ...))
#'
#' @param mask_val float value to assign to genes that are in the mask set (complement of the DE gene set)
#'
#' @return infercnv_obj
#' 

.mask_DE_genes <- function(infercnv_obj, all_DE_results, mask_val) {

    
    all_DE_genes_matrix = matrix(data=FALSE, nrow=nrow(infercnv_obj@expr.data),
                                 ncol=ncol(infercnv_obj@expr.data),
                                 dimnames = list(rownames(infercnv_obj@expr.data),
                                                 colnames(infercnv_obj@expr.data) ) )



    all_DE_genes_matrix[,] = FALSE
    
    ## turn on all normal genes
    all_DE_genes_matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)] = TRUE
    
    for (DE_results in all_DE_results) {
        tumor_type = DE_results$tumor
        genes = DE_results$de_genes
        
    
        all_DE_genes_matrix[rownames(all_DE_genes_matrix) %in% genes,
                            infercnv_obj@observation_grouped_cell_indices[[ tumor_type ]] ] = TRUE
    
    }
        
    infercnv_obj@expr.data[ ! all_DE_genes_matrix ] = mask_val
    
    return(infercnv_obj)
    

}



#' @title seurat_from_infercnv_obj()
#'
#' Generate a Seurat object from an infercnv_obj
#'
#' @param infercnv_obj
#'
#' @param normalize_factor
#'
#' 


seurat_from_infercnv_obj <- function(infercnv_obj, normalize_factor) {

    seurat_obj = CreateSeuratObject(raw.data = infercnv_obj@count.data, project="infcnv")
    
    seurat_obj = NormalizeData(object=seurat_obj, normalization.method = "LogNormalize",
                               scale.factor = normalize_factor)

    ## store cell type info as metadata

    all_cell_classes = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

    cell_classes = names(all_cell_classes)

    cell_column_annotations = rep(NA, ncol(infercnv_obj@expr.data))

    for (cell_class in cell_classes) {
        cell_column_annotations[ all_cell_classes[[cell_class]] ] <- cell_class
    }
    cell_column_annotations = as.factor(cell_column_annotations)
    names(cell_column_annotations) = colnames(seurat_obj@data)


    patient.ident = seurat_obj@ident

    seurat_obj = StashIdent(seurat_obj, save.name='patient.ident')

    seurat_obj@ident <- as.factor(cell_column_annotations)

    return(seurat_obj)


}



#' @title get_DE_genes_via_Seurat()
#'
#' 
#' Seurat FindMarkers() is used to define differentially expressed genes.
#' All tumor types are each compared pairwise to the normal cell type groups.
#' Any genes not found differentially expressed according to the defined parameters in any
#' of the tumor/normal comparisons are set to the center_val, effectively masking them.
#'
#' @param infercnv_obj
#'
#' @param thresh.use  from Seurat::FindMarkers : "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing thresh.use speeds up the function, but can miss weaker signals."
#'
#' @param min.pct   from Seurat::FindMarkers : "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed."
#'
#' @param p_val_threshold maximum allowed p-value
#'
#' @param test.use from Seurat::FindMarkers : "Denotes which test to use. Seurat currently implements "bimod" (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'. The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution"
#'
#' @param normalize_factor normalization factor used by Seurat in total sum count scaling. Typically set to 1e4 or 1e5. Default: estimates from the input data.
#'
#' @return DE_results_list, structured like so:
#'
#'             all_DE_results[[condition_pair]] = list(tumor=tumor_type,
#'                                                     normal=normal_type,
#'                                                     markers=markers,
#'                                                     de_genes=genes)
#' 
 

get_DE_genes_via_Seurat <- function(infercnv_obj,
                                    thresh.use=0.25,
                                    min.pct=0.50,
                                    p_val_thresh = 0.05,
                                    test.use='bimod',
                                    normalize_factor=1e5 ) {
    

    seurat_obj = seurat_from_infercnv_obj(infercnv_obj, normalize_factor)
    
    ## Find DE genes by comparing the mutant types to normal types
    normal_types = names(infercnv_obj@reference_grouped_cell_indices)
    tumor_types = names(infercnv_obj@observation_grouped_cell_indices)

    all_DE_results = list()
    
    ## turn on only DE genes in tumors
    for (tumor_type in tumor_types) {
        for (normal_type in normal_types) {
            flog.info(sprintf("Finding DE genes between %s and %s", tumor_type, normal_type))
            markers = FindMarkers(seurat_obj,
                                  ident.1=tumor_type,
                                  ident.2=normal_type,
                                  thresh.use=thresh.use,
                                  min.pct=min.pct,
                                  only.pos=FALSE,
                                  test.use=test.use)
            

            
            genes = rownames(markers[markers$p_val<0.05,])
            flog.info(sprintf("Found %d genes as DE", length(genes)))

            condition_pair = paste(tumor_type, normal_type, sep=",")
            
            all_DE_results[[condition_pair]] = list(tumor=tumor_type,
                                                    normal=normal_type,
                                                    markers=markers,
                                                    de_genes=genes)
            
            
        }
    }
    
    return(all_DE_results)

}



#'
#' @title mask_non_DE_genes_basic()
#'
#' Mask gene expression in infercnv_obj based on expression values that are found
#' statistically significantly different between tumor/normal based on specified
#' statistical test.
#'
#' Uses the expression data stored in infercnv_obj@expr.data
#'
#' @param infercnv_obj
#'
#' @param p_val_thresh  p-value threshold for assigning DE statistically significant 
#'
#' @param test.use statistical test to use.  (default: "t" for t-test.) **note, other tests still need to be implemented here. basic plumbing provided.
#'
#' @param center_val value to assign to those genes that are not found to be statistically DE.
#'
#' @return infercnv_obj
#' 


mask_non_DE_genes_basic <- function(infercnv_obj,
                                    p_val_thresh = 0.05,
                                    test.use="t",
                                    center_val=mean(infercnv_obj@expr.data) ) {
    

    all_DE_results = get_DE_genes_basic(infercnv_obj,
                                        p_val_thresh=p_val_thresh,
                                        test.use=test.use)

    infercnv_obj <- .mask_DE_genes(infercnv_obj,
                                   all_DE_results,
                                   mask_val=center_val)

    return(infercnv_obj)
}



#'
#' @title get_DE_genes_basic
#'
#' Retrieves genes identified as significantly DE based on the corresponding statistical test as applied
#' to the infercnv_obj@expr.data with tumor/normal comparison
#'
#' @param p_val_thresh maximum p-value for defining statistical signficance
#'
#' @param test.use  statistical test to use. (default: "t" for t-test)
#' 
#' @param return             all_DE_results[[condition_pair]] = list(tumor=tumor_type,
#'                                                                   normal=normal_type,
#'                                                                   pvals=pvals,
#'                                                                   de_genes=genes)
#'


get_DE_genes_basic <- function(infercnv_obj,
                               p_val_thresh = 0.05,
                               test.use="t" # other options:  (TBD)
                               ) {
    

    
    ## Find DE genes by comparing the mutant types to normal types
    normal_types = names(infercnv_obj@reference_grouped_cell_indices)
    tumor_types = names(infercnv_obj@observation_grouped_cell_indices)

    all_DE_results = list()

    statfxns = list()
    statfxns[[ "t" ]] <- function(x, idx1, idx2) {
        vals1 = x[idx1]
        vals2 = x[idx2]

        # useful way of handling tests that may fail:
        # https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html
        res = try(t.test(vals1, vals2), silent=TRUE)
        
        if (is(res, "try-error")) return(NA) else return(res$p.value)
        
    }
    
    statfxn = statfxns[[ test.use ]]
    
    ## turn on only DE genes in tumors
    for (tumor_type in tumor_types) {

        tumor_indices = infercnv_obj@observation_grouped_cell_indices[[ tumor_type ]] 
        
        for (normal_type in normal_types) {
            flog.info(sprintf("Finding DE genes between %s and %s", tumor_type, normal_type))

            normal_indices = infercnv_obj@reference_grouped_cell_indices[[ normal_type ]]


            pvals = apply(infercnv_obj@expr.data, 1, statfxn, idx1=normal_indices, idx2=tumor_indices)
            pvals = unlist(pvals)
            
            names(pvals) = rownames(infercnv_obj@expr.data)

            genes = names(pvals)[pvals<p_val_thresh]
            
            flog.info(sprintf("Found %d genes as DE", length(genes)))

            condition_pair = paste(tumor_type, normal_type, sep=",")
            
            all_DE_results[[condition_pair]] = list(tumor=tumor_type,
                                                    normal=normal_type,
                                                    pvals=pvals,
                                                    de_genes=genes)
            
            
        }
    }
    
    return(all_DE_results)

}

