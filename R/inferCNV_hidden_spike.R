


.get_hspike_chr_info <- function(num_genes_each, num_total) {

    num_remaining = num_total - 10*num_genes_each
    if (num_remaining < num_genes_each) {
        num_remaining = num_genes_each
    }

    ## design for fake chr
    chr_info = list(list(name='chrA',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_0',
                         cnv=0.01,
                         ngenes=num_genes_each),
                    list(name='chr_B',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_0pt5',
                         cnv=0.5,
                         ngenes=num_genes_each),
                    list(name='chr_C',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_1pt5',
                         cnv=1.5,
                         ngenes=num_genes_each),
                    list(name='chr_D',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_2pt0',
                         cnv=2.0,
                         ngenes=num_genes_each),
                    list(name='chr_E',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_3pt0',
                         cnv=3,
                         ngenes=num_genes_each),
                    list(name='chr_F',
                         cnv=1,
                         ngenes=num_remaining)
                    )

    return(chr_info)

}

.build_and_add_hspike <- function(infercnv_obj, sim_method=c('splatter', 'simple', 'meanvar')) {

    sim_method = match.arg(sim_method)

    flog.info("Adding h-spike")

    normal_cells_idx = infercnv::get_reference_grouped_cell_indices(infercnv_obj)

    params = .estimateSingleCellParamsSplatterScrape(counts=infercnv_obj@count.data[,normal_cells_idx])


    ## build a fake genome with fake chromosomes, alternate between 'normal' and 'variable' regions.

    num_cells = 100
    num_genes_per_chr = 400

    num_total_genes = nrow(infercnv_obj@expr.data)

    chr_info <- .get_hspike_chr_info(num_genes_per_chr, num_total_genes)

    gene_order = do.call(rbind, lapply(chr_info, function(x) { data.frame(chr=x$name, start=1:x$ngenes, end=1:x$ngenes) }))
    num_genes = nrow(gene_order)
    rownames(gene_order) <- paste0("gene_", 1:num_genes)

    ## sample gene info from the normal data
    normal_cells_expr = infercnv_obj@expr.data[,normal_cells_idx]
    gene_means = rowMeans(normal_cells_expr)
    ## remove zeros (might not be zero in the non-normal cells)
    gene_means = gene_means[gene_means>0]

    gene_means = sample(x=gene_means, size=num_genes, replace=T)

    names(gene_means) = rownames(gene_order)

    mean_p0_table <- NULL
    if (sim_method == 'simple') {
        mean_p0_table <- infercnv:::.get_mean_vs_p0_from_matrix(normal_cells_expr)
    }
    
    ## simulate normals:
    params[["nGenes"]] <- num_genes
    params[["nCells"]] <- num_cells

    if (sim_method == 'splatter') {
        sim.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(params, use.genes.means=gene_means)

        sim_normal_matrix = counts(sim.scExpObj)
    } else if (sim_method == 'simple') {
        ## using simple
        
        sim_normal_matrix <- infercnv:::.get_simulated_cell_matrix(gene_means,
                                                                   mean_p0_table=mean_p0_table,
                                                                   num_cells=num_cells,
                                                                   common_dispersion=0.1)
    } else if (sim_method == 'meanvar') {
        ## using mean,var trend
        sim_normal_matrix <- infercnv:::.get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj, gene_means, num_cells, include.dropout=TRUE)

    }

    colnames(sim_normal_matrix) = paste0('simnorm_cell_', 1:num_cells)
    rownames(sim_normal_matrix) = rownames(gene_order)

    ## apply spike-in multiplier vec
    hspike_gene_means = gene_means
    for (info in chr_info) {
        chr_name = info$name
        cnv = info$cnv
        if (cnv != 1) {
            gene_idx = which(gene_order$chr == chr_name)
            hspike_gene_means[gene_idx] =  hspike_gene_means[gene_idx] * cnv
        }
    }


    if (sim_method == 'splatter') {
        sim_spiked_cnv.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(params,
                                                                            use.genes.means=hspike_gene_means)
        sim_spiked_cnv_matrix = counts(sim_spiked_cnv.scExpObj)
    } else if (sim_method == 'simple') {
        ## using simple
        sim_spiked_cnv_matrix <- infercnv:::.get_simulated_cell_matrix(hspike_gene_means,
                                                                       mean_p0_table=mean_p0_table,
                                                                       num_cells=num_cells,
                                                                       common_dispersion=0.1)
    } else if (sim_method == 'meanvar') {

        sim_spiked_cnv_matrix <- infercnv:::.get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj, hspike_gene_means, num_cells, include.dropout=TRUE)
    }

    colnames(sim_spiked_cnv_matrix) = paste0('spike_cell_', 1:num_cells)
    rownames(sim_spiked_cnv_matrix) = rownames(gene_order)

    sim.counts.matrix = cbind(sim_normal_matrix, sim_spiked_cnv_matrix)
    reference_grouped_cell_indices = list('simnormal'=1:num_cells)
    observation_grouped_cell_indices = list('SPIKE'=(num_cells+1):(2*num_cells))

    .hspike <- new( Class="infercnv",
                   expr.data=sim.counts.matrix,
                   count.data=sim.counts.matrix,
                   gene_order=gene_order,
                   reference_grouped_cell_indices=reference_grouped_cell_indices,
                   observation_grouped_cell_indices=observation_grouped_cell_indices)

    validate_infercnv_obj(.hspike)

    .hspike <- normalize_counts_by_seq_depth(.hspike, median(colSums(normal_cells_expr))) # make same target counts/cell as normals

    infercnv_obj@.hspike <- .hspike

    return(infercnv_obj)
}


##
.sim_foreground <- function(infercnv_obj, sim_method) {

    flog.info("## simulating foreground")

    expr.matrix <- infercnv_obj@expr.data

    samples <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)


    normal_data <- expr.matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)]
    target_total_counts <- median(colSums(normal_data))

    params <- NULL
    if (sim_method == 'splatter') {
        params <- infercnv:::.estimateSingleCellParamsSplatterScrape(infercnv_obj@count.data[, unlist(infercnv_obj@reference_grouped_cell_indices)])
    }


    mean_p0_table <- NULL
    if (sim_method == 'simple') {
        mean_p0_table <- infercnv:::.get_mean_vs_p0_from_matrix(normal_data)
    }
    
    for (sample_name in names(samples)) {

        cell_idx = samples[[ sample_name ]]

        sample_expr_data = expr.matrix[, cell_idx]
        gene_means = rowMeans(sample_expr_data)
        gene_means[gene_means==0] <- 1e-3 #avoid zeros, breaks splatter sim

        num_cells = ncol(sample_expr_data)

        ## sim the tumor matrix
        if (sim_method == 'simple') {
            
            sim_matrix <- infercnv:::.get_simulated_cell_matrix(gene_means,
                                                                mean_p0_table=mean_p0_table,
                                                                num_cells=num_cells,
                                                                common_dispersion=0.1)
        } else if (sim_method == 'splatter') {
            params[['nCells']] <- num_cells
            sim_matrix <- infercnv:::.simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
            sim_matrix <- counts(sim_matrix)
        } else if (sim_method == 'meanvar') {
            sim_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj,
                                                                         gene_means,
                                                                         num_cells,
                                                                         include.dropout=TRUE)
        } else {
            stop(sprintf("not recognizing --sim_method: %s", args$sim_method))
        }

        expr.matrix[, cell_idx] <- sim_matrix

    }

    infercnv_obj@expr.data <- expr.matrix

    infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj, target_total_counts)

    return(infercnv_obj)
}

