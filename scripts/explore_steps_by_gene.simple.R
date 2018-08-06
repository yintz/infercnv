#!/usr/bin/env Rscript

args<-commandArgs(TRUE)

if (length(args) == 0) {
    stop("\n\n\tUsage:\texplore_steps_by_gene.simple.R gene_name\n\n");
}

gene_name = args[1]

options(error = function() traceback(2))

main = function() {

    data_bundles = load_all_data()

    pdf_filename = paste("by_step.", gene_name, ".pdf", sep="")
    pdf(pdf_filename)
    
    make_plots(gene_name, data_bundles)

    quit(save = "no", status = 0, runLast = FALSE)

}


load_data = function(matrix_filename) {

    message(paste("loading:", matrix_filename))
    
    data = read.table(matrix_filename, header=T, row.names=1, sep='\t')
    
    ref = read.table("references.txt")
    ref_cells = colnames(ref)
    
    obs = read.table("observations.txt")
    obs_cells = rownames(obs)
    
    ref_matrix = data[,colnames(data) %in% ref_cells]
    obs_matrix = data[,colnames(data) %in% obs_cells]
    

    data_bundle = list()

    data_bundle[['filename']] = matrix_filename
    
    data_bundle[['ref_cells']] = ref_cells
    data_bundle[['obs_cells']] = obs_cells

    data_bundle[['ref_matrix']] = ref_matrix
    data_bundle[['obs_matrix']] = obs_matrix

    return(data_bundle)
       
}

plot_gene_dist = function(gene_name, data_bundle, drop_zeros=F) {

    ref_gene_expr = na.omit(as.numeric(data_bundle$ref_matrix[gene_name,]))
    obs_gene_expr = na.omit(as.numeric(data_bundle$obs_matrix[gene_name,,drop=T]))

    if (drop_zeros) {
        ref_gene_dens = density(ref_gene_expr[ref_gene_expr!=0])
        obs_gene_dens = density(obs_gene_expr[obs_gene_expr!=0])
    }
    else {
        ref_gene_dens = density(ref_gene_expr)
        obs_gene_dens = density(obs_gene_expr)
    }
    xrange = range(ref_gene_dens$x, obs_gene_dens$x)
    yrange = range(ref_gene_dens$y, obs_gene_dens$y)

    plot(ref_gene_dens, xlim=xrange, ylim=yrange, t='l', col='blue',
         main=data_bundle$filename)
    points(obs_gene_dens, t='l', col='green')
}



make_plots = function(gene_name, data_bundles) {

    #dev.new()
    par(mfrow=c(3,3))

    count = 0
    for (data_bundle in data_bundles) {
        count = count + 1
        remove_zeros = F
        if (count < 3) { remove_zeros = T }
        plot_gene_dist(gene_name, data_bundle, remove_zeros)
    }

}

load_all_data = function() {
    
    S01 = load_data("01_incoming_data.pdf.txt")
    S03 = load_data("03_reduced_by_cutoff.pdf.txt")
    S04 = load_data("04_center_with_threshold.pdf.txt")
    
    S05 = load_data("05_smoothed.pdf.txt")
    S06 = load_data("06_recentered.pdf.txt")
    S07 = load_data("07_remove_average.pdf.txt")
    
    S08 = load_data("08_remove_ends.pdf.txt")
    S09 = load_data("09_denoise.pdf.txt")
    S10B = load_data("10B_remove_outlier.pdf.txt")


    retlist = list(S01, S03, S04,
                   S05, S06, S07,
                   S08, S09, S10B)

    names(retlist) = c(
        'S01', 'S03', 'S04',
        'S05', 'S06', 'S07',
        'S08', 'S09', 'S10B')
    
    return(retlist)
}


plot_ref_obs_comparison = function(data_bundle, gene_name) {
    S0x = data_bundle
    obs = na.omit(as.numeric(S0x$obs_matrix[gene_name,]))
    ref = na.omit(as.numeric(S0x$ref_matrix[gene_name,]))

    den_ref = density(ref)
    den_obs = density(obs)

    xrange = range(den_ref$x, den_obs$x)
    yrange = range(den_ref$y, den_obs$y)
    
    plot(density(ref), t='l', xlim=xrange, ylim=yrange);
    points(density(obs), col='blue', t='l');
    
}


main()
