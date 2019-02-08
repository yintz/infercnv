#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
    
parser = ArgumentParser()
parser$add_argument("--infercnv_obj", help="infercnv_obj file", required=TRUE, nargs=1)
args = parser$parse_args()

library(infercnv)

infercnv_obj_file = args$infercnv_obj

infercnv_obj = readRDS(infercnv_obj_file)

## write counts matrix
write.table(infercnv_obj@count.data, file='sc.counts.matrix', quote=F, sep="\t")

cellnames = colnames(infercnv_obj@count.data)

groupings = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

## write cell annotation file
cell.annots = do.call(rbind, lapply(names(groupings), function(groupname) {
    cell_idx = groupings[[ groupname ]]
    group.cellnames = cellnames[cell_idx]

    return(data.frame(cells=group.cellnames, type=groupname))
}))

cell.annots = cell.annots[ cell.annots$cells %in% colnames(infercnv_obj@count.data), ]

write.table(cell.annots, file="cell_annots.txt", quote=F, row.names=F, col.names=F, sep="\t")

## write infercnv runner:

cat(file='run.infercnv.R', sprintf("#!/usr/bin/env Rscript

options(error = function() { traceback(2); q(status = 1) } )

library(\"infercnv\")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=\"sc.counts.matrix\",
                                    annotations_file=\"cell_annots.txt\",
                                    delim=\"\t\",
                                    gene_order_file=\"gencode_v19_gene_pos.txt\",
                                    ref_group_names=c(\'%s\'))

out_dir=\"output_dir\"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             plot_steps=T,
                             HMM=T,
                             #HMM_mode='subclusters',
                             HMM_mode='samples',
                             sim_method='meanvar'
                             )
", paste(names(infercnv_obj@reference_grouped_cell_indices),collapse="','") ) )

Sys.chmod('run.infercnv.R', mode = "0775")


