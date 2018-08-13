# assumes infercnv is installed.

../scripts/inferCNV.R \
    --ref_groups "1:50,51:95" \
    --cutoff 1 \
    --noise_filter 0.2 \
    --output_dir quickstart \
    --vis_bound_threshold " -1,1" \
    --ref normal_cells.csv \
      oligodendroglioma.tp100k.expr.matrix gencode_v19_gene_pos.txt
