# assumes infercnv is installed.

../scripts/inferCNV.R \
    --delim " " \
    --ref_groups "1:235,236:280" \
    --cutoff 1 \
    --noise_filter 0.25 \
    --output_dir quickstart \
    --vis_bound_threshold " -1,1" \
    --ref example_refs.txt \
      example_expression.txt example_genomic_positions.txt
