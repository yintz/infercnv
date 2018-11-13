

window_size = 11
half_window = (window_size - 1) / 2

normal_types = names(infercnv_obj@reference_grouped_cell_indices)
tumor_groupings = infercnv_obj@observation_grouped_cell_indices

# tumor_groupings = subcluster_tumors(infercnv_obj, tumor_groupings, cut_tree_height_ratio=2.5, hclust_method="ward.D")


infercnv_obj_result = infercnv_obj

for (tumor_type in names(tumor_groupings)) {
    tumor_indices = tumor_groupings[[ tumor_type ]] 
    
    working_data = infercnv_obj@expr.data[, tumor_indices]
    xdim = dim(working_data)[1]
    ydim = dim(working_data)[2]
    # results = matrix(, xdim, ydim)
    results = working_data
    if (xdim >= window_size & ydim >= window_size) {
        for (posx in ((half_window + 1):(xdim - (half_window + 1)))) {
            for ( posy in ((half_window + 1):(ydim - (half_window + 1)))) {
                results[posx, posy] = median(working_data[(posx - half_window):(posx + half_window), (posy - half_window):(posy + half_window)])
            }
        }
    }
    infercnv_obj_result@expr.data[, tumor_indices] = results
    #for (normal_type in normal_types) {
    #    normal_indices = infercnv_obj@reference_grouped_cell_indices[[ normal_type ]]
    #}
}

infercnv_obj_final = infercnv::remove_spike(infercnv_obj_result)

plot_cnv(infercnv_obj=infercnv_obj_final,
         cluster_by_groups=TRUE,
         out_dir="noise_test/",
         color_safe_pal=FALSE,
         x.center=mean(infercnv_obj_result@expr.data),
         x.range="auto",
         title="median_filtering",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         output_filename="median_filtering",
         write_expr_matrix=TRUE
)


#infercnv_obj_final_noreduction = infercnv_obj_final
#infercnv_obj_final_reductionbefore = infercnv_obj_final

plot_cnv(infercnv_obj=infercnv_obj_final_noreduction,
         cluster_by_groups=TRUE,
         out_dir="noise_test/",
         color_safe_pal=FALSE,
         x.center=mean(infercnv_obj_result@expr.data),
         x.range="auto",
         title="median_filtering_noreduc",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         output_filename="median_filtering_noreduc",
         write_expr_matrix=TRUE
)

plot_cnv(infercnv_obj=infercnv_obj_final_reductionbefore,
         cluster_by_groups=TRUE,
         out_dir="noise_test/",
         color_safe_pal=FALSE,
         x.center=mean(infercnv_obj_result@expr.data),
         x.range="auto",
         title="median_filtering_noisereducbefore",
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         output_filename="median_filtering_noisereducbefore",
         write_expr_matrix=TRUE
)


