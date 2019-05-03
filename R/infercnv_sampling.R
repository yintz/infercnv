sample_object <- function(infercnv_obj,
    n_cells=100,
    on_references=TRUE,
    on_observations=TRUE) {


    new_obj <- new(
        Class = "infercnv",
        expr.data = matrix(), 
        count.data = matrix(),
        gene_order = infercnv_obj@gene_order,
        reference_grouped_cell_indices = list(),
        observation_grouped_cell_indices = list(),
        tumor_subclusters = infercnv_obj@tumor_subclusters,  # mainly copied for structure, will get updated
        .hspike = infercnv_obj@.hspike)

#### TODO add option to sample 1 every n instead of n_cells

    col1 = length(unlist(infercnv_obj@reference_grouped_cell_indices))
    col2 = length(unlist(infercnv_obj@observation_grouped_cell_indices))
    if (on_references == TRUE) {
    	col1 = length(infercnv_obj@reference_grouped_cell_indices) * n_cells
    }
    if (on_observations == TRUE) {
    	col2 = length(infercnv_obj@observation_grouped_cell_indices) * n_cells
    }
    new_obj@expr.data = matrix(nrow=nrow(infercnv_obj@expr.data), ncol=(col1 + col2))
    row.names(new_obj@expr.data) = row.names(infercnv_obj@expr.data)
    futur_colnames = rep("", ncol(new_obj@expr.data))
    new_obj@tumor_subclusters$subclusters=list()

    i = 1
    if (on_references == TRUE) {

        for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {

            if (length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) >= n_cells) {  # downsample
                sampled_indices = sample(infercnv_obj@reference_grouped_cell_indices[[sample_name]], size=n_cells, replace=FALSE)
                futur_colnames[(i:(i + n_cells - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices])
            }
            else {  # upsample
                n_copies = floor(n_cells / length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]))
                to_sample = n_cells %% length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])

                pre_sampled_indices = sample(seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]]), size=to_sample, replace=FALSE)
                sampled_indices = sort(c(pre_sampled_indices, rep(seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]]), n_copies)))

                new_data_order = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                        if (x %in% pre_sampled_indices) {
                            rep(x, (n_copies + 1))
                        }
                        else {
                            rep(x, n_copies)
                        }
                    }))

                ### add check in case there is no hclust for references because it was not set in run() options
                str_newick_to_alter = write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[sample_name]]))
                labels = colnames(infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE])
                for (k in seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]])) {
                    if (k %in% pre_sampled_indices) {
                        # need n_copies + 1 duplicates, so need n_copies added
                        current_label = labels[k]
                        to_replace = current_label
                        for (l in seq_len(n_copies)) {
                            replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                            str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                            to_replace = paste(current_label, "_", (l+1), sep="")
                        }
                    }
                    else if (n_copies > 1) {
                        # make n_copies duplicates, so need n_copies - 1 added
                        current_label = labels[k]
                        to_replace = current_label
                        for (l in seq_len(n_copies - 1)) {
                            replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                            str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                            to_replace = paste(current_label, "_", (l+1), sep="")
                        }
                    }
                }
                new_obj@tumor_subclusters$hc[[sample_name]] = as.hclust(read.tree(text=str_newick_to_alter))

    			futur_colnames[(i:(i + n_cells - 1))] = paste(colnames(infercnv_obj@expr.data[, sampled_indices]), seq_along(sampled_indices), sep="_")  # paste with seq_along to ensure unique labels
            }
            new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, sampled_indices]
            new_obj@reference_grouped_cell_indices[[sample_name]]=c(i:(i + n_cells - 1))
			# futur_colnames[(i:(i + n_cells - 1))] = paste(sample_name, seq_len(n_cells), sep="_")
			i = i + n_cells
		}
	}
	else {
		for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {
            new_obj@expr.data[, (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]]]
            new_obj@reference_grouped_cell_indices[[sample_name]] = (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))
            futur_colnames[(i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]]])
            i = i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])
        }
    }

    if (on_observations == TRUE) {
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {

            if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) >= n_cells) {
				## prune tree

                sampled_indices = sample(seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), size=n_cells, replace=FALSE)
                to_prune = colnames(infercnv_obj@expr.data)[infercnv_obj@observation_grouped_cell_indices[[sample_name]][-sampled_indices]]

                new_obj@tumor_subclusters$hc[[sample_name]] = prune(infercnv_obj@tumor_subclusters$hc[[sample_name]], to_prune)
                new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][sampled_indices], drop=FALSE]
                futur_colnames[(i:(i + n_cells - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][sampled_indices], drop=FALSE])


			}
	        # flog.info(paste(sample_name, " i=", i))
	        # new_obj@expr.data[, c(i:(i+j-1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][1:j]]
	        # a <- list()  # initialize empty object
	        # a$merge <- matrix(c(-1, -2), nc=2, byrow=TRUE)
	        # for (k in c(1:(j-2))) {
	        #     a$merge <- rbind(a$merge, c(k, -(k+2)))
	        # }
	        # a$height <- rep(1, (j-1))    # define merge heights
	        # a$order <- 1:j              # order of leaves(trivial if hand-entered)
	        # a$labels <- paste(sample_name, c(1:j), sep="_")    # labels of leaves
	        # class(a) <- "hclust"        # make it an hclust object
	        
	        # new_obj@tumor_subclusters$hc[[sample_name]] = a
	        # new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = c(i:(i+j-1))
	        
	        # new_obj@observation_grouped_cell_indices[[sample_name]]=c(i:(i+j-1))
	        # futur_colnames[c(i:(i+j-1))] = paste(sample_name, c(1:j), sep="_")
	        # i = i + j

	        else {  ## upsample
                n_copies = floor(n_cells / length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]))
                to_sample = n_cells %% length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])

                #pre_sampled_indices = sample(infercnv_obj@observation_grouped_cell_indices[[sample_name]], size=to_sample, replace=FALSE)
                pre_sampled_indices = sample(seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), size=to_sample, replace=FALSE)
                sampled_indices = sort(c(pre_sampled_indices, rep(seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), n_copies)))  ##

                new_data_order = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                        if (x %in% pre_sampled_indices) {
                            rep(x, (n_copies + 1))
                        }
                        else {
                            rep(x, n_copies)
                        }
                    }))

                # hc_to_alter = infercnv_obj@tumor_subclusters$hc[[sample_name]]
                # fake_order = rep(NA, n_cells)
                # fake_labels = rep(NA, n_cells)
                # fake_merge = matrix(, nrow=(n_cells - 1), ncol=2)
                # fake_height = c(rep(0, ((n_copies - 1) * length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) + to_sample)), hc_to_alter$height)
                # old_merge_to_new_merge = rep(NA, (n_cells - 1))

                str_newick_to_alter = write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[sample_name]]))
                # write.tree(base_newick)
                # plot(read.tree(text=write.tree(base_newick)))
                # gsub(to_replace, replacement, str_newick_to_alter)

                labels = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE])
                for (k in seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]])) {
                    if (k %in% pre_sampled_indices) {
                        # need n_copies + 1 duplicates, so need n_copies added
                        # current_label = row.names(infercnv_obj@expr.data[, k, drop=FALSE])
                        current_label = labels[k]
                        to_replace = current_label
                        for (l in seq_len(n_copies)) {
                            replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                            str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                            to_replace = paste(current_label, "_", (l+1), sep="")
                        }
                    }
                    else if (n_copies > 1) {
                        # make n_copies duplicates, so need n_copies - 1 added
                        # current_label = row.names(infercnv_obj@expr.data[, k, drop=FALSE])
                        current_label = labels[k]
                        to_replace = current_label
                        for (l in seq_len(n_copies - 1)) {
                            replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                            str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                            to_replace = paste(current_label, "_", (l+1), sep="")
                        }
                    }
                    # else {
                    #     # do nothing
                    # }
                }
                new_obj@tumor_subclusters$hc[[sample_name]] = as.hclust(read.tree(text=str_newick_to_alter))
                new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order], drop=FALSE]
                futur_colnames[(i:(i + n_cells - 1))] = paste(colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order]]), seq_along(new_data_order), sep="_")

                # for (k in seq_len(nrow(hc_to_alter$merge))) {
                #     if (hc_to_alter$merge[k, 1] < 0) {
                #         if (hc_to_alter$merge[k, 1] %in% pre_sampled_indices) {
                #             # make n_copies + 1 copies
                #         }
                #         else if (n_copies > 1) {
                #             # make n_copies
                #         }
                #         else {
                #             # no copy needed
                #         }
                #     }
                #     if (hc_to_alter$merge[k, 2] < 0) {
                #     }
                # }
                # infercnv_obj@observation_grouped_cell_indices[[sample_name]][hc_to_alter$order]
            }


            new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = seq_len(n_cells)  ## don't care about the actual values for plotting purpose, only the length
            new_obj@observation_grouped_cell_indices[[sample_name]]=c(i:(i + n_cells - 1))
            i = i + n_cells

            # new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, sampled_indices]
            # new_obj@observation_grouped_cell_indices[[sample_name]]=c(i:(i + n_cells - 1))
            # futur_colnames[(i:(i + n_cells - 1))] = paste(sample_name, seq_len(n_cells), sep="_")


            # for (k in infercnv_obj@observation_grouped_cell_indices[[sample_name]]) {

            #     if (k %in% pre_sampled_indices) {
            #         times = n_copies + 1
            #     }
            #     else {
            #         times = n_copies
            #     }
        }
    }
    colnames(new_obj@expr.data) = futur_colnames

    return(new_obj)
}


