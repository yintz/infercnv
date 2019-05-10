sample_object <- function(infercnv_obj,
    n_cells=100,
    every_n=NULL,
    above_m=NULL,
    on_references=TRUE,
    on_observations=TRUE) {

    do_every_n = FALSE
    if (!is.null(every_n) && !is.null(above_m)) {
        if (every_n < 2) {
            flog.error("every_n needs to be at least 2, otherwise nothing will be done.")
            stop("every_n < 2")
        }
        else if (floor(every_n) != every_n) {
            flog.error("every_n needs to be an integer.")
            stop("every_n not an integer")
        }
        else if(!("numeric" %in% is(above_m))) {
            flog.error("above_m needs to be numeric.")
            stop("above_m not numeric")
        }
        do_every_n = TRUE
        n_cells = NULL # mostly to make sure it's not used by inadvertence
    }
    else if (!is.null(every_n) || !is.null(above_m)) {
        flog.error("To use object sampling with every_n and above_m options, please set both. If you want to use random sampling, please do not specify either, only set n_cells.")
        stop("every_n and above_m both need to be set or neither.")
    }
    else if (is.null(n_cells) || n_cells < 1) {
        flog.error("Please provide a valid number of cells to sample to.")
        stop("invalid n_cells")
    }


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

    if (do_every_n) {
        if (on_references == TRUE) {
            col1 = sum(vapply(names(infercnv_obj@reference_grouped_cell_indices), function(sample_name) {
                if (length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) > above_m) {
                    as.integer(ceiling(length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) / every_n))
                }
                else {
                    length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])
                }
            }, integer(1)))
        }
        if (on_observations == TRUE) {
            col2 = sum(vapply(names(infercnv_obj@observation_grouped_cell_indices), function(sample_name) {
                if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) > above_m) {
                    as.integer(ceiling(length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) / every_n))
                }
                else {
                    length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
                }
            }, integer(1)))
        }
        # make the result matrix slightly bigger than expected in case of subclusters not getting sampled 
        # this allows to add a cell from each non represented subcluster without having to extend the matrix
        # we can then remove the empty columns from this extension based on i at the end (should never be more than this extension)
        # which will be more efficient than extending the matrix multiple times
        col2 = col2 + length(unlist(infercnv_obj@tumor_subclusters$subclusters, recursive=FALSE)) - length(infercnv_obj@tumor_subclusters$subclusters)
    }
    else {
        if (on_references == TRUE) {
        	col1 = length(infercnv_obj@reference_grouped_cell_indices) * n_cells
        }
        if (on_observations == TRUE) {
        	col2 = length(infercnv_obj@observation_grouped_cell_indices) * n_cells
        }
    }
    new_obj@expr.data = matrix(nrow=nrow(infercnv_obj@expr.data), ncol=(col1 + col2))
    row.names(new_obj@expr.data) = row.names(infercnv_obj@expr.data)
    futur_colnames = rep("", ncol(new_obj@expr.data))
    new_obj@tumor_subclusters$subclusters=list()

    i = 1
    if (on_references == TRUE) {

        for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {
            if ((!is.null(n_cells) && length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) >= n_cells)
                 || do_every_n) {  # downsample
                if (!do_every_n) {
                    sampled_indices = sample(infercnv_obj@reference_grouped_cell_indices[[sample_name]], size=n_cells, replace=FALSE)
                    # futur_colnames[(i:(i + n_cells - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices])
                }
                else if (length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) > above_m) {  # select 1 every_n because above_m
                    ## sort(unlist(infercnv_obj_subsample@tumor_subclusters$subclusters[[sample_name]], use.names = FALSE)) == infercnv_obj_subsample@reference_grouped_cell_indices[[sample_name]]
                    seq_indices = seq(1, length(unlist(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])), every_n)

                    sampled_labels = (infercnv_obj@tumor_subclusters$hc[[sample_name]]$labels[infercnv_obj@tumor_subclusters$hc[[sample_name]]$order])[seq_indices]
                    sampled_indices = match(sampled_labels, colnames(infercnv_obj@expr.data))

                    # check that all subclusters are represented
                    for (subcluster_id in names(infercnv_obj@tumor_subclusters$subcluster[[sample_name]])) {
                        if (!any(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]] %in% sampled_indices)) {
                            sampled_indices = c(sampled_indices, infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]][1])

                            ##### need to extend by 1 size of expr data matrix
                        }
                    }

                    # sampled_indices = unlist(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])[seq_indices]
                    # futur_colnames[(i:(i + length(seq_indices) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices])
                }
                else { # keep everything because we would select 1 every_n but we are not above_m
                    sampled_indices = infercnv_obj@reference_grouped_cell_indices[[sample_name]]
                    # futur_colnames[(i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices])
                }
                futur_colnames[(i:(i + length(sampled_indices) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE])
            }
            else if (!is.null(n_cells)) {  # upsample
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

    			futur_colnames[(i:(i + length(sampled_indices) - 1))] = paste(colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE]), seq_along(sampled_indices), sep="_")  # paste with seq_along to ensure unique labels
            }
            new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = c(i:(i + length(sampled_indices) - 1))
            new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, sampled_indices, drop=FALSE]
            new_obj@reference_grouped_cell_indices[[sample_name]]=c(i:(i + length(sampled_indices) - 1))
			# futur_colnames[(i:(i + n_cells - 1))] = paste(sample_name, seq_len(n_cells), sep="_")
			i = i + length(sampled_indices)
		}
	}
	else { # do nothing on references
		for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {
            new_obj@expr.data[, (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE]
            new_obj@reference_grouped_cell_indices[[sample_name]] = (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))
            new_obj@tumor_subclusters$hc[[sample_name]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]
            new_obj@tumor_subclusters$subclusters[[sample_name]] = infercnv_obj@tumor_subclusters$subclusters[[sample_name]]
            futur_colnames[(i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE])
            i = i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])
        }
    }

    if (on_observations == TRUE) {
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {
            # if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) >= n_cells) {
            if ((!is.null(n_cells) && length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) >= n_cells)
                 || do_every_n) { ## downsample 
                if (!do_every_n) {  # basic random sampling of n_cells
                    sampled_indices = sample(seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), size=n_cells, replace=FALSE)
                }
                else if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) > above_m) {  # select 1 every_n because above_m
                    ## sort(unlist(infercnv_obj_subsample@tumor_subclusters$subclusters[[sample_name]], use.names = FALSE)) == infercnv_obj_subsample@reference_grouped_cell_indices[[sample_name]]
                    seq_indices = seq(1, length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), every_n)
                    sampled_labels = (infercnv_obj@tumor_subclusters$hc[[sample_name]]$labels[infercnv_obj@tumor_subclusters$hc[[sample_name]]$order])[seq_indices]
                    sampled_indices = match(sampled_labels, colnames(infercnv_obj@expr.data))

                    # check that all subclusters are represented
                    for (subcluster_id in names(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])) {
                        if (!any(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]] %in% sampled_indices)) {
                            sampled_indices = c(sampled_indices, infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]][1])
                        }
                    }
                }
                else { # keep everything because we would select 1 every_n but we are not above_m
                    sampled_indices = infercnv_obj@observation_grouped_cell_indices[[sample_name]]
                }
                ## prune tree
                to_prune = colnames(infercnv_obj@expr.data)[setdiff(infercnv_obj@observation_grouped_cell_indices[[sample_name]], sampled_indices)]
                if (length(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order) > (length(to_prune) + 1) ) {  # more than 2 leaves left, so hclust can exist
                    new_obj@tumor_subclusters$hc[[sample_name]] = prune(infercnv_obj@tumor_subclusters$hc[[sample_name]], to_prune)
                }
                else { # 1 or no leaves left, can't have such an hclust object
                    new_obj@tumor_subclusters$hc[[sample_name]] = NULL
                }
                # new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][sampled_indices], drop=FALSE]
                new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, sampled_indices, drop=FALSE]
                # futur_colnames[(i:(i + length(sampled_indices) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][sampled_indices], drop=FALSE])
                futur_colnames[(i:(i + length(sampled_indices) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE])
			}

	        else if (!is.null(n_cells)) {  ## upsample
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

                new_suffixes = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                        if (x %in% pre_sampled_indices) {
                            seq_len(n_copies + 1)
                        }
                        else {
                            seq_len(n_copies)
                        }
                    }))

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
                # new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order], drop=FALSE]
                new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, new_data_order, drop=FALSE]
                # futur_colnames[(i:(i + n_cells - 1))] = paste(colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order]]), seq_along(new_data_order), sep="_")
                # futur_colnames[(i:(i + length(sampled_indices) - 1))] = paste(colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order]]), new_suffixes, sep="_")
                futur_colnames[(i:(i + length(sampled_indices) - 1))] = paste(colnames(infercnv_obj@expr.data[, new_data_order, drop=FALSE]), new_suffixes, sep="_")

            }
            # else {
                ## should never happen
            # }

            new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = c(i:(i + length(sampled_indices) - 1))
            # new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = seq_along(sampled_indices)  ## don't care about the actual values for plotting purpose, only the length
            new_obj@observation_grouped_cell_indices[[sample_name]]=c(i:(i + length(sampled_indices) - 1))
            i = i + length(sampled_indices)
        }
    }
    else { # do nothing on observations
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {
            new_obj@expr.data[, (i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE]
            new_obj@observation_grouped_cell_indices[[sample_name]] = (i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))
            new_obj@tumor_subclusters$hc[[sample_name]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]
            new_obj@tumor_subclusters$subclusters[[sample_name]] = infercnv_obj@tumor_subclusters$subclusters[[sample_name]]
            futur_colnames[(i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE])
            i = i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
        }
    }
    colnames(new_obj@expr.data) = futur_colnames

    # now drop empty columns/cells that were made in case of every_n sampling missing some subclusters
    if (i <= dim(new_obj@expr.data)[2]) {
        new_obj@expr.data = new_obj@expr.data[, -(i:dim(new_obj@expr.data)[2]), drop=FALSE]
    }

    return(new_obj)
}




############ make method that returns a list of infercnv_obj, 1 for each group of observations/references (to plot on their own)
#### try to add support for k_obs_groups based splitting and not only cluster_by_groups=TRUE

