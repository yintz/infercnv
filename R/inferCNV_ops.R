
#' Function doing the actual analysis before calling the plotting functions.
#'
#' @title run() : Invokes a routine inferCNV analysis to Infer CNV changes given a matrix of RNASeq counts.
#'
#' @param infercnv_obj An infercnv object populated with raw count data
#'
#' @param cutoff Cut-off for the min average read counts per gene among reference cells. (default: 1)
#'
#' @param min_cells_per_gene minimum number of reference cells requiring expression measurements to include the corresponding gene.
#'                           default: 3
#'
#' @param out_dir path to directory to deposit outputs (default: '.')
#'
#' ## Smoothing params
#' @param window_length Length of the window for the moving average
#'                          (smoothing). Should be an odd integer. (default: 101)#'
#'
#' @param smooth_method  Method to use for smoothing: c(runmeans,pyramidinal)  default: pyramidinal
#'
#' #####
#' 
#' @param num_ref_groups The number of reference groups or a list of
#'                           indices for each group of reference indices in
#'                           relation to reference_obs. (default: NULL)
#'
#' @param ref_subtract_use_mean_bounds   Determine means separately for each ref group, then remove intensities within bounds of means (default: TRUE)
#'                                       Otherwise, uses mean of the means across groups.
#'
#' #############################
#'
#' @param cluster_by_groups   If observations are defined according to groups (ie. patients), each group
#'                            of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
#'
#'
#' @param k_obs_groups Number of groups in which to break the observations. (default: 1)
#'
#'
#'
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' default("ward.D")
#'
#' @param analysis_mode options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions.
#'                                                          default: samples (fastest, but subclusters is ideal)
#'
#' @param max_centered_threshold The maximum value a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value. (default: 3),
#'                               can set to a numeric value or "auto" to bound by the mean bounds across cells.
#'                               Set to NA to turn off.
#'
#' @param scale_data  perform Z-scaling of logtransformed data (default: FALSE).  This may be turned on if you have
#'                    very different kinds of data for your normal and tumor samples. For example, you need to use GTEx
#'                    representative normal expression profiles rather than being able to leverage normal single cell data
#'                    that goes with your experiment.
#' 
#' #########################################################################
#' ## Downstream Analyses (HMM or non-DE-masking) based on tumor subclusters
#'
#' @param HMM  when set to True, runs HMM to predict CNV level (default: FALSE)
#'
#' @param HMM_transition_prob   transition probability in HMM (default: 1e-6)
#' 
#' @param tumor_subcluster_pval max p-value for defining a significant tumor subcluster (default: 0.01)
#'
#'
#' @param HMM_report_by   cell, consensus, subcluster (default: subcluster)  Note, reporting is performed entirely separately from the HMM prediction.  So, you can predict on subclusters, but get per-cell level reporting (more voluminous output).
#'
#'
#' @param HMM_type  HMM model type. Options: (i6 or i3):
#'                          i6: infercnv 6-state model (0, 0.5, 1, 1.5, 2, >2) where state emissions are calibrated based on simulated CNV levels.
#'                          i3: infercnv 3-state model (del, neutral, amp) configured based on normal cells and HMM_i3_z_pval
#'
#' @param HMM_i3_z_pval     p-value for HMM i3 state overlap (default: 0.05)
#' 
#' 
#' #############################
#' ## de-noising parameters ####
#'
#'
#' @param denoise       If True, turns on denoising according to options below
#'
#' @param noise_filter  Values +- from the reference cell mean will be set to zero (whitening effect)
#'                      default(NA, instead will use sd_amplifier below.
#'
#' @param sd_amplifier  Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier
#'                      default: 1.0
#'
#' @param noise_logistic use the noise_filter or sd_amplifier based threshold (whichever is invoked) as the midpoint in a
#'                       logistic model for downscaling values close to the mean. (default: TRUE)
#'
#'
#' #######################
#' ## Experimental options
#'
#'
#' @param remove_genes_at_chr_ends experimental option: If true, removes the window_length/2 genes at both ends of the chromosome.
#'
#'
#' @param prune_outliers  Define outliers loosely as those that exceed the mean boundaries among all cells.  These are set to the bounds.
#'
#' ## Masking non-DE genes parameters #########
#'
#' @param mask_nonDE_genes If true, sets genes not significantly differentially expressed between tumor/normal to
#'                          the mean value for the complete data set (default: 0.05)
#'
#' @param mask_nonDE_pval  p-value threshold for defining statistically significant DE genes between tumor/normal
#
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.'
#'
#' @param require_DE_all_normals If mask_nonDE_genes is set, those genes will be masked only if they are are found as DE according to test.use and mask_nonDE_pval in each of the comparisons to normal cells options: {"any", "most", "all"} (default: "any")
#'
#' ##################
#' ## Outlier pruning
#'
#' @param outlier_method_bound Method to use for bounding outlier values. (default: "average_bound")
#'                             Will preferentially use outlier_lower_bounda and outlier_upper_bound if set.
#' @param outlier_lower_bound  Outliers below this lower bound will be set to this value.
#' @param outlier_upper_bound  Outliers above this upper bound will be set to this value.
#'
#'
#' #############################################
#'
#' @param plot_steps If true, saves infercnv objects and plots data at the intermediate steps.
#'
#' ##########################
#'
#' @param final_scale_limits The scale limits for the final heatmap output by the run() method. Default "auto". Alt, c(low,high)
#'
#' @param final_center_val   Center value for final heatmap output by the run() method.
#'
#' @param debug If true, output debug level logging.
#'
#' @param num_threads (int) number of threads for parallel steps (default: 4)
#' 
#' @return infercnv_obj containing filtered and transformed data
#'
#' @export
#'

run <- function(infercnv_obj,

                # gene filtering settings
                cutoff=1,
                min_cells_per_gene=3,

                out_dir=".",

                analysis_mode=c('samples', 'subclusters', 'cells'), # for filtering and HMM
                        #
                ## smoothing params
                window_length=101,
                smooth_method=c('pyramidinal', 'runmeans'),
                
                num_ref_groups=NULL,
                ref_subtract_use_mean_bounds=TRUE,
                
                max_centered_threshold=3, # or set to a specific value or "auto", or NA to turn off

                ## tumor subclustering options
                tumor_subcluster_pval=0.05,
                tumor_subcluster_partition_method=c('random_trees', 'qnorm', 'pheight', 'qgamma', 'shc'),


                ## HMM opts
                HMM=FALSE, # turn on to auto-run the HMM prediction of CNV levels
                ## tumor subclustering opts
                HMM_transition_prob=1e-6,
                HMM_report_by=c("subcluster","consensus","cell"),
                HMM_type=c('i6', 'i3'),
                HMM_i3_z_pval=0.05,
                BayesMaxPNormal=0,
                
                
                ## some experimental params
                #sim_method=c('meanvar', 'simple', 'splatter'), ## only meanvar supported, others experimental
                sim_method='meanvar',
                sim_foreground=FALSE, ## experimental

                scale_data=FALSE,
                
                ## noise settings
                denoise=FALSE,
                noise_filter=NA,
                sd_amplifier = 1.5,
                noise_logistic=TRUE, # if false, does complete 'noise' elimination.

                # observation cell clustering settings
                cluster_by_groups=FALSE,
                k_obs_groups=1,

                # outlier adjustment settings
                outlier_method_bound="average_bound",
                outlier_lower_bound=NA,
                outlier_upper_bound=NA,

                hclust_method='ward.D2',
                                
                remove_genes_at_chr_ends=FALSE,

                mask_nonDE_genes=FALSE,
                mask_nonDE_pval=0.05, # use permissive threshold
                test.use='wilcoxon',
                require_DE_all_normals="any",


                plot_steps=FALSE,

                debug=FALSE, #for debug level logging

                prune_outliers=FALSE,

                final_scale_limits = NULL,
                final_center_val = NULL,
                
                reuse_subtracted = TRUE,

                num_threads = 4,
                
                hspike_aggregate_normals = FALSE
                
                ) {


    smooth_method = match.arg(smooth_method)
    HMM_report_by = match.arg(HMM_report_by)
    analysis_mode = match.arg(analysis_mode)
    tumor_subcluster_partition_method = match.arg(tumor_subcluster_partition_method)
    #sim_method = match.arg(sim_method)
    HMM_type = match.arg(HMM_type)
    
    
    if (debug) {
        flog.threshold(DEBUG)
    } else {
        flog.threshold(INFO)
    }
    
    flog.info(paste("::process_data:Start", sep=""))

    infercnv.env$GLOBAL_NUM_THREADS <- num_threads
    if(out_dir != "." & !file.exists(out_dir)){
        dir.create(out_dir)
    }

    step_count = 0;

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))

    infercnv_obj_file=file.path(out_dir, sprintf("%02d_incoming_data.infercnv_obj", step_count))
    if (! (reuse_subtracted && file.exists(infercnv_obj_file)) ) {
        saveRDS(infercnv_obj, infercnv_obj_file)
    }

    ###################################################
    ## Step: removing insufficiently expressed genes
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", step_count))

    # Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    # Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    #  with mean values below cutoff.

    infercnv_obj_file = file.path(out_dir, sprintf("%02d_reduced_by_cutoff.infercnv_obj",step_count))

    if (reuse_subtracted && file.exists(infercnv_obj_file) ) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj = readRDS(infercnv_obj_file)
    } else {

        infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, cutoff)

        ## require each gene to be present in a min number of cells for ref sets

        infercnv_obj <- require_above_min_cells_ref(infercnv_obj, min_cells_per_gene=min_cells_per_gene)

        saveRDS(infercnv_obj, file=infercnv_obj_file)
    }


    ###########################################
    ### STEP: normalization by sequencing depth

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", step_count))

    
    resume_file_token = ifelse( (HMM), paste0("HMM",HMM_type), "")
    
    infercnv_obj_file = file=file.path(out_dir, sprintf("%02d_normalized_by_depth%s.infercnv_obj", step_count, resume_file_token))

    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)

    } else {
        infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj)

        if (HMM && HMM_type == 'i6') {
            ## add in the hidden spike needed by the HMM
            infercnv_obj <- .build_and_add_hspike(infercnv_obj, sim_method=sim_method, aggregate_normals=hspike_aggregate_normals)
            
            if (sim_foreground) {
                infercnv_obj <- .sim_foreground(infercnv_obj, sim_method=sim_method)
            }
        }
        
        saveRDS(infercnv_obj, infercnv_obj_file)
    }


    ###########################
    ## Step: log transformation
    
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", step_count))

    infercnv_obj_file = file.path(out_dir, sprintf("%02d_logtransformed%s.infercnv_obj", step_count, resume_file_token))

    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    } else {

        infercnv_obj <- log2xplus1(infercnv_obj)

        saveRDS(infercnv_obj,
                file=infercnv_obj_file)


        ## Plot incremental steps.
        if (plot_steps){

            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_log_transformed_data",step_count),
                     output_filename=sprintf("infercnv.%02d_log_transformed",step_count),
                     write_expr_matrix=TRUE
                     )
        }
    }


    if (scale_data) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: scaling all expression data\n", step_count))

        infercnv_obj_file=file.path(out_dir, sprintf("%02d_scaled%s.infercnv_obj", step_count, resume_file_token))

        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {
            
            infercnv_obj <- scale_infercnv_expr(infercnv_obj)
        
            saveRDS(infercnv_obj, file=infercnv_obj_file)

            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_scaled",step_count),
                         output_filename=sprintf("infercnv.%02d_scaled",step_count),
                         write_expr_matrix=TRUE)
                
            }
        }
    }
    
    
    ## #################################################
    ## Step: Split the reference data into groups if requested

    if (!is.null(num_ref_groups)) {

        if (! has_reference_cells(infercnv_obj)) {
            stop("Error, no reference cells defined. Cannot split them into groups as requested")
        }
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", step_count, num_ref_groups))

        infercnv_obj_file = file.path(out_dir, sprintf("%02d_split_%02d_refs%s.infercnv_obj", step_count, resume_file_token, num_ref_groups))

        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {
            infercnv_obj <- split_references(infercnv_obj,
                                             num_groups=num_ref_groups,
                                             hclust_method=hclust_method)
            saveRDS(infercnv_obj, file=infercnv_obj_file)
        }

    }
    
    
    
    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method == 'random_trees') {
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))

        resume_file_token = paste0(resume_file_token, ".rand_trees")
        infercnv_obj_file = file.path(out_dir, sprintf("%02d_tumor_subclusters%s.%s.infercnv_obj", step_count, resume_file_token, tumor_subcluster_partition_method))
        
        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {
            infercnv_obj <- infercnv:::define_signif_tumor_subclusters_via_random_smooothed_trees(infercnv_obj,
                                                                                                  p_val=tumor_subcluster_pval,
                                                                                                  hclust_method=hclust_method)
            saveRDS(infercnv_obj, file=infercnv_obj_file)

            if (plot_steps) {

                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         output_filename=sprintf("infercnv.%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         write_expr_matrix=TRUE)

            }
        }
    }




    ## ##################################
    ## Step: Subtract average reference
    ## Since we're in log space, this now becomes log(fold_change)

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (before smoothing)\n", step_count))

    infercnv_obj_file = file.path(out_dir,
                                  sprintf("%02d_remove_ref_avg_from_obs_logFC%s.infercnv_obj", step_count, resume_file_token))

    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    } else {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)

        saveRDS(infercnv_obj, file=infercnv_obj_file)

        if (plot_steps) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("infercnv.%02d_remove_average", step_count),
                     write_expr_matrix=TRUE)
        }
    }


    if (! is.na(max_centered_threshold)) {

        ## #####################################################
        ## Apply maximum centered expression thresholds to data
        ## Cap values between threshold and -threshold, retaining earlier center

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold: %s\n", step_count, max_centered_threshold))

        infercnv_obj_file=file.path(out_dir, sprintf("%02d_apply_max_centered_expr_threshold%s.infercnv_obj", step_count, resume_file_token))

        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {

            threshold = max_centered_threshold
            if (is.character(max_centered_threshold) && max_centered_threshold == "auto") {
                threshold = mean(abs(get_average_bounds(infercnv_obj)))
                flog.info(sprintf("Setting max centered threshoolds via auto to: +- %g", threshold))
            }

            infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, threshold=threshold)

            saveRDS(infercnv_obj, file=infercnv_obj_file)


            ## Plot incremental steps.
            if (plot_steps){

                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_apply_max_centered_expr_threshold",step_count),
                         output_filename=sprintf("infercnv.%02d_apply_max_centred_expr_threshold",step_count),
                         write_expr_matrix=TRUE)

            }
        }
    }







    ###########################################################################
    ## Step: For each cell, smooth the data along chromosome with gene windows

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", step_count))

    infercnv_obj_file = file.path(out_dir, sprintf("%02d_smoothed_by_chr%s.infercnv_obj", step_count, resume_file_token))

    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    } else {
        
        if (smooth_method == 'runmeans') {

            infercnv_obj <- smooth_by_chromosome_runmeans(infercnv_obj, window_length)
        } else if (smooth_method == 'pyramidinal') {

            infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=TRUE)
        } else {
            stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
        }
        
        saveRDS(infercnv_obj, file=infercnv_obj_file)

        ## Plot incremental steps.
        if (plot_steps){

            plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_smoothed_by_chr",step_count),
                     output_filename=sprintf("infercnv.%02d_smoothed_by_chr", step_count),
                     write_expr_matrix=TRUE)
        }
    }


    ##
    ## Step:
    ## Center cells/observations after smoothing. This helps reduce the
    ## effect of complexity.

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", step_count))

    infercnv_obj_file = file.path(out_dir, sprintf("%02d_recentered_cells_by_chr%s.infercnv_obj", step_count, resume_file_token))
    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    } else {
        infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj, method="median")

        saveRDS(infercnv_obj, file=infercnv_obj_file)

        ## Plot incremental steps.
        if (plot_steps) {

            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_centering_of_smoothed",step_count),
                     output_filename=sprintf("infercnv.%02d_centering_of_smoothed", step_count),
                     write_expr_matrix=TRUE)

        }
    }



    ## ##################################
    ## Step: Subtract average reference (adjustment)

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (after smoothing)\n", step_count))

    infercnv_obj_file = file.path(out_dir,
                           sprintf("%02d_remove_ref_avg_from_obs_adjust%s.infercnv_obj", step_count, resume_file_token))

    if (reuse_subtracted && file.exists(infercnv_obj_file)) {
        flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
        infercnv_obj <- readRDS(infercnv_obj_file)
    } else {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)

        saveRDS(infercnv_obj, file=infercnv_obj_file)

        if (plot_steps) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("infercnv.%02d_remove_average", step_count),
                     write_expr_matrix=TRUE)
        }
    }


    ## Step: Remove Ends

    if (remove_genes_at_chr_ends == TRUE) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", step_count))

        infercnv_obj_file = file.path(out_dir, sprintf("%02d_remove_gene_at_chr_ends%s.infercnv_obj", step_count, resume_file_token))

        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {
            infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, window_length)

            saveRDS(infercnv_obj, file=infercnv_obj_file)

            ## Plot incremental steps.
            if (plot_steps){

                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_remove_genes_at_chr_ends",step_count),
                         output_filename=sprintf("infercnv.%02d_remove_genes_at_chr_ends",step_count),
                         write_expr_matrix=TRUE)

            }
        }
    }

    invert_logFC=TRUE
    if (invert_logFC) {
        ## ###########################
        ## Step: invert log transform  (convert from log(FC) to FC)

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", step_count))

        infercnv_obj_file = file.path(out_dir, sprintf("%02d_invert_log_transform%s.infercnv_obj", step_count, resume_file_token))

        if (reuse_subtracted && file.exists(infercnv_obj_file)) {
            flog.info(sprintf("-restoring infercnv_obj from %s", infercnv_obj_file))
            infercnv_obj <- readRDS(infercnv_obj_file)
        } else {

            infercnv_obj <- invert_log2(infercnv_obj)

            saveRDS(infercnv_obj, file=infercnv_obj_file)

            if (plot_steps) {
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_invert_log_transform log(FC)->FC",step_count),
                         output_filename=sprintf("infercnv.%02d_invert_log_FC",step_count),
                         write_expr_matrix=TRUE)

            }
        }
    }

    ## ###################################################################
    ## Done restoring infercnv_obj's from files now under reuse_subtracted
    ## ###################################################################

    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method != 'random_trees') {

        resume_file_token = paste0(resume_file_token, '.', tumor_subcluster_partition_method)
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))

        infercnv_obj_file = file.path(out_dir, sprintf("%02d_tumor_subclusters%s.infercnv_obj", step_count, resume_file_token))

        infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                                                        p_val=tumor_subcluster_pval,
                                                        hclust_method=hclust_method,
                                                        partition_method=tumor_subcluster_partition_method)

        saveRDS(infercnv_obj, file=infercnv_obj_file)

        if (plot_steps) {

            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_tumor_subclusters",step_count),
                     output_filename=sprintf("infercnv.%02d_tumor_subclusters",step_count),
                     write_expr_matrix=TRUE)
        }

    }


    ## This is a milestone step and results should always be examined here.
    infercnv_obj_prelim <- infercnv_obj
    infercnv_obj_file = file.path(out_dir, "preliminary.infercnv_obj")
    saveRDS(infercnv_obj_prelim, file=infercnv_obj_file)
    plot_cnv(infercnv_obj_prelim,
             k_obs_groups=k_obs_groups,
             cluster_by_groups=cluster_by_groups,
             out_dir=out_dir,
             title=sprintf("Preliminary infercnv (pre-noise filtering)",step_count),
             output_filename=sprintf("infercnv.preliminary",step_count),
             write_expr_matrix=TRUE)


    ## Below represent optional downstream analysis steps:

    if (prune_outliers) {

        ##################################
        # STEP: Remove outliers for viz

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", step_count))


        infercnv_obj = remove_outliers_norm(infercnv_obj,
                                            out_method=outlier_method_bound,
                                            lower_bound=outlier_lower_bound,
                                            upper_bound=outlier_upper_bound)

        saveRDS(infercnv_obj,
                file=file.path(out_dir, sprintf("%02d_remove_outlier.infercnv_obj", step_count)))


        ## Plot incremental steps.
        if (plot_steps) {

            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_removed_outliers",step_count),
                     output_filename=sprintf("infercnv.%02d_removed_outliers", step_count),
                     write_expr_matrix=TRUE)
        }
    }


    if (HMM) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: HMM-based CNV prediction\n", step_count))

        if (analysis_mode == 'subclusters') {

            if (HMM_type == 'i6') {
                hmm.infercnv_obj <- predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                             t=HMM_transition_prob)
            } else if (HMM_type == 'i3') {
                hmm.infercnv_obj <- ZHMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                                  z_p_val=HMM_i3_z_pval,
                                                                                  t=HMM_transition_prob)
            } else {
                stop("Error, not recognizing HMM_type")
            }
            
        } else if (analysis_mode == 'cells') {

            if (HMM_type == 'i6') {
                hmm.infercnv_obj <- predict_CNV_via_HMM_on_indiv_cells(infercnv_obj, t=HMM_transition_prob)
            } else if (HMM_type == 'i3') {
                hmm.infercnv_obj <- ZHMM_predict_CNV_via_HMM_on_indiv_cells(infercnv_obj,
                                                                            z_p_val=HMM_i3_z_pval,
                                                                            t=HMM_transition_prob)
            } else {
                stop("Error, not recognizing HMM_type")
            }
            

        } else {
            ## samples mode

            if (HMM_type == 'i6') {
                hmm.infercnv_obj <- predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, t=HMM_transition_prob)
            } else if (HMM_type == 'i3') {
                hmm.infercnv_obj <- ZHMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                                  z_p_val=HMM_i3_z_pval,
                                                                                  t=HMM_transition_prob)
            } else {
                stop("Error, not recognizing HMM_type")
            }
            
        }
        
        ## ##################################
        ## Note, HMM invercnv object is only leveraged here, but stored as file for future use:
        ## ##################################

        hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", analysis_mode)
        
        hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj", step_count, hmm_resume_file_token))
        saveRDS(hmm.infercnv_obj, file=hmm.infercnv_obj_file)
        
        ## report predicted cnv regions:
        generate_cnv_region_reports(hmm.infercnv_obj,
                                    output_filename_prefix=sprintf("%02d_HMM_preds", step_count),
                                    out_dir=out_dir,
                                    by=HMM_report_by)

        
        
        
        if (plot_steps) {
            
            ## Plot HMM pred img
            plot_cnv(infercnv_obj=hmm.infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_HMM_preds",step_count),
                     output_filename=sprintf("infercnv.%02d_HMM_pred",step_count),
                     write_expr_matrix=TRUE,
                     x.center=3,
                     x.range=c(0,6)
                     )
        }
        
        ##############################################################
        # Bayesian Network Mixture Model 
        ##############################################################
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Run Bayesian Network Model on HMM predicted CNV's\n", step_count))
        if (HMM_type == 'i6' & BayesMaxPNormal > 0) {
            hmm.infercnv_obj <- infercnv::inferCNVBayesNet( infercnv_obj    = infercnv_obj_prelim,
                                                            HMM_obj         = hmm.infercnv_obj,
                                                            BayesMaxPNormal = BayesMaxPNormal,
                                                            file_dir        = out_dir,
                                                            postMcmcMethod  = "removeCNV",
                                                            out_dir         = file.path(out_dir, "BayesNetOutput"),
                                                            quietly = TRUE,
                                                            CORES = num_threads)
            mcmc.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.infercnv_obj",
                                                               step_count, hmm_resume_file_token))
            saveRDS(hmm.infercnv_obj, file=mcmc.infercnv_obj_file)
            
            if (plot_steps) {
                ## Plot HMM pred img after cnv removal
                plot_cnv(infercnv_obj=hmm.infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds_Bayes_Net",step_count),
                         output_filename=sprintf("infercnv.%02d_HMM_pred.Bayes_Net",step_count),
                         write_expr_matrix=TRUE,
                         x.center=3,
                         x.range=c(0,6)
                )
            }    
        }
        
        
        ## convert from states to representative  intensity values

        ## 
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Converting HMM-based CNV states to repr expr vals\n", step_count))

        if (HMM_type == 'i6') {
            hmm.infercnv_obj <- assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
        } else if (HMM_type == 'i3') {
            hmm.infercnv_obj <- ZHMM_assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
        }
        
        hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.infercnv_obj",
                                                           step_count, hmm_resume_file_token))
        saveRDS(hmm.infercnv_obj, file=hmm.infercnv_obj_file)
        
        ## Plot HMM pred img
        plot_cnv(infercnv_obj=hmm.infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 title=sprintf("%02d_HMM_preds.repr_intensities",step_count),
                 output_filename=sprintf("infercnv.%02d_HMM_pred.repr_intensities",step_count),
                 write_expr_matrix=TRUE,
                 x.center=1,
                 x.range=c(-1,3)
                 )
    }
    
    
    ## all processes that are alternatives to the HMM prediction wrt DE analysis and/or denoising
    
    ## Step: Filtering significantly DE genes
    if (mask_nonDE_genes) {

        if (!has_reference_cells(infercnv_obj)) {
            stop("Error, cannot mask non-DE genes when there are no normal references set")
        }
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n", step_count))
        
        
        infercnv_obj <- mask_non_DE_genes_basic(infercnv_obj,
                                                p_val_thresh=mask_nonDE_pval,
                                                test.use = test.use,
                                                center_val=mean(infercnv_obj@expr.data),
                                                require_DE_all_normals=require_DE_all_normals)
        
        saveRDS(infercnv_obj,
                file=file.path(out_dir, sprintf("%02d_mask_nonDE%s.infercnv_obj", step_count, resume_file_token)))
        
        ## Plot incremental steps.
        if (plot_steps) {
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     title=sprintf("%02d_mask_nonDE",step_count),
                     output_filename=sprintf("infercnv.%02d_mask_nonDE", step_count),
                     write_expr_matrix=TRUE)
            
        }
    }


    if (denoise) {
        
        ## ##############################
        ## Step: de-noising
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))
        
        if (! is.na(noise_filter)) {
            
            if (noise_filter > 0) {
                flog.info(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
                infercnv_obj <- clear_noise(infercnv_obj,
                                            threshold=noise_filter,
                                            noise_logistic=noise_logistic)
            }
            else {
                ## noise == 0 or negative...
                ## don't remove noise.
            }
            
        }
        else {
            ## default, use quantiles, if NA
            flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", sd_amplifier))
            infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj,
                                                        sd_amplifier = sd_amplifier,
                                                        noise_logistic=noise_logistic)
        }
        
        saveRDS(infercnv_obj,
                file=file.path(out_dir, sprintf("%02d_denoise%s.infercnv_obj", step_count, resume_file_token)))
        
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 title=sprintf("%02d_denoised", step_count),
                 output_filename=sprintf("infercnv.%02d_denoised", step_count),
                 write_expr_matrix=TRUE)
        
        
    
    }
    
    saveRDS(infercnv_obj, file=file.path(out_dir, "run.final.infercnv_obj"))

    if (is.null(final_scale_limits)) {
        final_scale_limits = "auto"
    }
    if (is.null(final_center_val)) {
        final_center_val = 1
    }

    flog.info("\n\n## Making the final infercnv heatmap ##")
    plot_cnv(infercnv_obj,
             k_obs_groups=k_obs_groups,
             cluster_by_groups=cluster_by_groups,
             out_dir=out_dir,
             x.center=final_center_val,
             x.range=final_scale_limits,
             title="inferCNV",
             output_filename="infercnv",
             write_expr_matrix=TRUE)


    return(infercnv_obj)

}

#' Function for Generating a next-generation heatmap
#'
#' @title ngchm() : generates next gen heatmap
#'
#' @param infercnv_obj An infercnv object
#'
#' @param out_dir  output directory (default: '.')
#'
#' @param title title of the interactive heatmap (default: "NGCHM")
#'
#' @param gene_symbol ##TODO  (default: NULL)
#'
#' @param path_to_shaidyMapGen path to the shaidyMapGen jar file (default: NULL)
#'
#' @param x.center (integer) Center expression value for heatmap coloring.
#'
#' @param x.range (integer) Values for minimum and maximum thresholds for heatmap coloring.
#'
#' @export
#'

ngchm <- function(infercnv_obj,
                       out_dir=".",
                       title="NGCHM",
                       gene_symbol=NULL,
                       path_to_shaidyMapGen=NULL,
                       x.range = NA,
                       x.center = NA) {

    if (!is.null(path_to_shaidyMapGen)) {
        shaidy.path <- unlist(strsplit(path_to_shaidyMapGen, split = .Platform$file.sep))
        if (!file.exists(path_to_shaidyMapGen) || tail(shaidy.path, n = 1L) != "ShaidyMapGen.jar"){
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using the parameter \"path_to_shaidyMapGen\".",
                                   "Check that the correct pathway is being used.")
            flog.error(error_message)
            stop(error_message)
        }
    } else {
        path_to_shaidyMapGen <- Sys.getenv("SHAIDYMAPGEN")
        if (!file.exists(path_to_shaidyMapGen)){ ## check if envionrmental variable is passed
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using SHAIDYMAPGEN.",
                                   "Check that the correct pathway is being used.")
            flog.error(error_message)
            stop(error_message)
        }
    }

    if (!requireNamespace("NGCHM", quietly=TRUE)) {
        stop("The \"NGCHM\" library is required to use \"-ngchm=TRUE\" but it is not available.", .call=FALSE)
    }

    flog.info("Creating NGCHM as infercnv.ngchm")
    Create_NGCHM(infercnv_obj = infercnv_obj,
                 path_to_shaidyMapGen = path_to_shaidyMapGen,
                 out_dir = out_dir,
                 title = title,
                 gene_symbol = gene_symbol,
                 x.range = x.range,
                 x.center = x.center)
}



#' Subtracting the mean of the reference expr distributions from the observed cells.
#'
#' @title subtract_ref_expr_from_obs()
#'
#' @description Remove the average of the genes of the reference observations from all
#' observations' expression. In the case there are multiple reference groupings,
#' the averages are computed separately for each reference grouping, and the min|max
#' of the averages are subtracted from the observation expression levels. Any values within the range
#' of the min,max of the group are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param inv_log mean values will be determined based on (2^x -1)
#'
#' @return infercnv_obj containing the reference subtracted values.
#'
#' @export
#'

subtract_ref_expr_from_obs <- function(infercnv_obj, inv_log=FALSE, use_bounds=FALSE) {
                                            # r = genes, c = cells
    flog.info(sprintf("::subtract_ref_expr_from_obs:Start inv_log=%s, use_bounds=%s", inv_log, use_bounds))
    

    if (has_reference_cells(infercnv_obj)) {
        ref_groups = infercnv_obj@reference_grouped_cell_indices
        flog.info("subtracting mean(normal) per gene per cell across all data")
    } else {
        ref_groups = list('proxyNormal' = unlist(infercnv_obj@observation_grouped_cell_indices))
        flog.info("-no reference cells specified... using mean of all cells as proxy")
    }

    ref_grp_gene_means <- .get_normal_gene_mean_bounds(infercnv_obj@expr.data, ref_groups, inv_log=inv_log)
    
    infercnv_obj@expr.data <- .subtract_expr(infercnv_obj@expr.data, ref_grp_gene_means, inv_log=inv_log, use_bounds=use_bounds)
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- subtract_ref_expr_from_obs(infercnv_obj@.hspike, inv_log=inv_log, use_bounds=use_bounds)
    }
    
    return(infercnv_obj)
    
}



.get_normal_gene_mean_bounds <- function(expr.data, ref_groups, inv_log=FALSE) {

    get_indiv_gene_group_means_bounds_fun <- function(x) {

        grp_means = c()

        for (ref_group in ref_groups) {

            if (inv_log) {
                ref_grp_mean = log2(mean(2^x[ref_group] - 1) + 1)
            } else {
                ref_grp_mean = mean(x[ref_group])
            }

            grp_means = c(grp_means, ref_grp_mean)

        }
        
        names(grp_means) <- names(ref_groups)
        
        return(as.data.frame(t(data.frame(grp_means))))
    }
    
    gene_ref_grp_means <- do.call(rbind, apply(expr.data, 1, get_indiv_gene_group_means_bounds_fun))

    rownames(gene_ref_grp_means) <- rownames(expr.data)
    
    return(gene_ref_grp_means)
}



.subtract_expr <- function(expr_matrix, ref_grp_gene_means, inv_log=FALSE, use_bounds=FALSE) {

    my.rownames = rownames(expr_matrix)
    my.colnames = colnames(expr_matrix)

    flog.info(sprintf("-subtracting expr per gene, use_bounds=%s", use_bounds))
    
    subtract_normal_expr_fun <- function(row_idx) {

        gene_means <- as.numeric(ref_grp_gene_means[row_idx, , drop=T])

        gene_means_mean <- mean(gene_means)
        
        x <- as.numeric(expr_matrix[row_idx, , drop=T])
        
        row_init = rep(0, length(x))
        
        if (use_bounds) {
            
            grp_min = min(gene_means)
            grp_max = max(gene_means)

            above_max = which(x>grp_max)
            below_min = which(x<grp_min)
                    
            row_init[above_max] <- x[above_max] - grp_max
            row_init[below_min] <- x[below_min] - grp_min
            ## note, in-between values are left at zero!
            
        } else {
            
            row_init <- x - gene_means_mean

        }
                
        return(row_init)
    }
    
    subtr_data <- do.call(rbind, lapply(1:nrow(expr_matrix), subtract_normal_expr_fun))
    rownames(subtr_data) <- my.rownames
    colnames(subtr_data) <- my.colnames
    
    return(subtr_data)

}


# Not testing, parameters ok.
# Helper function allowing greater control over the steps in a color palette.
# Source:http://menugget.blogspot.com/2011/11/define-color-steps-for-
#               colorramppalette.html#more

# Args:
# steps Vector of colors to change use in the palette
# between: Steps where gradients change
#
# Returns:
# Color palette
color.palette <- function(steps,
                          between=NULL, ...){

    if (is.null(between)){
        between <- rep(0, (length(steps) - 1))
    }
    if (length(between) != length(steps) - 1){
        stop("Must have one less \"between\" value than steps")
    }

    fill.steps <- cumsum(rep(1, length(steps)) + c(0, between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[, fill.steps] <- col2rgb(steps)
    for(i in which(between > 0)){
        col.start <- RGB[, fill.steps[i]]
        col.end <- RGB[, fill.steps[i + 1]]
        for(j in seq(3)){
            vals <- seq(col.start[j],
                        col.end[j],
                        length.out=between[i] + 2)
            vals <- vals[2:(2 + between[i] - 1)]
            RGB[j, (fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue=255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Create a sepList forthe heatmap.3 plotting function given integer vectors
# of rows and columns where speration should take place.
# The expected input to the heatmap function is a list of 2 lists.
# The first list are column based rectangles, and the second row.
# To define a rectagle the index of the row or column where the line of the rectagle
# should be placed is done with a vector of integers, left, bottom, right and top line.
# Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
#
# Args:
# row_count Total number of rows
# col_count Total number of columns
# row_seps Vector of integers indices for row breaks
# col_seps Vector of integer indices for column breaks
#
# Returns
# List of lists of vectors
create_sep_list <- function(row_count,
                            col_count,
                            row_seps=NULL,
                            col_seps=NULL){
    sepList <- list()
    # Heatmap.3 wants a list of boxes for seperating columns
    # Column data
    if(!is.null(col_seps) &&
       !is.na(col_seps) &&
       (length(col_seps)>0) &&
       col_count > 0){
        colList <- list()
        for(sep in 1:length(col_seps)){
            colList[[sep]] <- c(col_seps[sep],0,col_seps[sep],row_count)
        }
        sepList[[1]] <- colList
    } else {
        sepList[[1]] <- list()
        sepList[[1]][[1]] <- c(0,0,0,0)
    }

    # Row data
    # This is measured from bottom to top
    # So you have to adjust the values of the data
    row_seps <- row_count-row_seps
    if(!is.null(row_seps) &&
       !is.na(row_seps) &&
       (length(row_seps)>0) &&
       row_count > 0){
        rowList <- list()
        for(sep in 1:length(row_seps)){
            rowList[[sep]] <- c(0,row_seps[sep],col_count,row_seps[sep])
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}


#' @title split_references()
#'
#' @description Split up reference observations in to k groups based on hierarchical clustering.
#'
#'
#' @param infercnv_obj infercnv_object
#'
#' @param num_groups (default: 2)
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return infercnv_obj
#'
#' @export
#'

split_references <- function(infercnv_obj,
                             num_groups=2,
                             hclust_method='complete') {


    flog.info(paste("::split_references:Start", sep=""))

    ref_expr_matrix = infercnv_obj@expr.data[ , get_reference_grouped_cell_indices(infercnv_obj) ]

    hc <- hclust(dist(t(ref_expr_matrix)), method=hclust_method)

    split_groups <- cutree(hc, k=num_groups)

    ref_groups <- list()

    grp_counter = 0
    for (cut_group in unique(split_groups)) {
        grp_counter = grp_counter + 1
        grp_name = sprintf("refgrp-%d", grp_counter)
        cell_names <- names(split_groups[split_groups == cut_group])
        ref_groups[[grp_name]] <- which(colnames(ref_expr_matrix) %in% cell_names)
    }

    infercnv_obj@reference_grouped_cell_indices <- ref_groups

    return(infercnv_obj)
}


#' @title remove_outliers_norm()
#'
#' @description Set outliers to some upper or lower bound.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param out_method method for computing the outlier bounds (default: "average_bound", involving
#'                   determining the range of values for each cell, and then taking the mean of those bounds.)
#'
#' @param lower_bound setting the lower bound for the data (default: NA, uses out_method above)
#'
#' @param upper_bound setting the upper bound for the data (default: NA, uses out_method above)
#'
#' @return infercnv_obj with data bounds set accordingly.
#'
#' @export
#'

remove_outliers_norm <- function(infercnv_obj,
                                 out_method="average_bound",
                                 lower_bound=NA,
                                 upper_bound=NA) {

    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))


    infercnv_obj@expr.data <- .remove_outliers_norm(data=infercnv_obj@expr.data,
                                                    out_method=out_method,
                                                    lower_bound=lower_bound,
                                                    upper_bound=upper_bound)

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- remove_outliers_norm(infercnv_obj@.hspike, out_method, lower_bound, upper_bound)
    }

    return(infercnv_obj)

}


.remove_outliers_norm <- function(data,
                                  out_method="average_bound",
                                  lower_bound=NA,
                                  upper_bound=NA) {

    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))



    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        flog.error("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }



    if ( (! is.na(lower_bound)) & (! is.na(upper_bound)) ) {

        ## using the values as specified.
        flog.info(paste("::remove_outlier_norm: using hard thresholds: ",
                        "lower_bound:" , lower_bound,
                        "upper_bound:", upper_bound) )

    } else if (! is.na(out_method)) {

        # using out_method instead of specified bounds.
        flog.info(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))

        if (out_method == "average_bound"){

            bounds = .get_average_bounds(data)
            lower_bound = bounds[1]
            upper_bound = bounds[2]

            flog.info(sprintf("outlier bounds defined between: %g - %g", lower_bound, upper_bound))

        } else {
            flog.error(paste("::remove_outlier_norm:Error, please",
                                    "provide an approved method for outlier",
                                    "removal for visualization."))
            stop(991)
        }
    } else {
        flog.error("::remove_outlier_norm:Error, must specify outmethod or define exact bounds")
        stop(992)
    }

    # apply bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound


    return(data)
}




#' @title center_cell_expr_across_chromosome()
#'
#' @description Centers expression data across all genes for each cell, using the cell mean or median expression
#' value as the center.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param method method to select the center of the cell expression value. (default: 'mean', options: 'mean,'median')
#'
#' @return infercnv_object
#'
#' @export
#'

center_cell_expr_across_chromosome <- function(infercnv_obj, method="mean") { # or median

    flog.info(paste("::center_smooth across chromosomes per cell"))

    infercnv_obj@expr.data <- .center_columns(infercnv_obj@expr.data, method)


    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- center_cell_expr_across_chromosome(infercnv_obj@.hspike, method)
    }


    return(infercnv_obj)
}

.center_columns <- function(expr_data, method) {

    # Center within columns (cells)
    if (method == "median") {
        row_median <- apply(expr_data, 2, function(x) { median(x, na.rm=T) } )

        expr_data <- t(apply(expr_data, 1, "-", row_median))
    }
    else {
        # by mean
        row_means <- apply(expr_data, 2, function(x) { mean(x, na.rm=T) } )

        expr_data <- t(apply(expr_data, 1, "-", row_means))
    }
    return(expr_data)
}




#' @title require_above_min_mean_expr_cutoff ()
#'
#' @description Filters out genes that have fewer than the corresponding mean value across all cell values.
#'
#' @param infercnv_obj  infercnv_object
#'
#' @param min_mean_expr_cutoff  the minimum mean value allowed for a gene to be retained in the expression matrix.
#'
#' @return infercnv_obj  the infercnv_object with lowly or unexpressed genes removed.
#'
#' @export
#'

require_above_min_mean_expr_cutoff <- function(infercnv_obj, min_mean_expr_cutoff) {

    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))


    indices <-.below_min_mean_expr_cutoff(infercnv_obj@expr.data, min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(sprintf("Removing %d genes from matrix as below mean expr threshold: %g",
                          length(indices), min_mean_expr_cutoff))

        infercnv_obj <- remove_genes(infercnv_obj, indices)

        expr_dim = dim(infercnv_obj@expr.data)
        flog.info(sprintf("There are %d genes and %d cells remaining in the expr matrix.",
                          expr_dim[1], expr_dim[2]))

    }

    return(infercnv_obj)

}


.below_min_mean_expr_cutoff <- function(expr_data, min_mean_expr) {

    average_gene <- rowMeans(expr_data)

    # Find averages above a certain threshold
    indices <- which(average_gene < min_mean_expr)

    return(indices)

}




#' @title require_above_min_cells_ref()
#'
#' @description Filters out genes that have fewer than specified number of cells expressing them.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#'
#' @return infercnv_obj infercnv_object with corresponding genes removed.
#'
#' @export
#'

require_above_min_cells_ref <- function(infercnv_obj, min_cells_per_gene) {

    genes_passed = which(apply(infercnv_obj@expr.data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))

    num_genes_total = dim(infercnv_obj@expr.data)[1]
    num_removed = num_genes_total - length(genes_passed)
    if (num_removed > 0) {

        flog.info(sprintf("Removed %d genes having fewer than %d min cells per gene = %g %% genes removed here",
                          num_removed, min_cells_per_gene, num_removed / num_genes_total * 100))


        if (num_removed == num_genes_total) {

            flog.warn(paste("::All genes removed! Must revisit your data..., cannot continue here."))
            stop(998)
        }


        infercnv_obj <- remove_genes(infercnv_obj, -1 * genes_passed)


    }
    else {

        flog.info("no genes removed due to min cells/gene filter")

    }

    return(infercnv_obj)

}


#' @title clear_noise()
#'
#' @description Remove values that are too close to the reference cell expr average and are considered noise.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param threshold values within reference mean +- threshold are set to zero.
#'
#' @return infercnv_obj
#'
#' @export
#'

clear_noise <- function(infercnv_obj, threshold, noise_logistic=FALSE) {

    flog.info(paste("********* ::clear_noise:Start. threshold: ", threshold, sep=""))

    if (threshold == 0) {
        return(infercnv_obj); # nothing to do
    }

    if (has_reference_cells(infercnv_obj)) {
        ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
        mean_ref_vals = mean(infercnv_obj@expr.data[,ref_idx])
    } else {
        ## no reference
        ## use mean of all data
        mean_ref_vals = mean(infercnv_obj@expr.data)
    }
    
    if (noise_logistic) {

        infercnv_obj <- depress_log_signal_midpt_val(infercnv_obj, mean_ref_vals, threshold)

    } else {

        infercnv_obj@expr.data <- .clear_noise(infercnv_obj@expr.data, threshold, center_pos=mean_ref_vals)
    }

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- clear_noise(infercnv_obj@.hspike, threshold, noise_logistic)
    }

    return(infercnv_obj)
}


.clear_noise <- function(expr_data, threshold, center_pos=0) {

    upper_bound = center_pos + threshold
    lower_bound = center_pos - threshold

    expr_data[expr_data > lower_bound & expr_data < upper_bound] = center_pos

    return(expr_data)
}



#' @title clear_noise_via_ref_mean_sd()
#'
#' @description Define noise based on the standard deviation of the reference cell expression data.
#' The range to remove noise would be mean +- sdev * sd_amplifier
#' where sd_amplifier expands the range around the mean to be removed as noise.
#' Data points defined as noise are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise
#'                     range (default: 1.5)
#'
#' @export
#'

clear_noise_via_ref_mean_sd <- function(infercnv_obj, sd_amplifier=1.5, noise_logistic=FALSE) {

    if (has_reference_cells(infercnv_obj)) {
        ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
        flog.info("denoising using mean(normal) +- sd_amplifier * sd(normal) per gene per cell across all data")
    }
    else {
        ref_idx = unlist(infercnv_obj@observation_grouped_cell_indices)
        flog.info("-no reference cells specified... using mean and sd of all cells as proxy for denoising")
    }
    vals = infercnv_obj@expr.data[,ref_idx]

    mean_ref_vals = mean(vals)

    mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=T))) * sd_amplifier

    upper_bound = mean_ref_vals + mean_ref_sd
    lower_bound = mean_ref_vals - mean_ref_sd

    flog.info(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                           lower_bound, "-", upper_bound, sep=" "))



    if (noise_logistic) {

        threshold = mean_ref_sd
        infercnv_obj <- depress_log_signal_midpt_val(infercnv_obj, mean_ref_vals, threshold)

    } else {
        smooth_matrix <- infercnv_obj@expr.data

        smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = mean_ref_vals

        infercnv_obj@expr.data <- smooth_matrix
    }

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- clear_noise_via_ref_mean_sd(infercnv_obj@.hspike, sd_amplifier, noise_logistic)
    }

    return(infercnv_obj)
}




# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr Indices of the chr in which the tails are to be removed.
# tail_length Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.
.remove_tails <- function(smooth_matrix, chr, tail_length){

    #flog.info(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
         tail_length <- floor(chr_length / 3)
    }
    remove_indices <- chr[1:tail_length]
    remove_indices <- c(remove_indices,
                        chr[ ( (chr_length + 1) - tail_length):
                             chr_length])

    return(remove_indices)
}


#' @title smooth_by_chromosome()
#'
#' @description Smooth expression values for each cell across each chromosome by using a
#' moving average with a window of specified length.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of window (number of genes) for the moving average
#'
#' @param smooth_ends perform smoothing at the ends of the chromosomes (default: TRUE)
#'
#' @return infercnv_obj
#'
#' @export
#'

smooth_by_chromosome <- function(infercnv_obj, window_length, smooth_ends=TRUE) {

    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))

    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))

        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=F]

        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window(data=chr_data,
                                               window_length=window_length)

            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))

            infercnv_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- smooth_by_chromosome(infercnv_obj@.hspike, window_length, smooth_ends)
    }


    return(infercnv_obj)
}


.smooth_window <- function(data, window_length) {

    flog.debug(paste("::smooth_window:Start.", sep=""))

    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }

    ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
    
    data_sm <- apply(data,
                     2,
                     .smooth_helper,
                     window_length=window_length)
    
    
    ## Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)


    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))

    
    return(data_sm)
}


# Helper function for smoothing with a moving average.
#
# Args:
# obs_data: Data to smooth
# window_length:  Length of the window to smooth.
#
# Returns:
# Data smoothed.
##.smooth_helper <- function(obs_data, tail_length) {

.smooth_helper <- function(obs_data, window_length) {
    # strip NAs out and replace after smoothing
    orig_obs_data = obs_data

    nas = is.na(obs_data)

    obs_data = obs_data[!nas]

    obs_length <- length(obs_data)
    end_data <- obs_data

    tail_length = (window_length - 1)/2
    if (obs_length >= window_length) {
        end_data <- .smooth_center_helper(obs_data, window_length)
    }

    # end_data will have the end positions replaced with mean values, smoothing just at the ends.

    obs_count <- length(obs_data)

    numerator_counts_vector = c(c(1:tail_length), tail_length + 1, c(tail_length:1))

    # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
    iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))

    for (tail_end in 1:iteration_range) {
        end_tail = obs_count - tail_end + 1

        d_left = tail_end - 1
        d_right = obs_count - tail_end
        d_right = ifelse(d_right > tail_length, tail_length, d_right)

        r_left = tail_length - d_left
        r_right = tail_length - d_right

        denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)

        left_input_vector_chunk = obs_data[1:(tail_end + d_right)]
        right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]

        numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]

        end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
        end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
    }

    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data

    return(orig_obs_data)
}

# Smooth vector of values over the given window length.
#
# Args:
# obs_data Vector of data to smooth with a moving average.
# window_length Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.
.smooth_center_helper <- function(obs_data, window_length){

    nas = is.na(obs_data)
    vals = obs_data[! nas]

    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(c(1:((window_length-1)/2)), ((window_length-1)/2)+1, c(((window_length-1)/2):1))

    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)

#    flog.info(paste("custom filter = ", custom_filter, "\n and window_length =", window_length, "\n and nrow(data) = ", nrow(obs_data), sep=""))

    #smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)
    smoothed = filter(vals, custom_filter, sides=2)

    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]

    obs_data[! nas] = vals

    return(obs_data)
}


smooth_by_chromosome_runmeans <- function(infercnv_obj, window_length) {

    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))

    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))

        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=F]

        if (nrow(chr_data) > 1) {
            chr_data = apply(chr_data, 2, caTools::runmean, k=window_length)

            infercnv_obj@expr.data[chr_genes_indices, ] <- chr_data
        }
    }

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- smooth_by_chromosome_runmeans(infercnv_obj@.hspike, window_length)
    }


    return(infercnv_obj)
}






#' @title get_average_bounds()
#'
#' @description Computes the mean of the upper and lower bound for the data across all cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return (lower_bound, upper_bound)
#'
#' @export
#'

get_average_bounds <- function (infercnv_obj) {

    return(.get_average_bounds(infercnv_obj@expr.data))

}

.get_average_bounds <- function(expr_matrix) {

    lower_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))

    return(c(lower_bound, upper_bound))
}

#' @title log2xplus1()
#'
#' @description Computes log(x+1), updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

log2xplus1 <- function(infercnv_obj) {

    flog.info("transforming log2xplus1()")

    infercnv_obj@expr.data <- log2(infercnv_obj@expr.data + 1)

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- log2xplus1(infercnv_obj@.hspike)
    }

    return(infercnv_obj)

}



#' @title invert_log2xplus1()
#'
#' @description Computes 2^x - 1
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

invert_log2xplus1 <- function(infercnv_obj) {

    flog.info("inverting log2xplus1()")

    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data - 1

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- invert_log2xplus1(infercnv_obj@.hspike)
    }

    return(infercnv_obj)
}


#' @title invert_log2()
#'
#' @description Computes 2^x
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

invert_log2 <- function(infercnv_obj) {

    flog.info("invert_log2(), computing 2^x")

    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- invert_log2(infercnv_obj@.hspike)
    }

    return(infercnv_obj)
}


#' @title make_zero_NA()
#'
#' @description Converts zero to NA
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

make_zero_NA <- function(infercnv_obj) {

    flog.info("make_zero_NA()")

    infercnv_obj@expr.data <- infercnv_obj@expr.data[infercnv_obj@expr.data == 0] <- NA


    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- make_zero_NA(infercnv_obj@.hspike)
    }

    return(infercnv_obj)

}


#' @title transform_to_reference_based_Zscores()
#'
#' @description Computes mean and standard deviation for the reference cells, then uses these values
#' to compute gene Z-scores for both the reference and observation cells.
#'
#' Note, reference cell gene expression values will be centered at zero after this operation.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

transform_to_reference_based_Zscores <- function(infercnv_obj) {

    # center and convert to z-scores
    flog.info(paste("::center_and_Zscore_conversion", sep=""))

    # remember, genes are rows, cells are cols

    # centering and z-scores based on the reference (normal) cells:

    # ref data represent the null distribution
    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)

    ref_data = infercnv_obj@expr.data[,ref_idx]

    gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=T)})
    gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=T)})

    # assume at least Poisson level variation
    gene_ref_sd = pmax(gene_ref_sd, gene_ref_mean)

    # center all genes at the ref (normal) center:
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_mean, FUN="-")

    # convert to z-scores based on the ref (normal) distribution
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- transform_to_reference_based_Zscores(infercnv_obj@.hspike)
    }

    return(infercnv_obj)

}


#' @title mean_center_gene_expr()
#'
#' @description mean-center all gene expression values across all cells
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @export
#'

mean_center_gene_expr <- function(infercnv_obj) {

    flog.info(paste("::centering", sep=""))

    infercnv_obj@expr.data <- sweep(infercnv_obj@expr.data, 1, rowMeans(infercnv_obj@expr.data, na.rm=T), FUN="-")

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- mean_center_gene_expr(infercnv_obj@.hspike)
    }

    return(infercnv_obj)
}


#' @title get_reference_grouped_cell_indices()
#'
#' @description Retrieves the matrix indices for the columns correspoinding to the reference cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return vector of column indices
#'
#' @export
#'

get_reference_grouped_cell_indices <- function(infercnv_obj) {

    return( unlist(infercnv_obj@reference_grouped_cell_indices) )

}


#' @title apply_max_threshold_bounds()
#'
#' @description Assumes centered at zero and sets bounds to +- threshold value.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param threshold value to threshold the data
#'
#' @export
#'

apply_max_threshold_bounds <- function(infercnv_obj, threshold) {

    flog.info(paste("::process_data:setting max centered expr, threshold set to: +/-: ", threshold))

    infercnv_obj@expr.data[infercnv_obj@expr.data > threshold] <- threshold
    infercnv_obj@expr.data[infercnv_obj@expr.data < (-1 * threshold)] <- -1 * threshold

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- apply_max_threshold_bounds(infercnv_obj@.hspike, threshold)
    }

    return(infercnv_obj)
}


#' @title remove_genes_at_ends_of_chromosomes()
#'
#' @description Removes genes that are within window_length/2 of the ends of each chromosome.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of the window to use.
#'
#' @export
#'

remove_genes_at_ends_of_chromosomes <- function(infercnv_obj, window_length) {

    contig_tail = (window_length - 1) / 2

    remove_indices <- c()
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    for (chr in chrs){
        #flog.info(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- .remove_tails(infercnv_obj@expr.data,
                                    which(gene_chr_listing == chr),
                                    contig_tail)

        #flog.debug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))

        remove_indices <- c(remove_indices, remove_chr)

    }
    if (length(remove_indices) > 0){

        infercnv_obj = remove_genes(infercnv_obj, remove_indices)

        flog.info(paste("::process_data:Remove genes at chr ends, ",
                        "new dimensions (r,c) = ",
                        paste(dim(infercnv_obj@expr.data), collapse=","),
                        " Total=", sum(infercnv_obj@expr.data, na.rm=TRUE),
                        " Min=", min(infercnv_obj@expr.data, na.rm=TRUE),
                        " Max=", max(infercnv_obj@expr.data, na.rm=TRUE),
                        ".", sep=""))

    }
    else {
        flog.error("No genes removed at chr ends.... something wrong here")
        stop(1234)
    }


    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- remove_genes_at_ends_of_chromosomes(infercnv_obj@.hspike, window_length)
    }

    return(infercnv_obj)

}


#' @title normalize_counts_by_seq_depth()
#'
#' @description Normalizes count data by total sum scaling
#'
#' For single cell data, a typical normalization factor is 1e5, providing counts per 100k total counts.
#' If a normalization factor is not provided, the median lib size is used.:
#'
#' @param infercnv_obj infercnv_object
#'
#' @param normalize_factor  total counts to scale the normalization to (default: NA, computed as described above)
#'
#' @export
#'

normalize_counts_by_seq_depth <- function(infercnv_obj, normalize_factor=NA) {

    data <- infercnv_obj@expr.data

    normalized_data <- .normalize_data_matrix_by_seq_depth(data, normalize_factor)

    
    infercnv_obj@expr.data <- normalized_data

    return(infercnv_obj)

}

.normalize_data_matrix_by_seq_depth <- function(counts.matrix, normalize_factor=NA) {

    flog.info("normalizing counts matrix by depth")
    
    data <- counts.matrix
    
    cs = colSums(data)

    print(cs)

    # make fraction of total counts:
    data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")

    if (is.na(normalize_factor)) {

        normalize_factor = median(cs)

        flog.info(sprintf("Computed total sum normalization factor as median libsize: %f", normalize_factor))

    } else {
        flog.info(sprintf("Using specified normalization factor: %f", normalize_factor))
    }

    if (is.na(normalize_factor)) {
        stop("Error, normalize factor not estimated")
    }

    data <- data * normalize_factor

    return(data)
    
}



#' @title anscombe_transform()
#'
#' @description Performs Anscombe's transformation:
#'    y = 2 * sqrt(x + 3/8)
#' as per
#' https://en.wikipedia.org/wiki/Anscombe_transform
#'
#' @param infercnv_obj infercnv_object
#'
#' @export
#'

anscombe_transform <- function(infercnv_obj) {

    infercnv_obj@expr.data <- 2 * sqrt(infercnv_obj@expr.data + 3/8)

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- anscombe_transform(infercnv_obj@.hspike)
    }

    return(infercnv_obj)

}

#' @keywords internal
#' @noRd
#'
add_pseudocount <- function(infercnv_obj, pseudocount) {

    flog.info(sprintf("Adding pseudocount: %g", pseudocount))

    infercnv_obj@expr.data = infercnv_obj@expr.data + pseudocount

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- add_pseudocount(infercnv_obj@.hspike, pseudocount)
    }

    return(infercnv_obj)
}


scale_infercnv_expr <- function(infercnv_obj) {

    flog.info("-scaling expr data")
    infercnv_obj@expr.data = t(scale(t(infercnv_obj@expr.data)))

    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- scale_infercnv_expr(infercnv_obj@.hspike)
    }

    return(infercnv_obj)
}


cross_cell_normalize <- function(infercnv_obj) {

    ## using upper quartile normalization

    flog.info("-cross cell normalization")

    upper_quart = apply(infercnv_obj@expr.data, 2, quantile, probs=0.75)
    mean_upper_quart = mean(upper_quart)
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 2, mean_upper_quart/upper_quart, "*")


    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- cross_cell_normalize(infercnv_obj@.hspike)
    }

    return(infercnv_obj)


}

