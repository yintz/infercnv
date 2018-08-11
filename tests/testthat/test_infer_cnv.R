# Global data

matrix_zeros <- matrix(rep(0,5), ncol=1)

matrix_one <- matrix(1:5, ncol=1)

matrix_one_long <- matrix(1:20, ncol=1)

matrix_one_long_2 <- matrix(c(1,2,4,7,9,11,12,14,17,19,16,14,
                              13,11,10,7,6,4,3,1), ncol=1)
matrix_two <- matrix(1:10, ncol=2)

matrix_two_long <- matrix(1:40, ncol=2)

matrix_two_long_2 <- matrix(c(1,2,4,7,9,11,12,14,17,19,16,14,13,11,10,7,6,4,3,
                              1,1,2,4,7,9,11,12,14,17,19,16,14,13,11,10,7,6,4,
                              3,1), ncol=2)

matrix_three <- matrix(1:15, ncol=3)

matrix_five <- matrix(1:25, ncol=5)

context("Test subtract_ref")

matrix_averef_five <- matrix(c(c(-101, -100, -100, -100, -99),
                               c(-101, -100, -99, -98, -99),
                               c(1, 1, 2, 3, 0),
                               c(110, 103, 90, 80, 70),
                               c(0, 0, 0, 0, 0),
                               c(100, 102, 100, 102, 102),
                               c(0, -1, -4, -1, -1),
                               c(105, 95, 80, 97, 80),
                               c(100, 99, 100, 101, 100),
                               c(0, 0, 0, 0, 0)),
                               ncol=5, byrow=FALSE)

avref_answer_1 <- matrix(0:4, ncol=1)
avref_answer_2 <- matrix(c(0:4,0:4), ncol=2)
avref_answer_3 <- matrix(c(-1:3, -1:3, -1:3), ncol=3)
avref_answer_4 <- matrix(rep(-3:1 + .5, 5),ncol=5)
avref_answer_5 <- matrix_zeros
matrix_averef_five_answer <- matrix(c(c(-1,0,0,0,0,-1,0,0,1,0),
                                      c(0,0,0,0,-1,40,33,20,10,0),
                                      c(0,0,0,0,0,0,0,0,0,0),
                                      c(0,0,-3,0,0,25,15,0,17,0),
                                      c(1,0,1,2,1,0,0,0,0,0)),
                                    ncol=10,
                                    byrow=TRUE)

test_that("subtract_ref works with one observation, one reference",{
    expect_equal(subtract_ref(average_data=t(matrix_one),
                                  #ref_observations=c(1),
                                  ref_groups=list(c(1))),
                 t(avref_answer_1))
          })
test_that("subtract_ref works with two observations, one reference",{
    expect_equal(subtract_ref(average_data=t(matrix_two),
                                  #ref_observations=c(1),
                                  ref_groups=list(c(1))),
                 t(avref_answer_2))
          })
#test_that("subtract_ref works with 3 observations, two reference",{
#    expect_equal(subtract_ref(average_data=t(matrix_three),
#                                  ref_observations=c(1,3),
#                                  ref_groups=list(c(1,2))),
#                 t(avref_answer_3))
#          })
test_that("subtract_ref updated works with 3 observaions, two reference",{
   expect_equal(subtract_ref(average_data=t(matrix_three),
				  #ref_observations=c(1,3),
                                  ref_groups=list(c(1,3))),
                 t(avref_answer_3))
	  })
test_that("subtract_ref works with 5 observations, two reference",{
    expect_equal(subtract_ref(average_data=t(matrix_five),
                                  #ref_observations=c(2,5),
                                  #ref_groups=list(c(1,2))),
				  ref_groups=list(c(2,5))),
                 t(avref_answer_4))
          })
test_that("subtract_ref works with 1 observation, 1 reference",{
    expect_equal(subtract_ref(average_data=t(matrix_zeros),
                                  #ref_observations=c(1),
                                  ref_groups=list(c(1))),
                 t(avref_answer_5))
          })
test_that("subtract_ref works with 10 obs, 5 references, 3 groups",{
    expect_equal(subtract_ref(average_data=t(matrix_averef_five),
                                  #ref_observations=c(2,4,6,8,10),
                                  #ref_groups=list(c(1),c(2,3,4),c(5))),
				  ref_groups=list(c(2),c(4,6,8),c(10))),
                 matrix_averef_five_answer)
          })


context("Test split_references")

split_matrix_one <- NULL
split_matrix_two <- matrix(1:10, ncol=1)
split_matrix_three <- matrix(1:10, ncol=1)
split_matrix_four <- matrix(c(1:10,2:11,31:40,32:41,33:42,0:9),
                            ncol=10, byrow=TRUE)
split_matrix_five <- matrix(c(1:10,2:11,31:40,32:41,33:42,0:9),
                            ncol=10, byrow=TRUE)

split_obs_one <- NULL
split_obs_two <- c(1)
split_obs_three <- c(2)
split_obs_four <- 1:3
split_obs_five <- c(2,4,6)

split_num_grp_one <- NULL
split_num_grp_two <- 1
split_num_grp_three <- 2
split_num_grp_four <- 1
split_num_grp_five <- 3

split_answer_one <- list()
split_answer_two <- list()
split_answer_two[[1]] <- split_obs_two
split_answer_three <- list()
split_answer_three[[1]] <- split_obs_three
split_answer_four <- list()
split_answer_four[[1]] <- split_obs_four
split_answer_five <- list()
split_answer_five[[1]] <- c(2)
split_answer_five[[2]] <- c(4)
split_answer_five[[3]] <- c(6)

test_that("split_references for null matrix.",{
    expect_equal(split_references(average_data=split_matrix_one,
                                  ref_obs=split_obs_one,
                                  num_groups=split_num_grp_one),
                 split_answer_one)
        })

test_that("split_references for one observation matrix, one group",{
    expect_equal(split_references(average_data=split_matrix_two,
                                  ref_obs=split_obs_two,
                                  num_groups=split_num_grp_two),
                split_answer_two)
        })

test_that("split_references for one observation matrix, two groups",{
    expect_equal(split_references(average_data=split_matrix_three,
                                  ref_obs=split_obs_three,
                                  num_groups=split_num_grp_three),
                split_answer_three)
        })

test_that("split_references for three observation matrix, one group",{
    expect_equal(split_references(average_data=split_matrix_four,
                                  ref_obs=split_obs_four,
                                  num_groups=split_num_grp_four),
                split_answer_four)
        })

test_that("split_references for three observation matrix, three groups",{
    expect_equal(split_references(average_data=split_matrix_five,
                                  ref_obs=split_obs_five,
                                  num_groups=split_num_grp_five),
                split_answer_five)
        })


context("Test center_with_threshold")

center_answer_1 <- matrix(rep(0,5), ncol=1)
center_answer_2 <- matrix(c(rep(-5,5),rep(0,5),rep(5,5)), ncol=3)
center_answer_3 <- matrix(c(rep(-6,5),rep(-5,5),rep(0,5),
                            rep(5,5),rep(6,5)), ncol=5)
center_answer_4 <- matrix(rep(0,25), ncol=5)

test_that(paste("center_with_threshold works with one observation,",
                "threshold too large to affect"),{
    expect_equal(center_with_threshold(center_data=matrix_one,
                                       threshold=100),
                 center_answer_1)
          })
test_that(paste("center_with_threshold works with three observations,",
                "threshold too large to affect"),{
    expect_equal(center_with_threshold(center_data=matrix_three,
                                       threshold=100),
                 center_answer_2)
          })
test_that(paste("center_with_threshold works with five observations,",
                "threshold affecting some"),{
    expect_equal(center_with_threshold(center_data=matrix_five,
                                       threshold=6),
                 center_answer_3)
         })
test_that(paste("center_with_threshold works with one observation,",
                "threshold of 0, affecting all"),{
    expect_equal(center_with_threshold(center_data=matrix_five,
                                       threshold=0),
                 center_answer_4)
         })


context("Test center_smoothed")

center_sm_1 <- matrix(1:10, ncol=1)
center_sm_1_answer <- matrix(c(-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5),
                             ncol=1)
center_sm_3 <- matrix(1:21, ncol=3)
center_sm_3_answer <- matrix(rep(c(-3,-2,-1,0,1,2,3),3), ncol=3)

test_that("center_smoothed works with 1 observations",{
    expect_equal(center_smoothed(data_smoothed=center_sm_3),
                 center_sm_3_answer)
    })

test_that("center_smoothed works with 3 observations",{
    expect_equal(center_smoothed(data_smoothed=center_sm_3),
                 center_sm_3_answer)
    })


context("Test above_min_mean_expr_cutoff")

above_answer_1 <- 1:5
above_answer_2 <- 1:5
above_answer_3 <- 4:5
above_answer_4 <- c(5)
above_answer_5 <- NULL
above_answer_6 <- NULL

test_that(paste("above_min_mean_expr_cutoff works with one observation,",
                "cutoff too large to affect"),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_one + 1),
                              cutoff=0),
                 above_answer_1)
         })
test_that(paste("above_min_mean_expr_cutoff works with three observations,",
                "threshold too large to affect"),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_three + 1),
                              cutoff=0),
                 above_answer_2)
         })
test_that(paste("above_min_mean_expr_cutoff works with one observation,",
                "threshold excluding two."),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_one + 1),
                              cutoff=2),
                 above_answer_3)
         })
test_that(paste("above_min_mean_expr_cutoff works with three observations,",
                "threshold excluding three."),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_three + 1),
                              cutoff=3.4),
                 above_answer_4)
         })
test_that(paste("above_min_mean_expr_cutoff works with one observation,",
                "threshold excluding all."),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_one + 1),
                              cutoff=100),
                 above_answer_5)
         })
test_that(paste("above_min_mean_expr_cutoff works with three observations,",
                "threshold excluding all."),{
    expect_equal(above_min_mean_expr_cutoff(data=log2(matrix_three + 1),
                              cutoff=100),
                 above_answer_6)
         })


context("Test remove_noise")

noise_answer_1 <- matrix_one
noise_answer_2 <- matrix(c(0,0,0,4,5), ncol=1)
noise_answer_3 <- matrix_zeros
noise_answer_4 <- matrix_three
noise_answer_5 <- matrix(c(rep(0,11),12:15), ncol=3)
noise_answer_6 <- matrix(rep(0,15), ncol=3)

test_that("remove_noise works with one observation, threshold 0",{
    expect_equal(remove_noise(smooth_matrix=matrix_one,
                              threshold=0),
                  noise_answer_1)
         })
test_that(paste("remove_noise works with one observation, one ref,",
                "threshold removing some"),{
    expect_equal(remove_noise(smooth_matrix=matrix_one,
                              threshold=4),
                  noise_answer_2)
         })
test_that(paste("remove_noise works with one observation,",
                "threshold removing all"),{
    expect_equal(remove_noise(smooth_matrix=matrix_one,
                              threshold=6),
                  noise_answer_3)
         })
test_that("remove_noise works with three observation, threshold 0",{
    expect_equal(remove_noise(smooth_matrix=matrix_three,
                              threshold=0),
                  noise_answer_4)
         })
test_that("remove_noise works with three observation, threshold some",{
    expect_equal(remove_noise(smooth_matrix=matrix_three,
                              threshold=12),
                  noise_answer_5)
         })
test_that("remove_noise works with three observation, threshold all",{
    expect_equal(remove_noise(smooth_matrix=matrix_three,
                              threshold=100),
                  noise_answer_6)
         })


context("Test remove_tails")

tail_answer_1 <- c()
tail_answer_2 <- -1 * c(1:5, 16:20)
tail_answer_3 <- -1 * c(2:6, 13:17)
tail_answer_4 <- -1 * c(5:9, 11:15)
tail_answer_5 <- -1 * c(1, 5)

test_that(paste("remove tails works with one contig,",
                "one observation, no tail length"),{
    expect_equal(remove_tails(smooth_matrix=matrix_one,
                              chr=1:length(matrix_one),
                              tail_length=0),
                 tail_answer_1)
         })
test_that("remove tails works with one contig, one observation, tail length 5",{
    expect_equal(remove_tails(smooth_matrix=matrix_one_long,
                              chr=1:length(matrix_one_long),
                              tail_length=5),
                 tail_answer_2)
         })
test_that("remove tails works with 3 contig, one observation, tail length 5",{
    expect_equal(remove_tails(smooth_matrix=matrix_one_long,
                              chr=2:17,
                              tail_length=5),
                 tail_answer_3)
         })
test_that("remove tails works with 3 contig, two observations, tail length 5",{
    expect_equal(remove_tails(smooth_matrix=matrix_two_long,
                              chr=5:15,
                              tail_length=5),
                 tail_answer_4)
         })
test_that(paste("remove tails works with one contig, one observation,",
                "tail length longer than contig"),{
    expect_equal(remove_tails(smooth_matrix=matrix_one,
                              chr=1:length(matrix_one),
                              tail_length=100),
                 tail_answer_5)
         })


context("smooth_window")

smooth_answer_1 <- matrix_one

smooth_answer_2 <- matrix_one

smooth_answer_3 <- matrix(c(1.00,2.53,4.60,6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.73,1.00), ncol=1)

smooth_answer_4 <- matrix(c(1.00,2.53,4.60,6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.73,1.00,
                            1.00,2.53,4.60,
                            6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.73,1.00), ncol=2)

smooth_answer_5 <- matrix_one

test_that("smooth_window works with one observation, window_length 0",{
    expect_equal(smooth_window(data=matrix_one,
                               window_length=0,
                               smooth_ends=F,
                               re_center=F),
                 smooth_answer_1)
         })
test_that("smooth_window works with one observation, window_length 1",{
    expect_equal(smooth_window(data=matrix_one,
                               window_length=1,
                               smooth_ends=F,
                               re_center=F),
                 smooth_answer_2)
         })
test_that("smooth_window works with one observation, window_length 5",{
    expect_equal(round(smooth_window(data=matrix_one_long_2,
                                     window_length=5,
                                     smooth_ends=T,
                                     re_center=F),2),
                 smooth_answer_3)
         })
test_that("smooth_window works with two observations, window_length 5",{
    expect_equal(round(smooth_window(data=matrix_two_long_2,
                                     window_length=5,
                                     smooth_ends=T,
                                     re_center=F),2),
                 smooth_answer_4)
         })
test_that(paste("smooth_window works with one observation,",
                "window_length longer than measurements"),{
                    expect_equal(smooth_window(data=matrix_one,
                                               window_length=100,
                                               smooth_ends=F,
                                               re_center=F),
                 smooth_answer_5)
         })


# test smooth_ends_helper
ends_data = rep(c(1,2,3), 10)
smooth_ends_tail_3_ans = c(1.0, 2.0, 1.8, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0,
                           3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0,
                           2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 2.2, 2.0, 3.0)

test_that("smooth_ends_helper works, tail length 3", {
    expect_equal(round(smooth_ends_helper(obs_data=ends_data, tail_length=3),2),
                 smooth_ends_tail_3_ans)
    })


smooth_ends_tail_7_ans = c(1.00, 2.00, 1.80, 1.86, 2.00, 1.91, 1.92, 2.00,
                           3.00, 1.00, 2.00, 3.00, 1.00, 2.00, 3.00,
                           1.00, 2.00, 3.00, 1.00, 2.00, 3.00, 1.00,
                           2.00, 2.08, 2.09, 2.00, 2.14, 2.20, 2.00, 3.00)

test_that("smooth_ends_helper works, tail length 7", {
    expect_equal(round(smooth_ends_helper(obs_data=ends_data, tail_length=7),2),
                 smooth_ends_tail_7_ans)
    })



context("create_sep_list")
create_sep_list_answer_1 <- list()
create_sep_list_answer_1[[1]] <- list()
create_sep_list_answer_1[[1]][[1]] <- c(0,0,0,0)
create_sep_list_answer_1[[2]] <- list()
create_sep_list_answer_1[[2]][[1]] <- c(0,0,0,0)

create_sep_list_answer_2 <- list()
# Column
create_sep_list_answer_2[[1]] <- list()
create_sep_list_answer_2[[1]][[1]] <- c(2,0,2,5)
create_sep_list_answer_2[[1]][[2]] <- c(8,0,8,5)
# Row
create_sep_list_answer_2[[2]] <- list()
create_sep_list_answer_2[[2]][[1]] <- c(0,4,9,4)
create_sep_list_answer_2[[2]][[2]] <- c(0,2,9,2)
 
test_that("create_sep_list works with 0 row and column and no seps",{
    expect_equal(create_sep_list(0, 0),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 0 row and column and seps",{
    expect_equal(create_sep_list(0, 0, 1:5, 3:6),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 10 row, 0 column and no seps",{
    expect_equal(create_sep_list(10, 0),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 0 row, 10 column and seps",{
    expect_equal(create_sep_list(0, 10),
                 create_sep_list_answer_1)
})
test_that("create_sep_list works with 5 row, 9 column and seps",{
    expect_equal(create_sep_list(row_count=5,
                                 col_count=9,
                                 row_seps=c(1,3),
                                 col_seps=c(2,8)),
                 create_sep_list_answer_2)
})
 
context("remove_outliers_norm")
remove_outlier_norm_in_1 <- matrix(1:20, ncol=4)
remove_outlier_norm_out_1 <- matrix(c(rep(5,5),6:14,rep(15,6)), ncol=4)
remove_outlier_norm_in_2 <- matrix(c(1:15,
                                     c(-5,-4,3:13,21,26),
                                     1:15,
                                     1:15), ncol=4)
remove_outlier_norm_out_2 <- matrix(c(1:15,
                                    c(-.5,-.5,3:13,17.75,17.75),
                                    1:15,
                                    1:15), ncol=4)

test_that("remove_outliers_norm for no change",{
    expect_equal(remove_outliers_norm( remove_outlier_norm_in_1 ),
                 remove_outlier_norm_in_1)
})
test_that("remove_outliers_norm for hard threshold outside of values",{
    expect_equal(remove_outliers_norm(remove_outlier_norm_in_1,
                                      lower_bound=-1,
                                      upper_bound=30),
                 remove_outlier_norm_in_1)
})
test_that("remove_outliers_norm for hard threshold within values",{
    expect_equal(remove_outliers_norm(remove_outlier_norm_in_1,
                                      lower_bound=5,
                                      upper_bound=15),
                 remove_outlier_norm_out_1)

})
test_that("remove_outliers_norm for average bound",{
    expect_equal(remove_outliers_norm(remove_outlier_norm_in_2,
                                      out_method="average_bound"),
                 remove_outlier_norm_out_2)
})

context("plot_cnv_observations")
obs_data <- matrix(1:20, ncol=4)
row.names(obs_data) <- paste("gene",1:5, sep="_")
colnames(obs_data) <- paste("Sample",1:4,sep="_")
obs_file_base_name <- "."
obs_contig_colors <- c("red","red","cyan","yellow","yellow")
names(obs_contig_colors) <- c(1,1,2,3,3)
obs_contig_label <- c("1","","2","3","") 
obs_contig_name <- c("1","1","2","2","3")
obs_col_pal <- rainbow
obs_contig_seps <- c(2,3)
obs_groups <- 2
names(obs_contig_seps) <- c(4,6)
obs_plot=FALSE
obs_answer_1 <-obs_data[c("gene_1","gene_2","gene_5","gene_3","gene_4"),]
obs_answer_3 <-obs_data[c("gene_1","gene_2"),]
obs_answer_4 <-obs_data[c("gene_3","gene_4","gene_5"),]
test_plot_cnv_file <- function(obs_data,
                              obs_col_pal,
                              obs_contig_seps,
                              obs_file_base_name,
                              obs_file){
                                  plot_cnv_observations(obs_data=obs_data,
                                                        file_base_name=obs_file_base_name,
                                                        contig_colors=obs_contig_colors,
                                                        contig_label=obs_contig_label,
                                                        contig_names=obs_contig_name,
                                                        col_pal=obs_col_pal,
                                                        contig_seps=obs_contig_seps,
                                                        num_obs_groups=obs_groups,
                                                        testing=!obs_plot);
                                  return(as.matrix(read.table(file.path(obs_file_base_name,obs_file),
                                                   row.names=1)))}

# TODO, turn on. The is creates a file, turned off normally as a result.
#test_that("plot_cnv_observations test for accurate output of observation file.",{
#    expect_equal(test_plot_cnv_file(obs_data=obs_data,
#                                   obs_col_pal=obs_col_pal,
#                                   obs_contig_seps=obs_contig_seps,
#                                   obs_file_base_name=obs_file_base_name,
#                                   obs_file="observations.txt"),
#                 obs_answer_1)
#})

# TODO, turn on. The is creates a file, turned off normally as a result.
#test_that("plot_cnv_observations test for accurate output of obs group 1 file.",{
#    expect_equal(test_plot_cnv_file(obs_data=obs_data,
#                                       obs_col_pal=obs_col_pal,
#                                       obs_contig_seps=obs_contig_seps,
#                                       obs_file_base_name=obs_file_base_name,
#                                       obs_file="General_HCL_1_members.txt"),
#                 obs_answer_3)
#})

# TODO, turn on. The is creates a file, turned off normally as a result.
#test_that("plot_cnv_observations test for accurate output of obs group 2 file.",{
#    expect_equal(test_plot_cnv_file(obs_data=obs_data,
#                                       obs_col_pal=obs_col_pal,
#                                       obs_contig_seps=obs_contig_seps,
#                                       obs_file_base_name=obs_file_base_name,
#                                       obs_file="General_HCL_2_members.txt"),
#                 obs_answer_4)
#})

context("plot_cnv_references")
references_ref_data <- matrix(1:20, ncol=4)
row.names(references_ref_data) <- paste("gene",1:5, sep="_")
colnames(references_ref_data) <- paste("Sample",1:4,sep="_")
references_ref_groups <- list()
references_ref_groups[[1]] <- c(3)
references_ref_groups[[1]] <- c(1,2,4)
references_col_pal <- rainbow
references_contig_seps <- c(2,3)
names(references_contig_seps) <- c(4,6)
references_file_base_name <- "."
ref_plot=FALSE
test_plot_cnv_references <- function(references_ref_data,
                                     references_ref_groups,
                                     references_col_pal,
                                     references_contig_seps,
                                     references_file_base_name){
                                         plot_cnv_references(ref_data=references_ref_data,
                                             ref_groups=references_ref_groups,
                                             col_pal=references_col_pal,
                                             contig_seps=references_contig_seps,
                                             file_base_name=references_file_base_name,
                                             testing=!ref_plot);
                                         return(as.matrix(read.table(file.path(references_file_base_name,"references.txt"),
                                                           row.names=1)))}

# TODO, turn on. The is creates a file, turned off normally as a result.
#test_that("plot_cnv_references test for accurate output reference file.",{
#    expect_equal(test_plot_cnv_references(references_ref_data=references_ref_data,
#                          references_ref_groups=references_ref_groups,
#                          references_col_pal=references_col_pal,
#                          references_contig_seps=references_contig_seps,
#                          references_file_base_name=references_file_base_name),
#                 references_ref_data[,ncol(references_ref_data):1])
#})

context("order_reduce")
order_reduce_data_1 <- matrix(rep(1:10,2), ncol=2)
colnames(order_reduce_data_1) <- c("Sample_1","Sample_2")
row.names(order_reduce_data_1) <- paste("gene",1:10,sep="_") 
order_reduce_exp_1 <- matrix(rep(c(10,5,8,3,4,9,1,7,6,2),2), ncol=2)
row.names(order_reduce_exp_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")
colnames(order_reduce_exp_1) <- c("Sample_1","Sample_2")
order_reduce_pos_1 <- data.frame(chr=factor(c(1,1,2,2,3,3,4,4,5,5),levels=1:5),
                                 start=c(1,5,1,5,1,5,1,5,1,5),
                                 stop=c(4,9,4,9,4,9,4,9,4,9))
row.names(order_reduce_pos_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")
order_reduce_chr_1 <- data.frame(chr=factor(c(1,1,2,2,3,3,4,4,5,5),levels=1:5))
row.names(order_reduce_chr_1) <- paste("gene",c(10,5,8,3,4,9,1,7,6,2),sep="_")

order_reduce_pos_2 <- data.frame(chr=factor(c(1,1,2,3,4,4),levels=1:4),
                                 start=c(1,5,5,5,1,5),
                                 stop=c(4,9,9,9,4,9))
row.names(order_reduce_pos_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")
order_reduce_exp_2 <- matrix(rep(c(10,5,3,9,1,7),2), ncol=2)
row.names(order_reduce_exp_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")
colnames(order_reduce_exp_2) <- c("Sample_1","Sample_2")
order_reduce_chr_2 <- data.frame(chr=factor(c(1,1,2,3,4,4),levels=1:4))
row.names(order_reduce_chr_2) <- paste("gene",c(10,5,3,9,1,7),sep="_")

order_reduce_pos_3 <- data.frame(chr=c(1,1,2,3,4,4),
                                 start=c(1,5,5,5,1,5),
                                 stop=c(4,9,9,9,4,9))
row.names(order_reduce_pos_3) <- paste("GENE",c(10,5,3,9,1,7),sep="_")

test_that("order_reduce for NULL input.",{
    expect_equal(order_reduce(NULL,NULL),
                 list(expr=NULL,order=NULL,chr_order=NULL))
})
test_that("order_reduce for happy path",{
    expect_equal(order_reduce(order_reduce_data_1, order_reduce_pos_1),
                 list(expr=order_reduce_exp_1,
                      order=order_reduce_pos_1,
                      chr_order=order_reduce_chr_1))
})
test_that("order_reduce for happy path but dropping genes",{
    expect_equal(order_reduce(order_reduce_data_1, order_reduce_pos_2),
                 list(expr=order_reduce_exp_2,
                      order=order_reduce_pos_2,
                      chr_order=order_reduce_chr_2))
})
test_that("order_reduce for no matching gene names",{
    expect_equal(order_reduce(order_reduce_data_1, order_reduce_pos_3),
                 list(expr=NULL,
                      order=NULL,
                      chr_order=NULL))
})
