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

context("Test average_over_ref")

avref_answer_1 <- matrix_one
avref_answer_2 <- matrix(c(rep(0,5), rep(5,5)), ncol=2)
avref_answer_3 <- matrix(c(rep(-5,5),
                           rep(0,5),
                           rep(5,5)), ncol=3)
avref_answer_4 <- matrix(c(rep(-12.5,5),
                           rep(-7.5,5),
                           rep(-2.5,5),
                           rep(2.5,5),
                           rep(7.5,5)), ncol=5)
avref_answer_5 <- matrix_zeros

test_that("average_over_ref works with one observation, one reference",{
    expect_equal(average_over_ref(average_data=matrix_one,
                                  ref_observations=c(1)),
                 matrix(rep(0, nrow(matrix_one)),ncol=1))
          })
test_that("average_over_ref works with two observations, one reference",{
    expect_equal(average_over_ref(average_data=matrix_two,
                                  ref_observations=c(1)),
                 avref_answer_2)
          })
test_that("average_over_ref works with 3 observations, two reference",{
    expect_equal(average_over_ref(average_data=matrix_three,
                                  ref_observations=c(1,3)),
                 avref_answer_3)
          })
test_that("average_over_ref works with 5 observations, two reference",{
    expect_equal(average_over_ref(average_data=matrix_five,
                                  ref_observations=c(2,5)),
                 avref_answer_4)
          })
test_that("average_over_ref works with 1 observation, 1 reference",{
    expect_equal(average_over_ref(average_data=matrix_zeros,
                                  ref_observations=c(1)),
                 avref_answer_5)
          })

context("Test split_references")

split_matrix_one <- NULL
split_matrix_two <- matrix(1:10, ncol=1)
split_matrix_three <- matrix(1:10, ncol=1)
split_matrix_four <- matrix(c(1:10,2:11,31:40,32:41,33:42,0:9), ncol=10, byrow=TRUE)
split_matrix_five <- matrix(c(1:10,2:11,31:40,32:41,33:42,0:9), ncol=10, byrow=TRUE)

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

# Add a test here to be the full test case on the demo data we are using.
# infer_cnv <- function(data, gene_order, cutoff, reference_obs,
#                       window_length, max_centered_threshold,
#                       noise_threshold, pdf_path){

context("Test above_cutoff")

above_answer_1 <- 1:5
above_answer_2 <- 1:5
above_answer_3 <- 3:5
above_answer_4 <- 4:5
above_answer_5 <- NULL
above_answer_6 <- NULL

test_that(paste("above_cutoff works with one observation,",
                "cutoff too large to affect"),{
    expect_equal(above_cutoff(data=log2(matrix_one+1),
                              cutoff=0),
                 above_answer_1)
         })
test_that(paste("above_cutoff works with three observations,",
                "threshold too large to affect"),{
    expect_equal(above_cutoff(data=log2(matrix_three+1),
                              cutoff=0),
                 above_answer_2)
         })
test_that(paste("above_cutoff works with one observation,",
                "threshold excluding two."),{
    expect_equal(above_cutoff(data=log2(matrix_one+1),
                              cutoff=2),
                 above_answer_3)
         })
test_that(paste("above_cutoff works with three observations,",
                "threshold excluding three."),{
    expect_equal(above_cutoff(data=log2(matrix_three+1),
                              cutoff=8),
                 above_answer_4)
         })
test_that(paste("above_cutoff works with one observation,",
                "threshold excluding all."),{
    expect_equal(above_cutoff(data=log2(matrix_one+1),
                              cutoff=100),
                 above_answer_5)
         })
test_that(paste("above_cutoff works with three observations,",
                "threshold excluding all."),{
    expect_equal(above_cutoff(data=log2(matrix_three+1),
                              cutoff=100),
                 above_answer_6)
         })

context("Test remove_noise")

noise_answer_1 <- matrix_one
noise_answer_2 <- matrix(c(0,0,0,4,5), ncol=1)
noise_answer_3 <- matrix_zeros
noise_answer_4 <- matrix_three
noise_answer_5 <- matrix(c(1:5,rep(0,6),12:15), ncol=3)
noise_answer_6 <- matrix(c(rep(0,10),11:15), ncol=3)

test_that("remove_noise works with one observation, one ref, threshold 0",{
    expect_equal(remove_noise(ref=1,
                              smooth_matrix=matrix_one,
                              threshold=0),
                  noise_answer_1)
         })
test_that(paste("remove_noise works with one observation, one ref,",
                "threshold removing some"),{
    expect_equal(remove_noise(ref=1,
                              smooth_matrix=matrix_one,
                              threshold=4),
                  noise_answer_2)
         })
test_that(paste("remove_noise works with one observation, one ref,",
                "threshold removing all"),{
    expect_equal(remove_noise(ref=1,
                              smooth_matrix=matrix_one,
                              threshold=6),
                  noise_answer_3)
         })
test_that("remove_noise works with three observation, one ref, threshold 0",{
    expect_equal(remove_noise(ref=2,
                              smooth_matrix=matrix_three,
                              threshold=0),
                  noise_answer_4)
         })
test_that("remove_noise works with three observation, one ref, threshold some",{
    expect_equal(remove_noise(ref=1,
                              smooth_matrix=matrix_three,
                              threshold=12),
                  noise_answer_5)
         })
test_that("remove_noise works with three observation, one ref, threshold all",{
    expect_equal(remove_noise(ref=3,
                              smooth_matrix=matrix_three,
                              threshold=100),
                  noise_answer_6)
         })

context("Test remove_tails")

tail_answer_1 <- matrix_one
tail_answer_2 <- matrix(c(rep(0,5),6:15,rep(0,5)), ncol=1)
tail_answer_3 <- matrix(c(1,rep(0,5),7:12,rep(0,5),18:20), ncol=1)
tail_answer_4 <- matrix(c(1:4,rep(0,5),10,rep(0,5),16:24,
                              rep(0,5),30,rep(0,5),36:40),ncol=2)
tail_answer_5 <- matrix_zeros

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
smooth_answer_3 <- matrix(c(1.00,2.33,4.60,6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.67,1.00), ncol=1)
smooth_answer_4 <- matrix(c(1.00,2.33,4.60,6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.67,1.00,1.00,2.33,4.60,
                            6.60,8.60,10.60,12.60,14.60,
                            15.60,16.00,15.80,14.60,12.80,11.00,9.40,
                            7.60,6.00,4.20,2.67,1.00), ncol=2)
smooth_answer_5 <- matrix_one

test_that("smooth_window works with one observation, window_length 0",{
    expect_equal(smooth_window(data=matrix_one,
                               window_length=0),
                 smooth_answer_1)
         })
test_that("smooth_window works with one observation, window_length 1",{
    expect_equal(smooth_window(data=matrix_one,
                               window_length=1),
                 smooth_answer_2)
         })
test_that("smooth_window works with one observation, window_length 5",{
    expect_equal(round(smooth_window(data=matrix_one_long_2,
                                     window_length=5),2),
                 smooth_answer_3)
         })
test_that("smooth_window works with two observations, window_length 5",{
    expect_equal(round(smooth_window(data=matrix_two_long_2,
                                     window_length=5),2),
                 smooth_answer_4)
         })
test_that(paste("smooth_window works with one observation,",
                "window_length longer than measurements"),{
    expect_equal(smooth_window(data=matrix_one,
                               window_length=100),
                 smooth_answer_5)
         })

#contest("order_reduce")
#order_reduce_pos_1 = 
#order_reduce_exp_1=
#test_that("order_reduce works with null data.",{
#    expect_equal(order_reduce(NULL, order_reduce_pos_1), NULL)
#})
#
#test_that("order_reduce works with null position.",{
#    expect_equal(order_reduce(order_reduce_exp_1, NULL), NULL)
#})


