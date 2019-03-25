#!/usr/bin/env Rscript 
# Script to validate inferCNV docker container instances.

#####
# Set up logging
#####

library(logging)
# Logging level choices
C_LEVEL_CHOICES <- names(loglevels)
logging::basicConfig(level='INFO') #initialize to info setting.

#####
# Data sources
#####

## input data for validation (provided in docker image)
infercnv_root <- '/inferCNV/'
validation_input_dir <- paste0(infercnv_root,'example/')
raw_counts_matrix <- paste0(validation_input_dir, 
  'oligodendroglioma_expression_downsampled.counts.matrix')
annotations_file <- paste0(validation_input_dir, 
  'oligodendroglioma_annotations_downsampled.txt')
gene_order_file <- paste0(validation_input_dir, 
  'gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt')
out_dir <- 'output_cli'

## reference output for validation
validation_reference <- paste0(validation_input_dir,'validation/',
  'reference-infercnv.observations.txt')

# Make sure the reference input data exists 
logging::loginfo(paste("Checking for inferCNV validation input files.", sep=""))
if (!file.exists(raw_counts_matrix) || 
    !file.exists(annotations_file)  || 
    !file.exists(gene_order_file)){
    logging::logerror(paste("Missing input file(s)", sep=""))
    stop(paste0('Error: expected input files cannot be found.'))
}


#####
# Run inferCNV for validation
#####
logging::loginfo(paste("Running inferCNV on validation input files.", sep=""))
inferCNV_exe <- paste0(infercnv_root, 'scripts/inferCNV.R')
validate_cmd <- paste0(inferCNV_exe,
                ' --raw_counts_matrix=', raw_counts_matrix,  
                ' --annotations_file=', annotations_file,
                ' --gene_order_file=', gene_order_file,
                ' --ref_group_names=',
                  '\"Microglia/Macrophage,Oligodendrocytes (non-malignant)\"',
                ' --cutoff=1',
                ' --out_dir=', out_dir,
                ' --cluster_by_groups',
                ' --denoise')
logging::loginfo(validate_cmd)
system(validate_cmd)

validation_input <- paste0(out_dir, '/infercnv.observations.txt')

if (!file.exists(validation_input)){
    logging::logerror(paste("Error: expected output file, infercnv.observations.txt, not found.", sep=""))
    stop('Validation aborted - inferCNV analysis on test data failed.\n')
}


#####
# Read in data for validation
#####

ref <- as.matrix(read.csv(validation_reference, header=T, sep = ' '))
obs <- as.matrix(read.csv(validation_input, header=T, sep = ' '))


#####
# Perform validation
#####
logging::loginfo(paste("Performing validation.", sep=""))
if (max ( abs(obs - ref)/abs(ref)) < 1.0e-8){
  unlink(out_dir, recursive=TRUE)
  logging::loginfo(paste("Successful validation - output passes similarity check.", sep=""))
} else { 
  logging::logerror(paste("Error: generated output fails similarity check", sep=""))
  logging::loginfo(paste("Saving validation files in current working directory", sep=""))
  file.copy(validation_reference, "./reference-infercnv.observations.txt")
  file.copy(validation_input, "./infercnv.observations.txt")
  unlink(out_dir, recursive=TRUE)
  stop('Validation failed - max relative error exceeds threshold.\n')
}
