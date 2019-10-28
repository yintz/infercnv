#!/usr/bin/env Rscript


library(optparse)
library(Matrix)
library(data.table)
library(logging)


logging::basicConfig(level='INFO')

pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--input data_matrix ",
                                            "--output sparse_matrix ",
                                            "--delim matrix_delimiter"
                                            ))

pargs <- optparse::add_option(pargs, c("--input"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="input",
                              metavar="input",
                              help=paste("Input raw counts matrix ",
                              			 "to prepare for infercnv run."))

pargs <- optparse::add_option(pargs, c("--output"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="output",
                              metavar="output",
                              help=paste("Output raw counts matrix ",
                              			 "as a sparseMatrix for infercnv run."))

pargs <- optparse::add_option(pargs, c("--delim"),
                              type="character",
                              action="store",
                              default="\t",
                              dest="delim",
                              metavar="delim",
                              help=paste("Delimiter for reading expression matrix",
                                         "[Default %default]"))

args <- optparse::parse_args(pargs)

if (is.null(args$input) || is.null(args$output)) {
	logging::logerror("Please provide input and output arguments")
}

logging::loginfo("Reading header.")

data_head = fread(input=args$input,
	  sep=args$delim,
	  header=FALSE,
	  nrows=1,
	  stringsAsFactors=FALSE,
	  check.names=FALSE,
	  nThread=1,
	  logical01=FALSE,
	  data.table=FALSE)

logging::loginfo("Done reading header.")
logging::loginfo("Reading matrix data.")

ddata = fread(input=args$input,
	  sep=args$delim,
	  header=FALSE,
	  skip=1,
	  stringsAsFactors=FALSE,
	  check.names=FALSE,
	  nThread=1,
	  logical01=FALSE,
	  data.table=FALSE)

logging::loginfo("Done reading matrix data.")

logging::loginfo("Backing up rownames.")
# store column names before dropping the column from the matrix
saved_names = as.vector(unlist(ddata[, 1]))
ddata = ddata[, -1, drop=FALSE]

in_size = object.size(ddata)

colnames(ddata) = as.vector(unlist(data_head))

logging::loginfo("Converting data.frame to Matrix.")
basic_matrix = as.matrix(ddata)
logging::loginfo("Done converting data.frame to Matrix.")
logging::loginfo("Freeing data.frame.")
rm(ddata)  # make memory available
gc()
logging::loginfo("Converting Matrix to sparseMatrix.")
sparse_matrix = Matrix(basic_matrix, sparse=T)
logging::loginfo("Done converting Matrix to sparseMatrix.")
logging::loginfo("Freeing Matrix.")
rm(basic_matrix)  # make memory available
gc()
logging::loginfo("Setting rownames.")
row.names(sparse_matrix) = saved_names

logging::loginfo("Saving sparseMatrix to RDS file.")
saveRDS(sparse_matrix, file=paste(args$output, "rds", sep="."))

out_size = object.size(sparse_matrix)

fileConn<-file("prepare_smallest.txt")
if (in_size < out_size) {
	writeLines(args$input, fileConn)
} else {
	writeLines(paste(args$output, "rds", sep="."), fileConn)	
}
close(fileConn)

