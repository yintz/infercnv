#!/usr/bin/env python


"""
Coverts a square expression matrix to an R-format compatible expression matrix
"""


# Import statements
import argparse
import csv
import os

__author__ = 'Jon Bistline'
__copyright__ = 'Copyright 2018'
__credits__ = ["Jon Bistline"]
__license__ = 'BSD-3'
__maintainer__ = 'Jon Bistline'
__email__ = 'bistline@broadinstitute.org'
__status__ = 'Development'

def convert_matrix_format(input_matrix, delimiter, output_name):
    """ Convert input expression matrix to R-compatible expression matrix (header line is 1 cell shorter than data lines)

    :param input_matrix: Path to input expression matrix
    :type input_matrix: String
    :param delimiter: delimiter to parse input matrix with (tab, comma, etc.)
    :type delimiter: String

    """

    if not input_matrix or not os.path.exists(input_matrix):
        print("".join(["check_matrix_format.py:: ",
                       "Could not find input matrix : " + input_matrix]))

    # read first line
    with open(input_matrix, "r") as exp_matrix:
        print("".join(["Opening input matrix and checking header format: ", input_matrix]))
        print("".join(["Using delimiter: ", delimiter]))
        rewrite_file = False
        matrix = csv.reader(exp_matrix, delimiter=delimiter)
        header_list = next(matrix)
        # check if first value in header_list needs to be removed
        headers_to_remove = ['GENE', 'gene', '']
        if header_list[0] in headers_to_remove:
            print("Input matrix is being converted to R format.")
            rewrite_file = True
            header_list.pop(0)
            with open(output_name, 'w+') as new_expression_matrix:
                writer = csv.writer(new_expression_matrix, delimiter=delimiter)
                writer.writerow(header_list)
                for line in matrix:
                    writer.writerow(line)

    if rewrite_file is True:
        print("".join(["Conversion complete, new output file: ", output_name]))
    else:
        os.rename(input_matrix, output_name)
        print("".join(["No conversion necessary, input matrix is in R format already, renamed to new output file: ", output_name]))

if __name__ == "__main__":

    # Parse arguments
    prsr_arguments = argparse.ArgumentParser(prog='check_matrix_format.py',
                                             description='Coverts a square expression matrix to an R-compatible expression matrix.',
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Add positional argument
    prsr_arguments.add_argument("--input_matrix",
                                metavar="input_matrix",
                                help="Path to the input expression matrix")
    prsr_arguments.add_argument("--delimiter",
                                metavar="delimiter",
                                default="\t",
                                help="delimiter to parse input matrix with (tab, comma, etc.)")
    prsr_arguments.add_argument("--output_name",
                                metavar="output_name",
                                default="expression.r_format.txt",
                                help="path to output expression matrix")
    args = prsr_arguments.parse_args()

    # Run Script
    convert_matrix_format(args.input_matrix, args.delimiter, args.output_name)
