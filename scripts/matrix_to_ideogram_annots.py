#!/usr/bin/env python


"""Converts clustered gene expression matrices to Ideogram.js annotations
"""

__author__ = 'Eric Weitz, Jonathan Bistline, Timothy Tickle'
__copyright__ = 'Copyright 2018'
__credits__ = ['Eric Weitz']
__license__ = 'BSD-3'
__maintainer__ = 'Eric Weitz'
__email__ = 'eweitz@bbroadinstitute.org'
__status__ = 'Development'

from argparse import ArgumentParser, RawDescriptionHelpFormatter


if __name__ == '__main__':

  # Parse command-line arguments
  ap = ArgumentParser(description=__doc__,  # Use text from file summary up top
                      formatter_class=RawDescriptionHelpFormatter)
  ap.add_argument('infercnv_output',
                  help='Path to pre_vis_transform.txt output from inferCNV')
  ap.add_argument('ordination_names',
                  help='List of cluster ordination names',
                  nargs='+')  # List must have one or more items
  ap.add_argument('ordination_paths',
                  help='List of cluster ordination paths or URLs',
                  nargs='+')

  args = ap.parse_args()

  infercnv_output = args.infercnv_output
  ordination_names = args.ordination_names
  ordination_paths = args.ordination_paths