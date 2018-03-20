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


def get_cluster_meta(names, paths):
  if len(names) != len(paths):
    raise ValueError('Number of cluster names must equal length of cluster paths')

  cluster_meta = {}

  for i, name in names:
    cluster_meta[name] = paths[i]

  return cluster_meta


def write_ideogram_annots(infercnv_output, cluster_meta):

  return


if __name__ == '__main__':

  # Parse command-line arguments
  ap = ArgumentParser(description=__doc__,  # Use text from file summary up top
                      formatter_class=RawDescriptionHelpFormatter)
  ap.add_argument('infercnv_output',
                  help='Path to pre_vis_transform.txt output from inferCNV')
  ap.add_argument('cluster_names',
                  help='List of cluster names',
                  nargs='+')  # List must have one or more items
  ap.add_argument('cluster_paths',
                  help='List of cluster paths or URLs',
                  nargs='+')

  args = ap.parse_args()

  infercnv_output = args.infercnv_output
  cluster_names = args.cluster_names
  cluster_paths = args.cluster_paths

  cluster_meta = get_cluster_meta(cluster_names, cluster_paths)

  write_ideogram_annots(infercnv_output, cluster_meta)
