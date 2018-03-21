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
import json


class MatrixToIdeogramAnnots:

    def __init__(self, infercnv_output, cluster_meta, output_path,
                 genomic_position_file_path):
        """Class and parameter docs in module summary and argument parser"""

        self.infercnv_output = infercnv_output
        self.cluster_meta = cluster_meta
        self.output_path = output_path
        self.genomic_position_file_path = genomic_position_file_path

        self.genes = self.get_genes()

        self.write_annots()

    def write_annots(self):
        """Write Ideogram.js annotations JSON file"""

        ideogram_annots = self.get_annots()

        ideogram_annots_json = json.dumps(ideogram_annots)

        with open(self.output_path) as f:
            f.write(ideogram_annots_json)

    def get_annots(self):
        """Get Ideogram.js annotations"""

        annots = []


        return annots

    def get_genes(self):
        """Convert inferCNV genomic position file into useful 'genes' dict"""

        genes = {}

        with open(self.genomic_position_file_path) as f:
            lines = f.readlines()

        for line in lines:
            columns = line.strip().split()
            id, chr, start, stop = columns
            genes[id] = {
                'id': id,
                'chr': chr,
                'start': start,
                'stop': stop
            }

        return genes

    def get_expression_matrix_dict(self):
        """Parse inferCNV output, return dict of cell expressions by gene"""

        em_dict = {}

        with open(self.infercnv_output) as f:
            lines = f.readlines()

        em_dict['cells'] = lines[0].strip().split(',')
        genes = {}

        for line in lines[1:]:
            columns = line.strip().split(',')
            gene = columns[0]
            expression_by_cell = columns[1:]
            genes[gene] = expression_by_cell

        em_dict['genes'] = genes

        return em_dict

    def compute_gene_expression_means(self):
        """Compute mean expression for each gene across all and each cluster"""

        scores_lists = []

        cluster_names = cluster_meta.keys()

        keys = ['name', 'all'] + cluster_names
        scores_lists.append(keys)

        em_matrix = self.get_expression_matrix_dict()

        gene_expression_lists = em_matrix['genes'].keys()

        for gene in gene_expression_lists:
            gene_exp_list = gene_expression_lists[gene]

        expression = 0.0

        return scores_lists


def get_cluster_meta(names, paths):
    if len(names) != len(paths):
        raise ValueError('Number of cluster names must equal length of cluster paths')

    cluster_meta = {}

    for i, name in names:
        cluster_meta[name] = paths[i]

    return cluster_meta


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
    ap.add_argument('output_path',
                    help='Path or URL to write output to')

    args = ap.parse_args()

    infercnv_output = args.infercnv_output
    cluster_names = args.cluster_names
    cluster_paths = args.cluster_paths
    output_path = args.output_path

    cluster_meta = get_cluster_meta(cluster_names, cluster_paths)

    MatrixToIdeogramAnnots(infercnv_output, cluster_meta, output_path)
