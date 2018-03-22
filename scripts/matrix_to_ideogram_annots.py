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
from statistics import mean


class MatrixToIdeogramAnnots:

    def __init__(self, infercnv_output, gen_pos_path, clusters,
                 output_path):
        """Class and parameter docs in module summary and argument parser"""

        self.infercnv_output = infercnv_output
        self.clusters = clusters
        self.output_path = output_path
        self.genomic_position_file_path = gen_pos_path

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

    def set_cells_by_cluster(self, cells):

        clusters = self.clusters

        cluster_names = clusters.keys()
        for name in clusters:
            clusters[name]['cells'] = []

        for i, cell in enumerate(cells):
            # TODO: Wire in data from paths for self.clusters for real data,
            # this is currently mock data.
            if i % 3 == 0:
                clusters[cluster_names[0]]['cells'].append(cell)
            else:
                clusters[cluster_names[1]]['cells'].append(cell)

        self.clusters = clusters

    def compute_gene_expression_means(self):
        """Compute mean expression for each gene across all and each cluster"""

        scores_lists = []

        cluster_names = self.clusters.keys()

        keys = ['name', 'all'] + cluster_names
        scores_lists.append(keys)

        em_matrix = self.get_expression_matrix_dict()

        cells = em_matrix['cells']
        self.set_cells_by_cluster(cells)
        clusters = self.clusters_meta

        gene_expression_lists = em_matrix['genes'].keys()

        for i, gene in enumerate(gene_expression_lists):
            scores_list = []
            gene_exp_list = gene_expression_lists[gene]
            mean_expression_all = mean(gene_exp_list)
            scores_list.append(gene, mean_expression_all)

            for cluster in clusters:
                expressions = cluster['genes'][gene]
                mean_expression_cluster = mean(expressions)
                scores_list.append(mean_expression_cluster)

            if i % 100 == 0 and i != 0:
                print('Processed ' + str(i) + ' of ' + str(len(gene_expression_lists)))

            scores_lists.append(scores_list)

        return scores_lists


def get_clusters(names, paths):
    if len(names) != len(paths):
        raise ValueError('Number of cluster names must equal length of cluster paths')

    clusters = {}

    for i, name in names:
        clusters[name] = paths[i]

    return clusters


if __name__ == '__main__':

    # Parse command-line arguments
    ap = ArgumentParser(description=__doc__,  # Use text from file summary up top
                        formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument('infercnv_output',
                    help='Path to pre_vis_transform.txt output from inferCNV')
    ap.add_argument('gen_pos_path',
                    help='Path to gen_pos.txt genomic positions file from inferCNV ',
                    nargs='+')
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
    gen_pos_path = args.gen_pos_path
    cluster_names = args.cluster_names
    cluster_paths = args.cluster_paths
    output_path = args.output_path

    clusters = get_clusters(cluster_names, cluster_paths)

    MatrixToIdeogramAnnots(infercnv_output, gen_pos_path, clusters, output_path)
