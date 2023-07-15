from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqUtils import GC
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import os


class FromGenes2Networks:

    def __init__(self, genebank_file_path, output_dir):
        self.genebank_path = genebank_file_path
        self.output_dir = output_dir

    def run(self):
        """
        1. pars CDSs sequences from gene bank data file
        2. calculate alignment matrix
        3. rank genes by PageRank, %GC and sequence length
        4. plot the dilute matrix as network
        5. plot correlation graphs for each pair of ranks
        """

        # pars data
        (seq_lst, name_lst, seq_sz_lst, gc_lst) = self.get_seq_data()
        lst_sz = len(seq_lst)
        mean_gc = sum(gc_lst) / lst_sz
        gc_lst = [abs(x - mean_gc) for x in gc_lst]

        # create alignment matrix
        alignment_mat = self.get_alignment_mat(seq_lst)

        # Analyze alignment matrix
        matrix_graph = nx.DiGraph(alignment_mat)
        pageRank = nx.pagerank(matrix_graph, alpha=0.9)
        pageRank = list(pageRank.values())

        # keep the features we want to explore
        d = {'seq_sz_lst': seq_sz_lst, 'pageRank': pageRank, 'gc_lst': gc_lst}
        df = pd.DataFrame(data=d, index=list(range(0, lst_sz)))

        self.plot_network(alignment_mat, name_lst)

        self.plot_correlation(df)

    def get_seq_data(self):
        """
        :return: function goal is to parse the input geneBank data file
        """

        seq_lst = list()
        name_lst = list()
        seq_sz_lst = list()
        gc_lst = list()

        # find the CDS records in the geneBank file
        for rec in SeqIO.parse(self.genebank_path, "genbank"):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        # store CDS data to lists
                        tmp_seq = feature.location.extract(rec).seq
                        seq_lst.append(str(tmp_seq)[:150])
                        name_lst.append(feature.qualifiers["protein_id"])
                        seq_sz_lst.append(len(str(tmp_seq)))
                        gc_lst.append(GC(tmp_seq))
        return seq_lst, name_lst, seq_sz_lst, gc_lst

    @staticmethod
    def get_alignment_mat(seq_lst):
        """
        :param seq_lst: pairs for genes (nodes)
        :return: function goal is to calculate alignment between pairs of genes and return the alignment matrix.
        """

        print('Please wait while I am calculating alignment Matrix...')

        genes_sz = len(seq_lst)
        alignment_mat = np.zeros((genes_sz, genes_sz))

        # calculate pairwise alignment between genes pairs
        for i in range(0, genes_sz - 1):
            for j in range(i + 1, genes_sz):
                tmp_score = pairwise2.align.globalxx(seq_lst[i], seq_lst[j], score_only=True)
                alignment_mat[i, j] = tmp_score

        # add the matrix conjugate transpose for Symmetry
        alignment_mat = alignment_mat.conj().T + alignment_mat
        alignment_mat = alignment_mat / alignment_mat.max()

        return alignment_mat

    def plot_network(self, alignment_mat, name_lst, percentile_threshold=99.9):
        """

        :param alignment_mat:
        :param name_lst:
        :param percentile_threshold: by default leave Highest 0.1%
        :return: display matrix as a network graph after diluting it by cut-off
        """

        # use original matrix for spring embedded layout
        matrix_graph = nx.DiGraph(alignment_mat)

        # create the network
        top_001 = np.percentile(alignment_mat, percentile_threshold)
        top_001_mat = (alignment_mat >= top_001) * alignment_mat
        top_001_graph = nx.DiGraph(top_001_mat)

        # add nodes labels
        label_dict = {}
        idx = 0
        for tmp_nm in name_lst:
            label_dict[idx] = tmp_nm
            idx += 1

        # plot the graph
        plt.figure()
        plt.title('Genes Alignment Network')
        nx.draw_networkx(top_001_graph, font_size=6, node_size=30, labels=label_dict, with_labels=True,
                         pos=nx.spring_layout(matrix_graph))
        plt.savefig(os.path.join(self.output_dir, 'GenesAlignmentNetwork'))
        plt.show()

    @staticmethod
    def plot_correlation(df):
        """
        :param df: the network df
        :return: plot for each pair of data sets the correlation graphs with gergesor and calculate the R & p-value.
        """

        g = sns.jointplot(y="seq_sz_lst", x="pageRank", kind="reg", data=df)
        g.fig.suptitle('Correlation between PageRank & Sequence Size')
        g = sns.jointplot(y="gc_lst", x="pageRank", kind="reg", data=df)
        g.fig.suptitle('Correlation between PageRank & %GC distance')
        g = sns.jointplot(x="seq_sz_lst", y="gc_lst", kind="reg", data=df)
        g.fig.suptitle('Correlation between %GC distance & Sequence Size')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='for a given set of genes, plot a visual '
                                                 'network based on different features')
    parser.add_argument('genebank_file_path', type=str, help='a file with a set of genes in a GeneBank format (.gb)')
    parser.add_argument('output_dir', type=str, help='a path to the output results directory')

    args = parser.parse_args()

    instance = FromGenes2Networks(args.genebank_file_path, args.output_dir)
    instance.run()
