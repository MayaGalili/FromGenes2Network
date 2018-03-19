# -*- coding: utf-8 -*-
"""
Created on Friday Jun 23 17:03:46 2017

@author: Maya Galili
"""
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqUtils import GC
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

         
def get_seq_data(geneBank_file):
#########################################
# 	function goal is to pars the input  
#   geneBank data file.                 
#########################################

    # init vars
    seq_lst = []
    name_lst = []
    seq_sz_lst = []
    gc_lst = []
    
	# find the CDS records in the geneBank file
    for rec in SeqIO.parse(geneBank_file, "genbank"):
        if rec.features:
            for feature in rec.features:
                if feature.type == "CDS":
					# store CDS data to lists
                    tmp_seq = (feature.location.extract(rec).seq)
                    seq_lst.append(str(tmp_seq)[:150])
                    name_lst.append(feature.qualifiers["protein_id"])
                    seq_sz_lst.append(len(str(tmp_seq)))
                    gc_lst.append(GC(tmp_seq))
    return (seq_lst, name_lst, seq_sz_lst, gc_lst)


def get_alignmnt_mat(seq_lst):
##############################################
# 	function goal is to calculate alignment 
#  	between paires of genes and return the  
#  	alignment matrix.                       
##############################################

    print('Please wait while I am calculating alignment Matrix.' )
    # init vars
    genes_sz = len(seq_lst)
    alignmnt_mat = np.zeros((genes_sz,genes_sz))
    
    # calculate pairwise alignment between genes paires 
    for i in range(0,genes_sz-1):
        for j in range(i+1, genes_sz):
            tmp_score=pairwise2.align.globalxx(seq_lst[i], seq_lst[j], score_only=True)
            alignmnt_mat[i, j] = tmp_score
    
    # add the matrix conjugate transpose for Symmetry
    alignmnt_mat = alignmnt_mat.conj().T+alignmnt_mat
    alignmnt_mat = alignmnt_mat/alignmnt_mat.max()

    return alignmnt_mat  	


def plot_network(alignmnt_mat,name_lst):
##############################################
# 	function goal is to display matrix as a 
#   network graph after diluting it by     
#   cut-off (Highest 0.1%)                 
##############################################
	
	# use original matrix for spring embedded layout
	matrix_graph =  nx.from_numpy_matrix(alignmnt_mat)

	# create the network
	top_001 = np.percentile(alignmnt_mat, 99.9)
	top_001_mat = (alignmnt_mat>=top_001)*alignmnt_mat
	top_001_graph =  nx.from_numpy_matrix(top_001_mat)

	# add nodes labels
	labeldict = {}
	idx = 0
	for tmp_nm in name_lst:
		labeldict[idx] = tmp_nm
		idx+=1

	# plot the graph
	plt.figure()
	plt.title('Genes Alignment Network')
	nx.draw_networkx(top_001_graph,font_size=6,node_size=30,labels=labeldict,with_labels = True,pos=nx.spring_layout(matrix_graph))
	plt.savefig('GenesAlignmentNetwork') # save as png
	plt.show()	
	
def plot_corelation(df):
##############################################
# 	function goal is to plot for each paire 
# 	of data sets the correlation graphs with 
# 	gergesor and calculate the R & p-value. 
##############################################

	g=sns.jointplot(y="seq_sz_lst", x="pageRank",kind="reg",data=df)
	g.fig.suptitle('Corelation between PageRank & Sequence Size')
	g=sns.jointplot(y="gc_lst", x="pageRank",kind="reg",data=df)
	g.fig.suptitle('Corelation between PageRank & %GC distance')
	g=sns.jointplot(x="seq_sz_lst", y="gc_lst",kind="reg",data=df)
	g.fig.suptitle('Corelation between %GC distance & Sequence Size') 
	
#############################################
#	main script - script will:               
#   1. pars CDSs sequences from gene bank   
#	   data file                            
#   2. calculate alignment matrix          			
#   3. rank genes by PageRank, %GC and     
#      sequence length                     
#   4. plot the dilute matrix as network   
#   5. plot corelation graphs for each pair
#       of ranks                           
#############################################

# input file
geneBank_file = 'sequence.gb'

# pars data
(seq_lst, name_lst, seq_sz_lst, gc_lst) = get_seq_data(geneBank_file)
lst_sz = len(seq_lst)
mean_gc = sum(gc_lst)/lst_sz
gc_lst = [abs(x-mean_gc) for x in gc_lst]

# create alignment matrix 
alignmnt_mat = get_alignmnt_mat(seq_lst)

# Analys alignment matrix
matrix_graph =  nx.from_numpy_matrix(alignmnt_mat)
pageRank = nx.pagerank(matrix_graph, alpha=0.9)
pageRank = list(pageRank.values())

# keep the features we want to explore
d = {'seq_sz_lst': seq_sz_lst, 'pageRank': pageRank,'gc_lst':gc_lst}
df = pd.DataFrame(data=d, index=list(range(0,lst_sz)))
             
# plot the alignmnt_mat in spring-embedded layout
plot_network(alignmnt_mat,name_lst)

# plot the corelation graphs
plot_corelation(df)    