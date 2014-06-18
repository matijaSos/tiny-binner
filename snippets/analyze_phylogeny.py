import sys, os
from collections import defaultdict
from itertools import chain
from pylab import *
import numpy as np
from matplotlib import cm
sys.path.append(os.getcwd())

import ncbi.taxonomy.tree as taxtree
import ncbi.taxonomy.ranks as rank

def iter_tax_nodes(fname):
	with open(fname) as fin:
		for line in fin:
			yield int(line.strip().split('\t')[1])

def iter_filtered_taxa(tree, rank2num, nodes, level):
	for node in nodes:
		rank = tree.nodes[node].rank
		if rank == 'None':
			continue
		ranknum = rank2num[rank]
		if ranknum <= level or ranknum == 28:
			yield node

def plot_data(composition, tree):
	N = len(composition)
	colors = cm.Set1(np.arange(N, 0, -1)/float(N))
	total = sum(composition.values())

	# make a square figure and axes
	figure(1, figsize=(6,6))
	ax = axes([0.1, 0.1, 0.8, 0.8])

	# The slices will be ordered and plotted counter-clockwise.
	#labels = map(lambda p: tree.nodes[p].organism_name, composition.keys())
	fracs = map(lambda p: float(p) / total, composition.values())
	labels = []
	for tid, num in composition.items():
		org_name = tree.nodes[tid].organism_name
		frac = float(num) / total * 100
		labels.append('%s (%2.2f %%)' % (org_name, frac))

	pie(fracs, colors=colors, startangle=90)
	legend(labels, loc="lower right")
	title('Bacterial class distribution')
	show()


def main():
	rank2num = rank.ranks
	if len(sys.argv) < 3:
		print 'Usage: python analyze_phylogeny.py [INPUT_FILE] [TAXONOMIC_DEPTH] [[OUTPUT_FIGURE]]'
		print 'Available taxonomic depths:'
		print '\n'.join(map(lambda r: '%3d\t%s' % (r[1], r[0]), sorted(rank2num.items(), key=lambda i: i[1])))
		sys.exit(-1)

	depth = int(sys.argv[2])
	print 'Loading NCBI taxonomy tree...'
	tree = taxtree.TaxTree()
	print 'DONE'

	tax_nodes = list(iter_tax_nodes(sys.argv[1]))

	composition = defaultdict(int)
	interesting_taxa = filter(lambda n: rank2num.get(tree.nodes[n].rank, -1) == depth, tree.nodes.keys())
	for taxid in interesting_taxa:
		for found_node in tax_nodes:
			if tree.is_child(found_node, taxid):
				composition[taxid] += 1

	plot_data(composition, tree)


if __name__ == '__main__':
	main()
