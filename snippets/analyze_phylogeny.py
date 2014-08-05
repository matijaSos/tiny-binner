import sys, os
sys.path.append(os.getcwd())
from itertools import chain
from collections import defaultdict

import numpy as np
from pylab import *
from matplotlib import cm

import ncbi.taxonomy.tree as taxtree
import ncbi.taxonomy.ranks as rank

def iter_tax_nodes(fname, tax_id_index):
	with open(fname) as fin:
		for line in fin:
			yield int(line.strip().split()[tax_id_index])

def iter_filtered_taxa(tree, rank2num, nodes, level):
	for node in nodes:
		rank = tree.nodes[node].rank
		if rank == 'None':
			continue
		ranknum = rank2num[rank]
		if ranknum <= level or ranknum == 28:
			yield node

def plot_data(composition, tree, _title=None):
	N = len(composition)
	colors = cm.Set1(np.arange(N, 0, -1)/float(N))
	total = sum(composition.values())

	# make a square figure and axes
	#figure(1, figsize=(6,6))
	figure(figsize=(6,6))
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
	if _title is None:
		title('Bacterial class distribution')
	else:
		title(_title)
	show()

def determine_composition(tax_nodes, tree, depth):
	rank2num = rank.ranks
	composition = defaultdict(int)
	interesting_taxa = filter(lambda n: rank2num.get(tree.nodes[n].rank, -1) == depth, tree.nodes.keys())
	for taxid in interesting_taxa:
		for found_node in tax_nodes:
			if tree.is_child(found_node, taxid):
				composition[taxid] += 1
	return composition


def main():
	rank2num = rank.ranks
	if len(sys.argv) < 4:
		print 'Usage: python analyze_phylogeny.py [INPUT_FILE] [TAXONOMIC_DEPTH] [TAX_ID_INDEX] [[OUTPUT_FIGURE]]'
		print 'Available taxonomic depths:'
		print '\n'.join(map(lambda r: '%3d\t%s' % (r[1], r[0]), sorted(rank2num.items(), key=lambda i: i[1])))
		sys.exit(-1)

	depth = int(sys.argv[2])
	tax_id_index = int(sys.argv[3])
	print 'Loading NCBI taxonomy tree...'
	tree = taxtree.TaxTree()
	print 'DONE'

	tax_nodes = list(iter_tax_nodes(sys.argv[1], tax_id_index))

	composition = determine_composition(tax_nodes, tree, depth)

	plot_data(composition, tree)


if __name__ == '__main__':
	main()
