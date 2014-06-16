import sys, os
from collections import defaultdict
from itertools import chain
import pdb
from pylab import *
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

def plot_data(phylo_levels):
	interesting_data = filter(lambda p: p[1] > 0, phylo_levels[7])
	total = sum(map(lambda p: p[1], interesting_data))

	# make a square figure and axes
	figure(1, figsize=(6,6))
	ax = axes([0.1, 0.1, 0.8, 0.8])

	# The slices will be ordered and plotted counter-clockwise.
	labels = map(lambda p: p[0], interesting_data)
	fracs = map(lambda p: float(p[1]) / total, interesting_data)
	explode=(0, 0.05, 0, 0)

	pie(fracs, labels=labels,
	                autopct='%1.1f%%', shadow=True, startangle=90)
	                # The default startangle is 0, which would start
	                # the Frogs slice on the x-axis.  With startangle=90,
	                # everything is rotated counter-clockwise by 90 degrees,
	                # so the plotting starts on the positive y-axis.

	title('Bacterial class distribution', bbox={'facecolor':'0.8', 'pad':5})

	interesting_data = filter(lambda p: p[1] > 0, phylo_levels[4])
	# make a square figure and axes
	figure(2, figsize=(6,6))
	ax = axes([0.1, 0.1, 0.8, 0.8])

	# The slices will be ordered and plotted counter-clockwise.
	labels = map(lambda p: p[0], interesting_data)
	fracs = map(lambda p: float(p[1]) / total, interesting_data)
	explode=(0, 0.05, 0, 0)

	pie(fracs, labels=labels,
	                autopct='%1.1f%%', shadow=True, startangle=90)
	                # The default startangle is 0, which would start
	                # the Frogs slice on the x-axis.  With startangle=90,
	                # everything is rotated counter-clockwise by 90 degrees,
	                # so the plotting starts on the positive y-axis.

	title('Phylum distribution', bbox={'facecolor':'0.8', 'pad':5})
	show()


def main():
	if len(sys.argv) < 3:
		print 'Usage: python analyze_phylogeny.py [INPUT_FILE] [ANALYSIS_DEPTH]'
		sys.exit(-1)

	rank2num = rank.ranks
	depth = int(sys.argv[2])
	print 'Loading NCBI taxonomy tree...'
	tree = taxtree.TaxTree()
	print 'DONE'

	tax_nodes = list(iter_tax_nodes(sys.argv[1]))

	level2taxid = defaultdict(list)
	interesting_taxa = set()
	level2taxid[0] = [1]
	for level in range(1, depth + 1):
		parents = level2taxid[level - 1]
		children = list(chain(map(lambda tid: tree.child_nodes[tid], parents)))
		children = [y for x in children for y in x]
		filt_children = list(iter_filtered_taxa(tree, rank2num, children, depth))
		level2taxid[level].extend(filt_children)
		interesting_taxa.update(filt_children)

	del level2taxid

	composition = defaultdict(int)
	for node in interesting_taxa:
		for tax in tax_nodes:
			if tree.is_child(tax, node):
				composition[node] += 1


	phylo_levels = defaultdict(list)
	for level in range(0, depth + 1):
		lvl_taxa = filter(lambda tid: rank2num[tree.nodes[tid].rank] == level, interesting_taxa)
		if not lvl_taxa:
			continue
		print tree.nodes[lvl_taxa[0]].rank
		for tid in lvl_taxa:
			phylo_levels[level].append((tree.nodes[tid].organism_name, composition[tid]))
			if composition[tid] > 0:
				print '\t(%d)' % tid, tree.nodes[tid].organism_name, composition[tid]

	plot_data(phylo_levels)


if __name__ == '__main__':
	main()
