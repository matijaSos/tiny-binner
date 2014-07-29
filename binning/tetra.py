import math

from itertools import product
from collections import defaultdict
from contextlib import nested

import numpy as np
from Bio import SeqIO

from ncbi.taxonomy.tree import TaxTree


VIRTUAL_ZERO = 1e-10

TETRAS = sorted(map(lambda l: ''.join(l), product('ATGC', repeat=4)))
TRIPLETS = sorted(map(lambda l: ''.join(l), product('ATGC', repeat=3)))

def get_tetranucl_z_scores(seq):
	freqs_quad = defaultdict(int)
	freqs_tri = defaultdict(int)
	freqs_duo = defaultdict(int)
	for i in range(0, len(seq) - 3):
		freqs_quad[seq[i:i+4]] += 1
	for i in range(0, len(seq) - 2):
		freqs_tri[seq[i:i+3]] += 1
	for i in range(0, len(seq) - 1):
		freqs_duo[seq[i:i+2]] += 1
	N = len(seq)
	z_scores = []
	for tetra in TETRAS:
		trinucl1 = tetra[:-1]
		trinucl2 = tetra[1:]
		dinucl = tetra[1:-1]
		exp = float(freqs_tri[trinucl1] * freqs_tri[trinucl2]) / max(freqs_duo[dinucl], VIRTUAL_ZERO)
		var = exp * (freqs_duo[dinucl] - freqs_tri[trinucl1]) * \
					(freqs_duo[dinucl] - freqs_tri[trinucl2]) / max(freqs_duo[dinucl]**2, VIRTUAL_ZERO)
		z_scores.append((freqs_quad[tetra] - exp) / max(math.sqrt(var), VIRTUAL_ZERO))
	return np.array(z_scores)

def pearson_corr(zscores1, zscores2):
	return sum(zscores1 * zscores2) / 256.

def assign_sequence(seq, db_freqs):
	seq_z_scores = get_tetranucl_z_scores(seq)
	corrs = {}
	for tax, z_score in db_freqs.iteritems():
		corrs[pearson_corr(seq_z_scores, z_score)] = tax
	max_corr = max(corrs.keys())
	print min(corrs.keys())
	return max_corr, corrs[max_corr]


def load_database(freq_db, gi_key=False):
	freqs = {}
	with open(freq_db) as fin:
		while True:
			line1 = fin.readline()
			if not line1:
				break
			if line1.startswith('#'):
				continue
			line1 = line1.strip()
			gi, tax = line1.split()
			gi, tax = int(gi), int(tax)
			line2 = fin.readline().strip()
			freq = np.array(map(lambda f: float(f), line2.split()))
			if gi_key:
				freqs[gi] = freq
			else:
				freqs[tax] = freq
	return freqs


def create_database(fasta_db, output_file, data_access):
	frequency_db = {}

	with nested(open(fasta_db), open(output_file, 'w')) as (fin, fout):

		fout.write('# SOURCE FILE: %s\n# FORMAT:\n' % fasta_db)
		fout.write('# gi taxid\n')
		fout.write('# frequencies (256 element vector from AAAA to TTTT)\n')

		i = 0
		records = SeqIO.parse(fin, 'fasta')
		for rec in records:
			seq = str(rec.seq)
			gi = int(rec.name.split('|')[1])
			tax = data_access.get_taxids([gi])
			if gi not in tax:
				continue
			else:
				tax = tax[gi]
			freqs = get_tetranucl_z_scores(seq)
			fout.write('%d %d\n' % (gi, tax))
			fout.write('%s\n' % ' '.join(map(lambda f: '%.4f' % f, freqs)))
			i += 1
			print 'Processing %5d. sequence...' % i
