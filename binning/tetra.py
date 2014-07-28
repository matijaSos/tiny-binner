import sys,os
sys.path.append(os.getcwd())

from itertools import product
from collections import defaultdict
from contextlib import nested

from Bio import SeqIO

from ncbi.db.data_access import DataAccess
from utils.argparser import DefaultBinnerArgParser


TETRAS = sorted(map(lambda l: ''.join(l), product('ATGC', repeat=4)))

def analyze_sequence(seq):
	frequencies = defaultdict(int)
	for i in range(0, len(seq) - 3):
		frequencies[seq[i:i+4]] += 1
	N = len(seq)
	Ef = 1./256
	for t in TETRAS:
		frequencies[t] = (float(frequencies[t])/N - Ef) * 100
	return frequencies

def create_database(fasta_db, output_file, data_access):
	frequency_db = {}

	with nested(open(fasta_db), open(output_file, 'w')) as (fin, fout):

		fout.write('# SOURCE FILE: %s\n# FORMAT:\n' % fasta_db)
		fout.write('# gi taxid\n')
		fout.write('# frequencies (256 element vector from AAAA to TTTT)')

		i = 0
		records = SeqIO.parse(fin, 'fasta')
		for rec in records:
			seq = str(rec.seq)
			gi = int(rec.name.split('|')[1])
			tax = data_access.get_taxids([gi])[gi]
			freqs = analyze_sequence(seq)
			fout.write('%d %d\n' % (gi, tax))
			fout.write('%s\n', ' '.join(map(lambda f: '%.4f' % f, [freqs[t] for t in TETRAS])))
			i += 1
			print i


def main():

	argparser = DefaultBinnerArgParser('Creates a database from input fasta file.')
	args = argparser.parse_args()
	data_access = DataAccess(args)
	fasta_db = '/home/ana/Data/metagenome/db/all_bacteria.fa'
	output_file = '/home/ana/Data/metagenome/db/allbact.tetra.db'
	create_database(fasta_db, output_file, data_access)

if __name__ == '__main__':
	main()