import argparse 
import random
import string
import sys

LARGE_COVERAGE = 2
SMALL_COVERAGE = 1
SEQ_LEN        = 100
GI             = 188558146 # E.coli genome index
DB_SOURCE      = 'tst'

def generate_random_sequence (seq_len):
    return ''.join(random.choice('ATGC') for x in range(seq_len))

class Coverage (object):
    def __init__ (self, coverage, seq_len):
        self.sequence = generate_random_sequence(seq_len)
        self.coverage = coverage
        self.reads = self._generate_reads()

    def _generate_reads(self):
        reads = []
        for i in range(self.coverage):
            seq_id = ''.join(random.choice(string.ascii_uppercase+string.digits) for x in range(20))
            sequence = self.sequence
            reads.append ((seq_id, sequence))
        return reads

class CDS (object):
    def __init__ (self, coverage_list, active_coverage_list):
        self.id = ''.join(random.choice(string.ascii_uppercase+string.digits) for x in range(20)) + '.1'
        self.header = "gi|%d|%s|%s|" % (GI, DB_SOURCE, self.id)
        self.sequence = ''
        for cov in coverage_list:
            self.sequence += cov.sequence + '\n'

        self.active_reads = []
        for cov in active_coverage_list:
            for (read_id, seq) in cov.reads:
                self.active_reads.append(read_id)


argparser = argparse.ArgumentParser(description='Uses parameters to work')

argparser.add_argument('cds_number', help='number of CDSs to be generated', type=int)
argparser.add_argument('read_output_file', help='path to read output fasta file', default='reads.fa', type=str, nargs='?')
argparser.add_argument('cds_output_file', help='path to CDS output fasta file', default='cds.fa', type=str, nargs='?')
argparser.add_argument('cds_ordinal_output_file', help='file to write the correct order of CDS choices', default='cds_ordering.txt', type=str, nargs='?')

args=argparser.parse_args()

init_cds_segments = args.cds_number
large_coverage = []
small_coverage = []
cds_ordinals   = {}

# number of CDSs
N = args.cds_number

# generate large coverage
num_large = (N-1)*N/2
for i in range(0, num_large):
    large_coverage.append(Coverage(LARGE_COVERAGE, SEQ_LEN))

# generate small_coverage
num_small = N
for i in range(0, num_small):
    small_coverage.append(Coverage(SMALL_COVERAGE, SEQ_LEN))


cds_ordinals = {}
large_cov_index = 0
active_reads = 0
for x in range (0, N):
    coverages = []
    for y in range (0, N):
        if (x+y < N-1):
            coverages.append (large_coverage[large_cov_index])
            large_cov_index += 1
        else:
            coverages.append (small_coverage[y])
    cds_ordinals[x+1] = CDS(coverages, coverages[0:N-active_reads])
    active_reads += 1

# write CDSs to a file
cds_output = open(args.cds_output_file, 'w')
cds_ordinal_output = open(args.cds_ordinal_output_file, 'w')
for (ordinal, cds) in cds_ordinals.items():
    cds_output.        write(">%s\n%s\n" % (cds.header, cds.sequence))
    cds_ordinal_output.write("%s\n"  % cds.id)
    cds_ordinal_output.write("%s\n" % ";".join(cds.active_reads))
cds_output.close()
cds_ordinal_output.close()

# write reads to a file
coverages = large_coverage + small_coverage
reads_output = open(args.read_output_file, 'w')
for cov in coverages:
    reads = cov.reads
    for (read_id, sequence) in reads:
        reads_output.write(">%s\n%s\n" % (read_id, sequence))
reads_output.close()

