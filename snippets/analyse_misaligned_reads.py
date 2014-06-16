
import argparse
import sys,os
import time
sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer

import filters.host as host_filter

#  For CDS loading
from data.containers.record import RecordContainer
from data.containers.cdsaln import CdsAlnContainer

from utils.location import Location

from binning.LCABinner import LCABinner

class ArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(ArgParser, self).__init__('''\
                Loads alignment file, \
                takes all given reads with \
                and determines species they are aligning to.
                .''')
        self.add_argument('alignment_file',
                help='path to the read alignment file',
                type=str)
        self.add_argument('read_ids_file',
                help='path to the file containing read ids',
                type=str)
        self.add_argument('out_file',
                help='output file',
                type=str)

def main():
    '''
    Script to analyse given reads and determine species they are
    aligning to.
    '''

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()

    # Access database
    dataAccess = DataAccess(args)

    # ------------------ #

    print '1. Loading tax tree...'
    start = time.time()

    tax_tree = TaxTree()

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    print '2. Loading alignment file...'
    start = time.time()

    read_container = ReadContainer()
    read_container.load_alignment_data(args.alignment_file)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    print '3. Loading read_ids file...'
    start = time.time()

    # Get file contents
    read_ids = [line.rstrip('\n') for line in open(args.read_ids_file)]
    
    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    print '4. Going over given reads...'
    start = time.time()

    tax_ids = {}
    total_alns = 0
    for read_id in read_ids:
        read = read_container.fetch_read(read_id)

        # Go through alignments
        for aln in read.get_alignments():
            tax_ids[aln.tax_id] = tax_ids.get(aln.tax_id, 0) + 1

            total_alns += 1

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # -------- Stats output ---------- #
    

    out_file = open(args.out_file, 'w')
    out_file.write("total_alns: {0}\n\n".format(total_alns))

    for tax_id in sorted(tax_ids, key=tax_ids.get, reverse=True):
        num     = tax_ids[tax_id]
        frac    = float(num)/total_alns * 100
        try:
            name = tax_tree.nodes[int(tax_id)].organism_name
        except:
            name = "unknown"
        out_file.write("{0:10} {1:80} {2:8d} {3:10.4f}%\n".format(tax_id, name, num, frac))

    out_file.close()

if __name__ == '__main__':
    main()
