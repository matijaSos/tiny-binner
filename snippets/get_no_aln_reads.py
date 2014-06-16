
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
                takes all reads assigned with zero alignments \
                and exports their ids to the given file.
                .''')
        self.add_argument('alignment_file',
                help='path to the read alignment file',
                type=str)
        self.add_argument('read_ids_out',
                help='path to the file where ids of reads with no\
                      alignments will be exported',
                type=str)


def main():
    '''
    Script to extract reads without any alignments
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

    # Loop through reads and take those with no alignments
    out_file = open(args.read_ids_out, 'w')
    
    no_aln_count = 0
    for read in read_container.fetch_all_reads(format=list):
        if not read.has_alignments():
            out_file.write("{0}\n".format(read.id))
            no_aln_count += 1

    out_file.close()
    
    total_read_count = read_container.get_read_count()

    print("total number of reads   : {0}".format( total_read_count ))
    print("reads without alignments: {0}".format(no_aln_count))
    print
    print("no aln percentage: {0:.2f}%".format(no_aln_count * 100 / float(total_read_count)))

if __name__ == '__main__':
    main()
