import guppy.heapy.RM
import time
import logging
import logging.config
import argparse
import sys,os
import operator
#sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer

import filters.readprocessing as rstate
from filters.binning import bin_reads
from utils import timeit
from utils.location import Location
from formats.xml_output import *

from testing.solution import Solution   
from testing.rank_accuracy import RankAccuracy

from binning.LCABinner import LCABinner

class TestRunArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(TestRunArgParser, self).__init__('''\
                Loads alignment file, \
                assigns tax id to each alignment and \
                using LCA algorithm assigns each read to a tax node.''')
        self.add_argument('alignment_file',
                help='path to the read alignment file',
                type=str)
def main():
    '''
    Script to perform LCA binning with a given read alignment file.
    '''

    # Input arguments
    argparser = TestRunArgParser()
    args  = argparser.parse_args()

    # Access database
    dataAccess = DataAccess(args)

    print '1. Loading tax tree...'
    tax_tree = TaxTree()
    print 'done.'

    print '2. Loading alignment file...'
    read_container = ReadContainer()
    read_container.load_alignment_data(args.alignment_file)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)
    print 'done'

    print '4. Creating LCA solution...'
    lca_binner = LCABinner(tax_tree)
    lca_sol = lca_binner.bin_reads(read_container)
    print 'done'

    print '5. Nice output of read assignment here'
    lca_sol.print_nicely(tax_tree)

if __name__ == '__main__':
    main()
