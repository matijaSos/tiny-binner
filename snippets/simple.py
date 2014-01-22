import guppy.heapy.RM
import time
import logging
import logging.config
import argparse
import sys,os
import operator
sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer
from data.containers.record import RecordContainer
from data.containers.cdsaln import CdsAlnContainer
import filters.host as host_filter
import filters.readprocessing as rstate
from filters.binning import bin_reads
from utils import timeit
from utils.location import Location
from formats.xml_output import *

from testing.solution import Solution   
from testing.rank_accuracy import RankAccuracy

class TestRunArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(TestRunArgParser, self).__init__('''\
                Loads input alignment file, fills read container, \
                filters host reads and fills record container \
                with read from specified database table.''')
        self.add_argument('input',
                help='input alignment file',
                type=str)
        self.add_argument('metasim_fasta',
                help='FASTA file generated by metasim - used to evaluate solution',
                type=str)
def main():
    '''
    Script to experiment with binner.
    '''
    print "Hello world!"

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
    read_container.load_alignment_data(args.input)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)
    print 'done'

    '''
    # Output species and number of alignments to each
    speciesCount = {}
    for read in read_container.fetch_all_reads():
        for aln in read.get_alignments():
            # print "\t" + str(aln.tax_id)
            # print "\t" + dataAccess.get_organism_name(aln.tax_id)
            if aln.tax_id in speciesCount:
                speciesCount[aln.tax_id] += 1
            else:
                speciesCount[aln.tax_id] = 1

    # Sort by occs
    sorted_speciesCount = sorted(speciesCount.iteritems(), key=lambda x: x[1], 
                                 reverse = True)
    for tax_id, count in sorted_speciesCount:
        if count < 150:
            continue
        print dataAccess.get_organism_name(tax_id) + ": \t" + str(count)
    '''

    print '3. Loading correct solution...'
    sol = Solution(args.metasim_fasta, dataAccess)
    print 'done'
    
    print '4. Evaluating our solution...'
    rankAcc = RankAccuracy(tax_tree, sol, read_container)
    print 'done'

    # Print test results
    rankAcc.printData()
    


if __name__ == '__main__':
    main()
