import guppy.heapy.RM
import time
import logging
import logging.config
import argparse
import sys,os
sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer
from data.containers.record import RecordContainer
from data.containers.cdsaln import CdsAlnContainer
import filters.host as host_filter
from filters.bin import bin as bin_reads
from utils import timeit
from utils.location import Location

class TestRunArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(TestRunArgParser, self).__init__('''\
                Loads input alignment file, fills read container, \
                filters host reads and fills record container \
                with read from specified database table.''')
        self.add_argument('input',
                help='input alignment file',
                type=str)
        self.add_argument('output',
                help='output file with rna information',
                type=str)



def log_start (log, args):
    log.info('BINNER RUN')
    log.info("Input: %s" % args.input)



def main():
    '''
    Script to run binner in one of the most common
    usage scenarios.
    * load alignment data
    * load taxonomy data
    * do basic alignment data filtering (remove host reads ecc)
    '''

    #----------------------------------#
    #------ INPUT ARGUMENTS -----------#
    argparser = TestRunArgParser()
    args  = argparser.parse_args()

    #----------------------------------#
    #------- STATIC DATA SOURCE -------#
    # CDS - GI2TAXID -- NAMES -- NODES #
    dataAccess = DataAccess(args)
    #raw_input('Data access created')
    #----------------------------------#

    #-------- TAXONOMY TREE -----------#
    print '1. Loading tax tree...'
    tax_tree = TaxTree()
    # tax_tree.load_taxonomy_data(dataAccess)
    print 'done.'

    #----------------------------------#
    #------- ALIGNMENT DATA SOURCE ----#
    print '2. Loading alignment file...'
    read_container = ReadContainer()
    read_container.load_alignment_data(args.input)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)
    print 'done'

    #------- FILTER HOST READS -------#
    print '3. Filtering host reads & alignments...'
    new_reads = host_filter.filter_potential_host_reads(
        read_container.fetch_all_reads(format=list),
        tax_tree.tax2relevantTax,
        tax_tree.potential_hosts,
        #delete_host_alignments =
        True,
        #filter_unassigned =
        True,
        #unassigned_taxid=
        -1,
        host_filter.is_best_score_host)
    dataAccess.clear_cache()    # deletes gi2taxid cache
    reads_with_no_host_alignments = host_filter.filter_potential_hosts_alignments(
        new_reads,
        tax_tree.tax2relevantTax,
        tax_tree.potential_hosts,
        True,   # delete host alignments
        True,   # filter unassigned
        -1)     # unassigned taxid
    read_container.set_new_reads(reads_with_no_host_alignments)
    print 'done'

    #----------------------------------#
    #------- LOAD ALL RECORDS   -------#
    print '4. Loading referenced records...'
    record_container = RecordContainer()
    record_container.set_db_access(dataAccess)
    record_container.populate(read_container.fetch_all_reads_versions(), table='cds')
    print 'done'
    #----------------------------------#
    #-- MAP ALIGNMENTS TO GENES   -----#
    print '5. Mapping alignments to genes...'
    read_container.populate_cdss(record_container)
    #----------------------------------#
    #- RECORD ALL ALIGNEMENTS TO GENE -#
    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))
    print 'done'

    print '6. Estimating organisms present in sample...'
    target_taxids = [633, 632, 263, 543, 86661, 1392, 55080, 1386]
    print 'done.'
   
    print '7. Binning reads...' 
    bin_reads(read_container.read_repository,
              cds_aln_container.cds_repository,
              cds_aln_container.read2cds,
              tax_tree,
              target_taxids,
              None,
              None,
              True) 
    print 'done.'

if __name__ == '__main__':
    main()
