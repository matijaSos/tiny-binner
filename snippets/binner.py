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
import filters.host as host_filter
from utils import timeit

class TestRunArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(TestRunArgParser, self).__init__('''\
                Loads input alignment file, fills read container, \
                filters host reads and fills record container \
                with read from specified database table.''')
        self.add_argument('input',
                help='input alignment file',
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
    #------- DATA LOGGING INIT --------#
    log = logging.getLogger(__name__)
    logging.config.fileConfig(args.log_configuration,
                              disable_existing_loggers=False)
    log_start(log, args)


    #----------------------------------#
    #------- STATIC DATA SOURCE -------#
    # CDS - GI2TAXID -- NAMES -- NODES #
    dataAccess = DataAccess(args)
    raw_input('Data access created')
    #----------------------------------#

    #-------- TAXONOMY TREE -----------#
    tax_tree = TaxTree()
    # tax_tree.load_taxonomy_data(dataAccess)
    raw_input('Tax tree loaded.')

    #----------------------------------#
    #------- ALIGNMENT DATA SOURCE ----#
    read_container = ReadContainer()
    timeit(read_container.load_alignment_data, args.input)
    raw_input('Alignment file loaded')
    print 'Number of PREFILTERED reads: ', len(read_container.fetch_all_reads(format=list))
    timeit(read_container.set_taxids, dataAccess)
    raw_input('Data access set')

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

    read_container.set_new_reads(new_reads)
    raw_input('Reads filtered.')
    print 'Number of FILTERED reads: ', len(read_container.fetch_all_reads(format=list))
    dataAccess.clear_cache()

    record_container = RecordContainer()
    record_container.set_db_access(dataAccess)
    record_container.populate(read_container.fetch_all_reads_versions(), table='rrna')
    raw_input('Records fetched')

    log.info('BINNER EXIT')



if __name__ == '__main__':
    main()
