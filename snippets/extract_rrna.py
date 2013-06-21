import time
import logging
import logging.config
import argparse
import sys,os
sys.path.append(os.getcwd())

from utils.argparser import *
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer
from data.containers.record import RecordContainer
import filters.host as host_filter
from utils import timeit

def log_start (log, args):
    log.info('BINNER RUN')
    log.info("Input: %s" % args.input)
    log.info("Xml template: %s" % args.descr)
    log.info("Output: %s" %  args.output)



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
    argparser = get_binner_argparser()
    args  = argparser.parse_args()
    error = validate_args(args)
    if error: exit(-1)

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

    for record in record_container.fetch_all_records():
        if record:
            print record.version



    log.info('BINNER EXIT')



if __name__ == '__main__':
    main()
