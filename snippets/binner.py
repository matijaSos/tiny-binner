import logging
import logging.config
import argparse
import sys,os
sys.path.append(os.getcwd())

from utils.argparser import *
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer
import filters.host as host_filter

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

    #----------------------------------#
    #------- ALIGNMENT DATA SOURCE ----#
    read_container = ReadContainer()
    read_container.load_alignment_data(args.input)
    read_container.set_taxids(dataAccess)

    #----------------------------------#
    #-------- TAXONOMY TREE -----------#
    tax_tree = TaxTree() 
    # tax_tree.load_taxonomy_data(dataAccess) 

    host_filter.filter_potential_hosts_alignments(
        read_container.fetch_all_reads(), 
        tax_tree.tax2relevantTax, 
        tax_tree.potential_hosts, 
        delete_host_alignments = False, 
        filter_unassigned = True, 
        unassigned_taxid=-1)



    log.info('BINNER EXIT')
    


if __name__ == '__main__':
    main()
