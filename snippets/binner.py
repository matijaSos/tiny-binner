import logging
import logging.config
import argparse
import sys,os
sys.path.append(os.getcwd())

from utils.argparser import *
from ncbi.db.data_access import DataAccess
from data.containers.read import ReadContainer

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
    for read in read_container.fetch_all_reads(format=list):
        for aln in read.get_alignments():
            print aln.taxid


    log.info('BINNER EXIT')
    


if __name__ == '__main__':
    main()