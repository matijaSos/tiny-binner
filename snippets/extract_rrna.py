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

class ExtractRnasParser(DefaultBinnerArgParser):
    def __init__(self):
        super(ExtractRnasParser, self).__init__('''Loads input alignment file,\
        ... filteres host reads and produces file with all records found \
        ...among rrnas.'''
        )
        self.add_argument('input', help='Input alignment file', 
                           type=str)
        self.add_argument('output', help='Output file with list of found rna accessions.',
                           type=str)
        self.add_argument('-dh', '--determine_host', help='determine host method',
            choices=['best_score_host', '50perc_host', 'all_host'],
            default = '50perc_host')


def log_start (log, args):
    log.info('EXTRACT RRNA RUN')
    log.info("Input: %s" % args.input)
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
    # 1. input alignment file
    # 2. output rrna accession file
    argparser = ExtractRnasParser()
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
    timeit(read_container.load_alignment_data, input_aln_file)
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

    output_rrna_accession_file = open(output_rrna_accession_file, 'w')
    for record in record_container.fetch_all_records():
        if record:
            output_rrna_accession_file.write('%s\n' % (record.version))
    output_rrna_accession_file.close()
    log.info('BINNER EXIT')



if __name__ == '__main__':
    main()
