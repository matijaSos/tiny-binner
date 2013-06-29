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
    tax_tree = TaxTree()
    # tax_tree.load_taxonomy_data(dataAccess)
    #raw_input('Tax tree loaded.')

    #----------------------------------#
    #------- ALIGNMENT DATA SOURCE ----#
    read_container = ReadContainer()
    read_container.load_alignment_data(args.input)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)

    #------- FILTER HOST READS -------#
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
    dataAccess.clear_cache()    # deletes gi2taxid cache

    #----------------------------------#
    #------- LOAD ALL RECORDS   -------#
    record_container = RecordContainer()
    record_container.set_db_access(dataAccess)
    record_container.populate(read_container.fetch_all_reads_versions(), table='rrna')
    #----------------------------------#
    #-- MAP ALIGNMENTS TO GENES   -----#
    read_container.populate_cdss(record_container)
    #----------------------------------#
    #- RECORD ALL ALIGNEMENTS TO GENE -#
    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))

    output_file = open(args.output, 'w')

    for cds_aln in cds_aln_container.fetch_all_cds_alns():
        output = cds_aln.cds.version + ','
        if cds_aln.cds.gene is not None:
            output += 'gene:%s,' % cds_aln.cds.gene
        location = Location.from_location_str(cds_aln.cds.location)
        if location.start is not None:
            aln_loc = '(%d,%d)' % (location.start, location.end)
        else:
            aln_loc = ''
            for subloc in location.sublocations:
                aln_loc += '(%d,%d),' % (subloc.start, subloc.end)
        output += '%s:' % aln_loc
        for aln_subloc in cds_aln.aligned_regions.values():
            if aln_subloc.location.start is None or aln_subloc.location.end is None:
                aln_str = ''
                for subloc in aln_subloc.location.sublocations:
                    aln_str += '(%d,%d),' % (subloc.start, subloc.end)
            else: aln_str = '(%d,%d)' % (aln_subloc.location.start, aln_subloc.location.end)
            output += "%s,%.1f,%s;" % (aln_subloc.read_id, aln_subloc.score, aln_str)
        output_file.write('%s\n' % output)

    output_file.close()
    log.info('BINNER EXIT')



if __name__ == '__main__':
    main()
