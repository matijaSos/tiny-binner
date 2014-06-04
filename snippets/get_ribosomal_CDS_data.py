
import argparse
import sys,os
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
                takes all reads assigned to given tax id and \
                identifies "hit" genes in that organism.
                .''')
        self.add_argument('alignment_file',
                help='path to the read alignment file',
                type=str)
        self.add_argument('export_folder',
                help='Folder to write ribosomal CDS analysis data to',
                type=str)


        # -------- Optional args -------- #

        self.add_argument('--remove_host', 
                help='remove host reads and alignments',
                action='store_true')

        self.add_argument('--export_charts', 
                help='folder export charts with CDS coverage')

def export_CDS_stats_data(cds_alns, export_folder, export_file):
    '''Write row containing data for each given CDS.
       Data: gene, protein.id, tax_id, stddev/mean_cov, std, mean_cov, length

    Args:
        cds_alns        ([CdsAlignment]):   Array of cds alignments objects
        export_folder   (string):           Folder to export the data to
        export_file     (string):           File name to save the data to
    '''

    # Create folder if does not exist
    if not os.path.exists(export_folder):
        os.makedirs(export_folder)

    export_path = os.path.join(export_folder, export_file)
    with open(export_path, 'w') as f:
        for cds_aln in cds_alns:

            gene        = cds_aln.cds.gene
            if gene is not None:
                gene = gene.replace(' ', '-') # Get rid of spaces

            protein_id  = cds_aln.cds.protein_id
            tax_id      = cds_aln.cds.taxon
            # Quality measure data
            measure = cds_aln.get_std_over_mean()
            std     = cds_aln.get_coverage_std_dev()
            mean    = cds_aln.get_mean_coverage()
            length  = cds_aln.get_cds_length()
            f.write('{0:4} {1:10} {2:10d} {3:10.4f} {4:8.4f} {5:8.4f} {6:8d}\n'.format(gene, 
                                                                        protein_id, tax_id, 
                                                                        measure, std, mean, length))

def export_CDS_graph_data(cds_alns, export_path):
    '''Export coverage data of each cds_alignment to the specified folder.

    This data is later used to create coverage graphs

    Args:
        cds_alns    ([CdsAlignment]):   cds alignments objects
        export_path (string):           folder to export data to
    '''

    # Create folder where data about CDSs charts will be stored
    if not os.path.exists(export_path):
        os.makedirs(export_path)

    # Export some amount of best CDSs
    i = 1
    for cds_aln in cds_alns:
        filename = "cds_" + str(i) + ".txt"
        coverage_path = os.path.join(export_path, filename)

        #print(str(i) + ": " + str(cds_aln.get_std_over_mean()))

        cds_aln.coverage_to_file(coverage_path)
        i += 1

def count_nones(cds_alns):

    nones = {}
    nones['gene'] = 0
    nones['protein_id'] = 0
    nones['product'] = 0

    for cds_aln in cds_alns:
        if cds_aln.cds.gene is None:
            nones['gene'] += 1
        if cds_aln.cds.protein_id is None:
            nones['protein_id'] += 1
        if cds_aln.cds.product is None:
            nones['product'] += 1
        
    return nones


def main():
    '''
    Script to analyse  genes expressed for given tax id.
    '''

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()

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

    # Create folder if does not exist
    if not os.path.exists(args.export_folder):
        os.makedirs(args.export_folder)
    # File for data analysis summary
    summary_path = os.path.join(args.export_folder, "CDSs_summary.txt")
    cds_summary = open(summary_path, 'w')

    if args.remove_host:
        print "Removing host..."
        #------- FILTER HOST READS -------#
        #print '3. Filtering host reads & alignments...'
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
            host_filter.perc_of_host_alignments_larger_than)

        dataAccess.clear_cache()    # deletes gi2taxid cache

        reads_with_no_host_alignments = host_filter.filter_potential_hosts_alignments(
            new_reads,
            tax_tree.tax2relevantTax,
            tax_tree.potential_hosts,
            True,   # delete host alignments
            True,   # filter unassigned
            -1)     # unassigned taxid

        read_count          = len(read_container.fetch_all_reads(format=list))
        host_read_count     = read_count - len(reads_with_no_host_alignments)
        non_host_read_count = read_count - host_read_count
        
        cds_summary.write("total   : {0:8d}\n".format(read_count))
        cds_summary.write("host    : {0:8d} {1:.2f}\n".format(host_read_count, 
                                      host_read_count / float(read_count)
                                      ))
        cds_summary.write("non-host: {0:8d} {1:.2f}\n".format(non_host_read_count, 
                                      non_host_read_count / float(read_count)
                                      ))

        read_container.set_new_reads(reads_with_no_host_alignments)
        print "done"

    # ------------------------------------- #

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
    print 'done'
    #----------------------------------#
    #- RECORD ALL ALIGNEMENTS TO GENE -#
    print '6. Populating CDS container...'
    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))
    print 'done'

    cds_alns = cds_aln_container.fetch_all_cds_alns(format=list)

    # Sort CDSs by their "good looks"!
    cds_alns = sorted(cds_alns,
                    key=lambda cds_aln: cds_aln.get_std_over_mean(),
                    reverse=False)

    print "Exporting phase 0 - all CDSs..."
    export_CDS_stats_data(cds_alns, args.export_folder, "0_all_CDSs.txt")
    print "done"

    # Count Nones in cds_alns
    nones = count_nones(cds_alns)
    cds_summary.write("gene None       : {0}\n".format(nones['gene']))
    cds_summary.write("protein_id  None: {0}\n".format(nones['protein_id']))
    cds_summary.write("product  None   : {0}\n".format(nones['product']))

    cds_summary.write("CDSs all: {0}\n".format(len(cds_alns)))

    print 'Filtering valid CDSs...'
    # Take only CDSs of given tax_id
    # Remove CDSs with too low mean coverage value
    min_mean_coverage   = 0
    min_length          = 0
    cds_alns_targeted = [cds_aln for cds_aln in cds_alns 
                         # Filters
                         if cds_aln.get_cds_length() > min_length
                         and cds_aln.get_mean_coverage() > min_mean_coverage]

    # Remove CDSs with no gene/product
    cds_alns_targeted = [cds_aln for cds_aln in cds_alns_targeted
                         if cds_aln.cds.product is not None]
                         #if  cds_aln.cds.gene != None
                         #and cds_aln.cds.product != None]
    print 'done'


    # All valid CDSs - Output coverage/length histogram data
    print "Exporting phase 1 - all CDSs..."
    export_CDS_stats_data(cds_alns_targeted, args.export_folder, "1_all_valid_CDSs.txt")
    print "done"

    # ------------------- CDSs filtered and ready to be analyzed ------------------- #

    print 'Extracting ribosomal CDSs...'
    # Number of targeted CDSs
    cds_summary.write("CDSs valid: {0}\n".format(len(cds_alns_targeted)))

    cds_alns_ribosomal = []
    for cds_aln in cds_alns_targeted:

        # If has word "ribosomal" in name, store coverage data for graph
        gene        = cds_aln.cds.gene
        product     = cds_aln.cds.product
        protein_id  = cds_aln.cds.protein_id

        if "ribosomal" in product:
            #print("{0} {1} {2}\n".format(gene, protein_id, product))
            cds_alns_ribosomal.append(cds_aln)

    print 'done'
    # ------------------- Ribosomal CDSs acquired! --------------------- #

    print 'Analysing ribosomals...'

    # Extract interesting data
    # Mean coverage, max coverage
    mm_cov  = 0
    max_cov = 0
    for cds_aln in cds_alns_ribosomal:
        mean_cov = cds_aln.get_mean_coverage()
        mm_cov += mean_cov
        max_cov = max(max_cov, mean_cov)
    if mm_cov > 0:
        mm_cov /= len(cds_alns_ribosomal)

    cds_summary.write("ribosomals all {0}\n".format(len(cds_alns_ribosomal)))
    cds_summary.write("mean coverage: {0}\n".format(mm_cov))
    cds_summary.write("max coverage : {0}\n".format(max_cov))
    print 'done'

    # Ribosomal CDSs only - Output coverage/length histogram
    print "Exporting phase 2 - ribosomal CDSs only..."
    export_CDS_stats_data(cds_alns_ribosomal, args.export_folder, "2_ribosomal_CDSs.txt")
    print "done"

    # ------------------- Making biological sense -------------------- #

    print 'Filtering under-average ribosomals...'
    # NOTE: take length into consideration?
    cds_alns_ribosomal = [cds_aln for cds_aln in cds_alns_ribosomal
                         # Filters
                         if cds_aln.get_mean_coverage() > mm_cov]
    print 'done'
    cds_summary.write("ribosomals over-mean: {0}\n".format(len(cds_alns_ribosomal)))

    # Filtered ribosomal CDSs only - Output coverage/length histogram

    # Species level resolution
    # See with species are present - dump ones with not enough CDSs
    # NOTE: So far done in determine_species_by_ribosomals.py

    print 'Phase 3 - filtered ribosomal CDSs...'
    export_CDS_stats_data(cds_alns_ribosomal, args.export_folder, "3_ribosomal_CDSs_filtered.txt")
    print 'done'
    
    '''
    for cds_aln in cds_alns_ribosomal:
        gene        = cds_aln.cds.gene
        protein_id  = cds_aln.cds.protein_id
        tax_id      = cds_aln.cds.taxon
        # Quality measure data
        measure = cds_aln.get_std_over_mean()
        std     = cds_aln.get_coverage_std_dev()
        mean    = cds_aln.get_mean_coverage()
        length  = cds_aln.get_cds_length()
        cds_data.write('{0:4} {1:10} {2:10d} {3:10.4f} {4:8.4f} {5:8.4f} {6:8d}\n'.format(gene, 
                                                                    protein_id, tax_id, 
                                                                    measure, std, mean, length))
    '''
    cds_summary.close()
    
    # Store graph data
    if args.export_charts:
        export_CDS_graph_data(cds_alns_ribosomal, args.export_charts)

    # ----------------------- Read assigning phase -------------------------- #

    # Go through all reads and look who can be assigned where

if __name__ == '__main__':
    main()
