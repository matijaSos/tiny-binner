
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
        self.add_argument('export_path',
                help='Folder to export data to',
                type=str)

def export_CDS_graph_data(cds_alns, export_path):
    '''Export coverage data of each cds_alignment to the specified folder.

    This data is later used to create coverage graphs

    Args:
        cds_alns    ([CdsAlignment]):   cds alignments objects
        export_path (string):           folder to export data to
    '''

    # Create folder where data about CDSs will be stored
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
    

def main():
    '''
    Script to analyse  genes expressed for given tax id.
    '''

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()

    # Access database
    dataAccess = DataAccess(args)

    #print '1. Loading tax tree...'
    tax_tree = TaxTree()
    #print 'done.'

    #print '2. Loading alignment file...'
    read_container = ReadContainer()
    read_container.load_alignment_data(args.alignment_file)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)
    #print 'done'

    '''
    # TODO: Here i should recognize host reads!
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
    
    print ("total   : {0:8d}".format(read_count))
    print ("host    : {0:8d} {1:.2f}".format(host_read_count, 
                                  host_read_count / float(read_count)
                                  ))
    print ("non-host: {0:8d} {1:.2f}".format(non_host_read_count, 
                                  non_host_read_count / float(read_count)
                                  ))
    print
    read_container.set_new_reads(reads_with_no_host_alignments)
    '''

    # ------------------------------------- #

    #----------------------------------#
    #------- LOAD ALL RECORDS   -------#
    #print '4. Loading referenced records...'
    record_container = RecordContainer()
    record_container.set_db_access(dataAccess)
    record_container.populate(read_container.fetch_all_reads_versions(), table='cds')
    #print 'done'
    #----------------------------------#
    #-- MAP ALIGNMENTS TO GENES   -----#
    #print '5. Mapping alignments to genes...'
    read_container.populate_cdss(record_container)
    #----------------------------------#
    #- RECORD ALL ALIGNEMENTS TO GENE -#
    #print '6. Populating CDS container...'
    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))
    #print 'done'

    #print("Loaded CDS container")
    
    # Take only CDSs of given tax_id
    # Remove CDSs with too low mean coverage value
    min_mean_coverage   = 0
    min_length          = 0

    cds_alns = cds_aln_container.fetch_all_cds_alns(format=list)
    print ( "CDSs all  : " + str(len(cds_alns)) )

    cds_alns_targeted = [cds_aln for cds_aln in cds_alns 
                         # Filters
                         if cds_aln.get_cds_length() > min_length
                         and cds_aln.get_mean_coverage() > min_mean_coverage]

    # Remove CDSs with no gene/product
    cds_alns_targeted = [cds_aln for cds_aln in cds_alns_targeted
                         if  cds_aln.cds.gene != None
                         and cds_aln.cds.product != None]

    # ------------------- CDSs filtered and ready to be analyzed ------------------- #

    # Number of targeted CDSs
    print ( "CDSs valid: " + str(len(cds_alns_targeted)) )
    print

    cds_alns_ribosomal = []
    for cds_aln in cds_alns_targeted:

        # If has word "ribosomal" in name, store coverage data for graph
        gene        = cds_aln.cds.gene
        product     = cds_aln.cds.product
        protein_id  = cds_aln.cds.protein_id

        if "ribosomal" in product:
            #print("{0} {1} {2}\n".format(gene, protein_id, product))
            cds_alns_ribosomal.append(cds_aln)

    # ------------------- Ribosomal CDSs acquired! --------------------- #
    # Sort it!
    cds_alns_ribosomal = sorted(cds_alns_ribosomal, 
                             key=lambda cds_aln: cds_aln.get_std_over_mean(),
                             reverse=False)

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

    # Print
    print("ribosomals: " + str(len(cds_alns_ribosomal)))
    print("mean coverage: " + str(mm_cov))
    print("max coverage : {0}".format(max_cov))
    print
    for cds_aln in cds_alns_ribosomal:
        gene        = cds_aln.cds.gene
        product     = cds_aln.cds.product
        protein_id  = cds_aln.cds.protein_id

        taxon       = cds_aln.cds.taxon
        name        = tax_tree.nodes[taxon].organism_name
        print("{0:4} {1:10} {2:50} {3:10d} {4:60}".format(gene, protein_id, product, taxon, name))

    # Store graph data
    export_CDS_graph_data(cds_alns_ribosomal, args.export_path)

    '''
    # Write to file stuff(gene, protein_id) for each cds_aln
    print("Writing data to file")
    path = args.export_path
    f = open(path, 'w')
    for cds_aln in cds_alns_targeted:

        gene        = cds_aln.cds.gene
        product     = cds_aln.cds.product
        protein_id  = cds_aln.cds.protein_id

        # If has word "ribosomal" in name, store coverage data for graph

        f.write("{0} {1} {2}\n".format(gene, protein_id, product))
    f.close()
    '''

if __name__ == '__main__':
    main()
