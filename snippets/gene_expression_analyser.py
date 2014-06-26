
import argparse
import sys,os
sys.path.append(os.getcwd()) # What's the purpose of this?

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree
from data.containers.read import ReadContainer

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
        self.add_argument('tax_id',
                help='Taxonomy id to analyse',
                type=int)
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
    for cds_aln in cds_alns_sorted:
        filename = "cds_" + str(i) + ".txt"
        coverage_path = os.path.join(export_path, filename)

        print(str(i) + ": " + str(cds_aln.get_std_over_mean()))

        cds_aln.coverage_to_file(coverage_path)

        if i == 50: # TODO: Define this somehow as parameter
            break

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

    print '1. Loading tax tree...'
    tax_tree = TaxTree()
    print 'done.'

    print '2. Loading alignment file...'
    read_container = ReadContainer()
    read_container.load_alignment_data(args.alignment_file)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)
    print 'done'

    # TODO: Here i should recognize host reads!

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
    #----------------------------------#
    #- RECORD ALL ALIGNEMENTS TO GENE -#
    print '6. Populating CDS container...'
    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))
    print 'done'

    print("Loaded CDS container")
    
    # Take only CDSs of given tax_id
    # Remove CDSs with too low mean coverage value
    min_mean_coverage   = 10
    min_length          = 20

    cds_alns = cds_aln_container.fetch_all_cds_alns(format=list)
    print ( "All CDSs (all organisms): " + str(len(cds_alns)) )

    cds_alns_targeted = [cds_aln for cds_aln in cds_alns 
                         if  cds_aln.get_tax_id() == args.tax_id
                         # Filters
                         and cds_aln.get_cds_length() > min_length
                         and cds_aln.get_mean_coverage() > min_mean_coverage]

    # Remove CDSs with no gene/product
    cds_alns_targeted = [cds_aln for cds_aln in cds_alns_targeted
                         if  cds_aln.cds.gene != None
                         and cds_aln.cds.product != None]

    # ------------------- CDSs filtered and ready to be analyzed ------------------- #

    # Number of targeted CDSs
    print ( "Targeted CDSs: " + str(len(cds_alns_targeted)) )

    '''
    print ("Sorting CDSs: stddev/mean")
    cds_alns_sorted = sorted(cds_alns_targeted, 
                             key=lambda cds_aln: cds_aln.get_std_over_mean(),
                             reverse=False)

    # TODO: Here I should somehow determine which CDSs are "expressed", and which are not?

    # Write to file stuff(gene, protein_id) for each cds_aln
    print("Writing data to file")
    path = args.export_path
    f = open(path, 'w')
    
    for cds_aln in cds_alns_targeted:
        gene        = cds_aln.cds.gene
        product     = cds_aln.cds.product
        protein_id  = cds_aln.cds.protein_id

        f.write("{0} {1}\n".format(gene, protein_id))
    f.close()

    print("Done")
    '''
    
    # -------------------- #
    

    '''
    # Analyse those CDSs
    print ( "Targeted CDSs: " + str(len(cds_alns_targeted)) )

    # See the mean length of CDS
    mean_cds_length = 0
    no_locs_num = 0
    for cds_aln in cds_alns_targeted:
        try:
            loc_length = cds_aln.get_cds_length()
            mean_cds_length += loc_length
        except:
            no_locs_num += 1
    # Get mean
    mean_cds_length /= float(len(cds_alns_targeted))

    print("---------------------------------------------")
    print("Mean CDS length: " + str(mean_cds_length))
    print("Nones: " + str(no_locs_num))
            
    print("---------------------------------------------")
    '''


    '''
    # Create folder where data about CDSs will be stored
    if not os.path.exists(args.export_path):
        os.makedirs(args.export_path)

    # Export some amount of best CDSs
    i = 1
    for cds_aln in cds_alns_sorted:
        filename = "cds_" + str(i) + ".txt"
        coverage_path = os.path.join(args.export_path, filename)

        print(str(i) + ": " + str(cds_aln.get_std_over_mean()))

        cds_aln.coverage_to_file(coverage_path)

        if i == 50: # TODO: Define this somehow as parameter
            break

        i += 1

    # Load CDS container
    '''

    # Analyse stuff
    # print("Analysing stuff!")

    '''
    print '4. Creating LCA solution...'
    lca_binner = LCABinner(tax_tree)
    lca_sol = lca_binner.bin_reads(read_container)
    print 'done'

    print '5. Nice output of read assignment here'
    lca_sol.print_nicely(tax_tree)
    '''

if __name__ == '__main__':
    main()
