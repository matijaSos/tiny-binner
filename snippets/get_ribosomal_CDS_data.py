
import argparse
import sys,os
import time
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
            product     = cds_aln.cds.product # TODO: this one should be optional!
            # Quality measure data
            # measure = cds_aln.get_std_over_mean()
            # std     = cds_aln.get_coverage_std_dev()
            mean    = cds_aln.get_mean_coverage()
            length  = cds_aln.get_cds_length()
            f.write('{0:4} {1:10} {2:10d} {3:8.4f} {4:8d} {5:50}\n'.format(gene, 
                                                                        protein_id, tax_id, 
                                                                        mean, length,
                                                                        product))

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

def is_ribosomal(product):
    '''Determines if given product is of ribosomal gene.

    Args:
        product (string):   Gene product
    Returns:
        True if product is considered to be ribosomal,
        False otherwise
    '''
    # May have any of
    any_of  = ["ribosomal", "5S", "16S", "23S",
                            "5s", "16s", "23s"]
    # Must not have
    none_of = ["nonribosomal"]

    for elem in none_of:
        if elem in product: return False

    for elem in any_of:
        if elem in product: return True

    return False

def count_nones(cds_alns):
    '''Counts how many data in cds alignment data is missing.

    If it is None, then it is considered to be missing.
    Data checked are 'gene', 'protein_id' and 'product' attributes of each
    cds alignments are checked.

    Args:
        cds_alns    ([CdsAlignment]):   cds alignments objects
    Returns:
        dict(attr_name(string) -> count(int)): Number of Nones for each attribute.
    '''

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

def perc_format(text, part, total):
    return "{0}: {1:.2f} ({2})\n".format(text, part/float(total), part)


def export_species_data(species_data, total_reads, export_path, tax_tree, CDS_count=None):

    with open(export_path, 'w') as f:
        for tax_id in sorted(species_data, key=species_data.get, reverse=True):
            count = species_data[tax_id]
            try:
                rank = tax_tree.nodes[int(tax_id)].rank
                name = tax_tree.nodes[int(tax_id)].organism_name
            except KeyError:
                rank = None
                name = None
            frac = count / float(total_reads) * 100

            line = "{0:10} {1:10d} {2:80} {3:10.2f}% {4:10d}".format(rank, tax_id, name, frac, count)
            if CDS_count:
                cds_count = CDS_count[tax_id]
                line += "  CDSs: {0:5d}".format(cds_count)

            f.write(line + "\n")


        
def assignment_analysis(tax_ids, reads, tax_tree, export_folder, CDS_count=None):
    '''Analyses how god given reads "fit" to given tax_ids.
       For those reads who cannot be assigned to any given tax_id,
       separate analysis is conducted.

    The analysis is conducted on the level of species.

    Args:
        tax_ids         (set(int)): Set of tax_ids of species
        reads           ([Read]):   Reads to analyse
        tax_tree        (TaxTree):  TaxTree used to get species level tax_ids
        export_folder   (String):   Path of folder where analysis data will be stored
        CDS_count       ({tax_id -> cds_num}): Holds number of CDSs of for each given tax_id
    Returns:
        None, but three files are creates as output
    '''

    species_given   = {}
    species_new     = {}

    # Stats
    reads_assigned  = 0
    no_aln_reads    = 0

    for read in reads:

        # Skip reads without alignments
        if not read.has_alignments():
            no_aln_reads += 1
            continue

        # Get given species to which this read "fits"
        fits = set()

        for aln in read.get_alignments():
            tax_id = tax_tree.get_parent_with_rank(aln.tax_id, 'species')
            if tax_id in tax_ids:
                fits.add(tax_id)

        if fits:
            reads_assigned += 1
            # Increase counter for "fit" species
            for tax_id in fits:
                species_given[tax_id] = species_given.get(tax_id, 0) + 1
        else:
            # Unassigned reads
            species_set = set()
            for aln in read.get_alignments():
                tax_id = tax_tree.get_parent_with_rank(aln.tax_id, 'species')
                species_set.add(tax_id)
            for tax_id in species_set:
                species_new[tax_id] = species_new.get(tax_id, 0) + 1

            
    # Store analysis data - for given species
    species_given_path = os.path.join(export_folder, "species_ribosomal.txt")
    export_species_data(species_given, len(reads), species_given_path, tax_tree, CDS_count)
    # For "new", unexpected species
    species_new_path = os.path.join(export_folder, "species_new.txt")
    export_species_data(species_new, len(reads), species_new_path, tax_tree)

    # Create summary
    assn_summary_path = os.path.join(export_folder, "assignment_summary.txt")
    with open(assn_summary_path, 'w') as f:

        total_reads = len(reads)
        f.write("Total reads: {0}\n".format(total_reads))
        f.write("\n")
        f.write("Number of ribosomal species: {0}\n".format(len(tax_ids)))
        f.write("Number of 'new' species: {0}\n".format(len(species_new)))
        f.write("\n")
        f.write(perc_format("Reads assigned to tax ", reads_assigned, total_reads))
        f.write("\n")
        f.write(perc_format("Reads assigned to other ", total_reads - no_aln_reads - reads_assigned, total_reads))
        f.write(perc_format("Reads with 0 alns ", no_aln_reads, total_reads))

        '''
        f.write("Mean aln number per read: {0:.2f}\n".format(avg_aln_num))
        f.write(perc_format("Reads with 0 alns ", no_aln_reads, total_reads))
        f.write("\n")
        f.write(perc_format("Reads assigned to tax ", readsAssignedToTax, total_reads))
        f.write(perc_format("Reads unassigned to tax ", readsUnassignedToTax, total_reads))
        '''



def main():
    '''
    Script to analyse  genes expressed for given tax id.
    '''

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()

    # Access database
    dataAccess = DataAccess(args)

    # ------------------ #

    print '1. Loading tax tree...'
    start = time.time()

    tax_tree = TaxTree()

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    print '2. Loading alignment file...'
    start = time.time()

    read_container = ReadContainer()
    read_container.load_alignment_data(args.alignment_file)
    #---SET TAXIDS FOR ALL ALIGNMENTS--#
    read_container.set_taxids(dataAccess)

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    # Create folder if does not exist
    if not os.path.exists(args.export_folder):
        os.makedirs(args.export_folder)
    # File for data analysis summary
    summary_path = os.path.join(args.export_folder, "CDSs_summary.txt")
    cds_summary = open(summary_path, 'w')

    if args.remove_host:
        print "Removing host..."
        start = time.time()

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
        # Set host-free reads
        read_container.set_new_reads(reads_with_no_host_alignments)

        end = time.time()
        print("done: {0:.2f} sec".format(end - start))

    #------- LOAD ALL RECORDS   -------#

    print '4. Loading referenced records...'
    start = time.time()

    record_container = RecordContainer()
    record_container.set_db_access(dataAccess)
    record_container.populate(read_container.fetch_all_reads_versions(), table='cds')

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    #-- MAP ALIGNMENTS TO GENES   -----#

    print '5. Mapping alignments to genes...'
    start = time.time()

    read_container.populate_cdss(record_container)

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    #- RECORD ALL ALIGNEMENTS TO GENE -#

    print '6. Populating CDS container...'
    start = time.time()

    cds_aln_container = CdsAlnContainer()
    cds_aln_container.populate(read_container.fetch_all_reads(format=list))

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------------------- #

    print 'Sorting CDSs ...DISABLED'
    start = time.time()

    # Sort CDSs by their "good looks"!
    cds_alns = cds_aln_container.fetch_all_cds_alns(format=list)
    '''
    cds_alns = sorted(cds_alns,
                    key=lambda cds_aln: cds_aln.get_std_over_mean(),
                    reverse=False)
    '''
    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------------------- #

    '''
    print "Exporting phase 0 - all CDSs..."
    export_CDS_stats_data(cds_alns, args.export_folder, "0_all_CDSs.txt")
    print "done"
    '''

    # Count Nones in cds_alns
    nones = count_nones(cds_alns)
    cds_summary.write("\n")
    cds_summary.write("gene None       : {0}\n".format(nones['gene']))
    cds_summary.write("protein_id  None: {0}\n".format(nones['protein_id']))
    cds_summary.write("product  None   : {0}\n".format(nones['product']))
    cds_summary.write("\n")

    cds_summary.write("CDSs all: {0}\n".format(len(cds_alns)))

    print 'Filtering valid CDSs...'
    start = time.time()

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

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))


    # All valid CDSs - Output coverage/length histogram data
    print "Exporting phase 1 - all CDSs..."
    start = time.time()

    export_CDS_stats_data(cds_alns_targeted, args.export_folder, "1_all_valid_CDSs.txt")

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

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



        # TODO: make a function for this
        if is_ribosomal(product):
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

    # ------------------- Making biological sense - choosing CDSs -------------------- #

    print 'Filtering under-average ribosomals...'
    # NOTE: take length into consideration?
    cds_alns_ribosomal = [cds_aln for cds_aln in cds_alns_ribosomal
                         # Filters
                         if cds_aln.get_mean_coverage() > mm_cov]
    print 'done'
    cds_summary.write("ribosomals over-mean: {0}\n".format(len(cds_alns_ribosomal)))
    cds_summary.close()

    print 'Phase 3 - filtered ribosomal CDSs...'
    export_CDS_stats_data(cds_alns_ribosomal, args.export_folder, "3_ribosomal_CDSs_filtered.txt")
    print 'done'
    
    # Store charts cov data - if selected so
    if args.export_charts:
        print "Exporting chart coverage data..."
        export_CDS_graph_data(cds_alns_ribosomal, args.export_charts)
        print "done."

    # --------------------- I have chosen CDSs - determine species and analyse ------------------------ #

    # Species level resolution
    # See which species are present - dump ones with not enough CDSs
    # NOTE: So far done in determine_species_by_ribosomals.py

    CDS_count   = {}    # Count CDSs of each species
    species_set = set() # Get estimated tax_ids
    for cds_aln in cds_alns_ribosomal:
        tax_id = cds_aln.cds.taxon

        # Put each tax_id up to the "species" level
        tax_id_species = tax_tree.get_parent_with_rank(tax_id, 'species')

        species_set.add(tax_id_species)
        CDS_count[tax_id_species] = CDS_count.get(tax_id_species, 0) + 1

    # Get reported CDSs ids
    reported_CDS_ids = set()
    for cds_aln in cds_alns_ribosomal:
        reported_CDS_ids.add(cds_aln.cds.id)

    # ------------ Read assignment analysis -------------- #

    print "Read assignment analysis..."

    reads = read_container.fetch_all_reads(format=list)
    assignment_analysis(species_set, reads, tax_tree, args.export_folder, CDS_count)

    '''
    readsAssignedToTax   = 0
    readsUnassignedToTax = 0

    readsAsnToCDS   = 0
    readsUnasnToCDS = 0

    readsAsnToCDS_first = 0

    # Average number of alignments per read
    avg_aln_num = 0
    # Reads without alns
    no_aln_reads = 0

    # Check each read and see if has any alignment to some tax_id from speciesSet
    reads = read_container.fetch_all_reads(format=list)
    for read in reads:
        assignedToTax = False
        assignedToCDS = False
        on_first_aln  = False # TODO: use array for histogram

        aln_num = len(read.get_alignments())
        avg_aln_num += aln_num
        # Count reads without alns
        if aln_num == 0:
            no_aln_reads += 1

        aln_idx = 0
        for aln in read.get_alignments():
            aln_idx += 1

            # Check if maps to the reported tax_id - put up to species
            tax_id = aln.tax_id
            tax_id_species = tax_tree.get_parent_with_rank(tax_id, 'species')
            if tax_id_species in species_set:
                assignedToTax = True

            # Check if maps to the reported CDS
            cdss = [pair[0] for pair in aln.aligned_cdss] # pair = (cds, location)
            for cds in cdss:
                if cds.id in reported_CDS_ids:
                    assignedToCDS = True
                    # Store aln index - is first?
                    if aln_idx == 1:
                        on_first_aln = True
                        
        # Count state
        # TODO: use total reads count
        if assignedToTax:
            readsAssignedToTax += 1
        else:
            readsUnassignedToTax += 1

        if assignedToCDS:
            readsAsnToCDS += 1
        else: 
            readsUnasnToCDS += 1

        if on_first_aln:
            readsAsnToCDS_first += 1

    avg_aln_num = float(avg_aln_num) / read_container.get_read_count()

    # Write to assignment summary to special file
    assn_summary_path = os.path.join(args.export_folder, "assignment_summary.txt")

    with open(assn_summary_path, 'w') as f:

        total_reads = read_container.get_read_count()
        f.write("Total reads: {0}\n".format(total_reads))
        f.write("Number of species: {0}\n".format(len(species_set)))

        f.write("Mean aln number per read: {0:.2f}\n".format(avg_aln_num))
        f.write(perc_format("Reads with 0 alns ", no_aln_reads, total_reads))
        f.write("\n")
        f.write(perc_format("Reads assigned to tax ", readsAssignedToTax, total_reads))
        f.write(perc_format("Reads unassigned to tax ", readsUnassignedToTax, total_reads))
        f.write("\n")
        f.write(perc_format("Reads assigned to CDS ", readsAsnToCDS, total_reads))
        f.write(perc_format("Asb to CDS on 1st aln", readsAsnToCDS, total_reads))

    print "done."
    '''



    # ----------------------- Read assigning phase -------------------------- #

    # Go through all reads and look who can be assigned where

    # TO BE IMPLEMENTED

    # read_container

if __name__ == '__main__':
    main()