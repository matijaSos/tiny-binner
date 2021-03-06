from Bio import SeqIO
import re
import csv

class Solution(object):
    '''
    Contains solution of metagenomic sample.
    To each read single tax_id is assigned.
    '''
    def __init__(self, id2taxon):
        '''Constructor

        Args:
            id2taxon (dict): read.id -> taxon
        '''
        self.id2taxon = id2taxon

    def get_tax_id(self, read_id):
        '''Fetches tax_id assigned to a given read

        Args:
            read_id (string): id of a read
        Returns:
            (int): tax_id
        '''
        return self.id2taxon[read_id]

    def add_assignment(self, read_id, tax_id):
        '''Adds assignment of read to taxon.

        Args:
            read_id (string): Read id
            tax_id (int): Taxon id
        '''
        self.id2taxon[read_id] = tax_id;

    # ------------------------- Output -------------------------- #

    def print_data(self):
        '''
        Outputs inner dictionary: read.id -> taxon.
        '''
        for idx,taxon in self.id2taxon.iteritems():
            print str(idx) + " " + str(taxon)

    def get_taxon_count(self):
        '''Returns for each taxon how much reads is assigned to it.

        Args:
            None
        Returns:
            (dict): tax_id -> count
            
        '''
        taxon_count = {}
        for read_id, taxon in self.id2taxon.iteritems():
            taxon_count[taxon] = taxon_count.get(taxon, 0) + 1

        return taxon_count

    def get_taxon_shares(self):
        '''Returns share(abundance) of each taxon in read assignment.

        Args:
            None
        Returns:
            (dict): tax_id -> share
        '''
        taxon_count = self.get_taxon_count()

        # Get shares - divide by number of reads
        num_reads = len(self.id2taxon)
        taxon_shares = {taxon:(float(count)/num_reads) 
                        for taxon, count in taxon_count.items()}
        
        return taxon_shares;

    def print_nicely(self, tax_tree):
        '''Outputs share and count of each taxon.
           
        Output is sorted, starting from the taxon with the highest share(count).

        Args:
            tax_tree (TaxTree): Taxonomy tree, containing tax_id -> name mapping
        Returns:
            None
        '''
        taxon_shares = self.get_taxon_shares()
        taxon_count  = self.get_taxon_count()

        # Sort by share, higher first
        sorted_tax_shares = sorted(taxon_shares.iteritems(), 
                                   key      = lambda t:t[1],
                                   reverse  = True)

        total_share = 0
        total_count = 0
        # Print data for each tax id
        for tax_id, share in sorted_tax_shares:
            count = taxon_count[tax_id]
            name  = tax_tree.nodes[tax_id].organism_name
            rank  = tax_tree.nodes[tax_id].rank

            print("{0:10d} {1:20} {2:60} {3:10.4f}% {4:10d}".format(
                    tax_id,
                    rank,
                    name,
                    share*100,
                    count)
                  )
            total_share += share
            total_count += count

        print "-------------"

        print("total_share: {0}%".format(total_share*100))
        print("total_count: {0}".format(total_count))


    # --------------------- Factory methods ---------------------- #

    @staticmethod
    def from_metasim_fasta(path, dataAccess):
        '''Loads read assignment from FASTA file generated by metasim.

        Args:
            path (string): Path to the metasim fasta file
            dataAccess (DataAccess): Db connection object
        Returns:
            (Solution): New instance containing loaded data
        '''
        # For each read, gi of origin is extracted from header 
        # and mapped to tax_id.
        id2taxon = {}

        handle = open(path, "rU")
        for read in SeqIO.parse(handle, "fasta"):
            try:
                gi = re.search("GI=(\d+)", read.description).group(1)
                taxon = dataAccess.get_taxids([gi], format=list)[0]
            except Exception:
                taxon = "none"
            id2taxon[read.id] = taxon

        handle.close()
        return Solution(id2taxon)

    @staticmethod
    def from_CSV(path):
        '''Loads read assignment from CSV file.

        Assumption: CSV is comma separated
        Row format: "read_id,tax_id"

        Args:
            path (string): Path to the CSV file
        Returns:
            (Solution): New Solution instance
        '''
        id2taxon = {}

        handle = open(path, "rb")
        reader = csv.reader(handle)
        for row in reader:
            read_id = row[0]
            tax_id  = int(row[1])
            id2taxon[read_id] = tax_id
        
        handle.close()
        return Solution(id2taxon)

    @staticmethod
    def create_empty():
        ''' Creates empty Solution with no assignments.
        '''
        return Solution({})

