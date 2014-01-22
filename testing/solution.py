from Bio import SeqIO
import re

class Solution(object):
    '''
    Containts correct solution of metagenomic sample.
    To each read assigns tax_id of its origin.
    '''
    def __init__(self, metasim_fasta_path, dataAccess):
        '''
        @param metasim_fasta_path   (string) Path to metasim fasta file
        @param dataAccess           (DataAccess) Db connection object
        '''
        self.id2taxon = {}

        handle = open(metasim_fasta_path, "rU")
        for read in SeqIO.parse(handle, "fasta"):
            try:
                gi = re.search("GI=(\d+)", read.description).group(1)
                taxon = dataAccess.get_taxids([gi], format=list)[0]
            except AttributeError:
                taxon = "none"
            self.id2taxon[read.id] = taxon

    def get_tax_id(self, read_id):
        return self.id2taxon[read_id]
