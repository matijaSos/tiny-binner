from data.read import Read
from utils.location import Location

import time

class ReadContainer (object):
    ''' Contains all the reads loaded from an 
        alignment file. Can be queried by read id.
    '''
    def __init__(self):
        """
        (dict) read_repository Dictionary where value is (Read)read and key is (str)read id.
        """
        self.read_repository = {}
        
    def load_alignment_data (self, read_alignment_file):
        ''' Adds all the reads in the alignment file to the
            read repository.
            This is the first stage of filling the read container.
        '''
        aln_file = open(read_alignment_file, 'r')
        for line in aln_file.readlines():
            self._add_read_from_str(line)

    def set_taxids (self, data_access):
        for read in self.fetch_all_reads(format=iter):
            for alignment in read.get_alignments():
                taxid = data_access.get_taxids([alignment.genome_index])
                setattr(alignment, 'tax_id', taxid)

    def get_protein_ids(self, exclude_host=False):
        protein_ids = set([])
        for read in self.read_repository.values():
            if exclude_host:
                if read.is_host_read:
                    continue
            for readAln in read.alignment_locations:
                for (cds, alignment_location) in readAln.aligned_cdss:
                    protein_ids.add(cds.protein_id)
        return protein_ids

    def populate_cdss (self, record_container):
        '''
        Coding sequences are determined and stored for every read alignment.
        Prerequisite: record container has been populated with all records
        mentioned in the alignment file
        @param record_container (RecordContainer)
        '''
        
        for read in self.fetch_all_reads(format=iter):
            for read_alignment in read.get_alignments(format=iter):
                record = record_container.fetch_existing_record(
                    read_alignment.nucleotide_accession)
                read_alignment.determine_coding_seqs_optimal(record)

    
    def fetch_read (self, read_id):
        if self.read_repository.has_key(read_id):
            return self.read_repository[read_id]
        else:
            raise KeyError("Read repository doesn't contain read associated with read ID: {0}".format(read_id))


    def fetch_all_reads (self, format=iter):
        return format(self.read_repository.values())
    
    def fetch_all_reads_versions(self):
        '''
        Returns an iterator returning all versions of all reads in this
        container
        
        .. note::
            Duplicate values are not filtered out
            
        :returns: iterator returning versions of reads
        '''
        reads = self.fetch_all_reads(format=iter)
        for read in reads:
            for read_alignment in read.get_alignments(format=iter):
                yield read_alignment.nucleotide_accession
        
    
    def _add_read_from_str (self, read_str):
        read = Read.from_read_str(read_str)
        assert (not self.read_repository.has_key(read.id))
        self.read_repository[read.id] = read
        
