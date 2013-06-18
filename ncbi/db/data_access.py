from utils import enum
from ncbi.db.ncbitax_from_file import *
from ncbi.db.access import DbQuery

class DataAccess ():
    '''
    Provides access to CDS data and NCBItax data.
    Encapsulates data access, since all data
    can be accessed in two ways: using MySql database 
    or loading raw data from files.
    This method must provide all methods that would be 
    available via standard database querying.
    '''

    load_type = enum (FILE=1,
                      DATABASE=2
                     )

    def __init__(self, args, lazy_load=False):
        '''
        :param Any object containing following arguments:
        * cds_fasta
        * cds_db_connection
        * gi2taxid
        * names
        * nodes
        * ncbitax_db_connection
        '''
        self._h_set_load_type(args)

        if self.cds_source_type == DataAccess.load_type.FILE:
            raise ValueError('Loading CDSs from fasta file currently not supported.')

        self._db_access = DbQuery()

        if self.ncbitax_source_type == DataAccess.load_type.FILE:
            self._h_load_ncbitax_data()

    def get_record(self, version):
        '''
        Returns the record associated with the given accession.version.
        
        :param version: string - GenBank/EMBL/DDBJ/RefSeq Accesion.Version
        :returns: UnityRecord - record associated with the given 
                  accession.version. None if no record is found
        '''
        return self._db_access.get_record(version)

    def get_taxids (self, gis, format=dict):
        '''
        Fetches taxonomy ID for each of the GIs.

        :param gis (list) list of integers representing GIs
        "param format (object type) list, set or dict.
        :rtype based on format parameter, returns either list of 
        tax IDs or a dictionary mapping gis to tax ids. List can 
        contain duplicates.
        '''

        if self.ncbitax_source == DataAccess.load_type.FILE:
            requested_gis = dict((gi, self.gi2taxid.get(gi,-1)) for gi in gis)
            if format == dict:
                return requested_gis
            else:
                return format(requested_gis.values())
        else:
            return self._db_access.get_taxids(gis, format)

    def get_organism_name (self, taxid, name_class='scientific name'):
        if self.ncbitax_source_type == DataAccess.load_type.FILE:
            return self.taxid2name.get(taxid, None)
        else:
            return self._db_access.get_organism_name(taxid, name_class)

    def get_organism_rank (self, query, by_name=False):
        if self.ncbitax_source_type == DataAccess.load_type.FILE:
            return self.taxid2rank.get(taxid, None)
        else:
            return self._db_access.get_organism_rank(query, by_name)


    def _h_load_ncbitax_data(self):
        '''
        Loads taxonomy data from NCBI taxonomy dump files.
        After invoking this method, you can query three dictionaries:
        * gi2taxid (by GI)
        * taxid2name (by taxid)
        * taxid2rank (by taxid)
        '''
        gi2taxid_fpath = self.ncbitax_source['gi2taxid']
        nodes_fpath = self.ncbitax_source['nodes']
        names_fpath = self.ncbitax_source['names']

        self._gi2taxid = loadGi2Taxid(gi2taxid_fpath)
        self._taxid2name = loadNcbiNames(names_fpath)
        self._taxid2rank = loadNcbiRanks(nodes_fpath)

    def _h_set_load_type(self, args):
        '''
        Determines how to load each information:
        * cdss info 
        * ncbitax info
        After invoking this method, following parameters will be set:
        * cds_source_type (enum: load_type)
        * cds_source (str: file_path/db_connection_str)
        * ncbitax_source_type (enum: load_type)
        * ncbitax_source ({ncbitaxfiletype:path_to_path})
        '''
        # determine CDS loading location
        if args.cds_fasta is not None:
            self.cds_source_type = DataAccess.load_type.FILE 
            self.cds_source      = args.cds_fasta
        else:
            assert (args.cds_db_connection is not None)
            self.cds_source_type = DataAccess.load_type.DATABASE
            self.cds_source      = args.cds_db_connection

        # determine ncbitax loading location
        if args.gi2taxid is None or args.nodes is None or args.names is None:
            self.ncbitax_source_type = DataAccess.load_type.DATABASE
            self.ncbitax_source      = args.ncbitax_db_connection
        else:
            assert (args.gi2taxid is not None)
            assert (args.nodes is not None)
            assert (args.names is not None)
            self.ncbitax_source_type = DataAccess.load_type.FILE
            self.ncbitax_source      = {'gi2taxid': args.gi2taxid,
                                        'nodes'   : args.nodes,
                                        'names'   : args.names}