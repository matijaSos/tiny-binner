
import logging
from ncbi.db.access import WrongTableError, DbQuery

class RecordContainer (object):
    ''' Serves as a local Record Repository.
        If a GenBank/EMBL/DDBJ record has already been
        fetched from the database, it can be fetched localy
        from the record repository.
    '''

    def __init__ (self):
        self.record_repository  = {}
        self.num_missing_records = 0
        self.log = logging.getLogger(__name__)

    def get_num_missing_records_stats(self):
        missed_recs = self.num_missing_records
        all_recs = len(self.record_repository)
        return {'num_missing_records':missed_recs,
                'num_all_records':all_recs,
                'missing_records_percentage':
                '{0:.2f}'.format(missed_recs/float(all_recs)*100)+'%'}

    def set_db_access(self, db_query):
        '''
        @param: db_query (DbQuery, MockDbQuery)
        '''
        assert (hasattr(db_query, 'get_record'))
        self.db_query = db_query

    def populate (self, versions, table='cds'):
        '''
        Populates the record container with all the records
        that have produced significant alignments

        :param list of NT (GenBank, EMBL, DDBJ) accession.versions
        :param table (str) available tables are cds, rrna, mrna, misc_rna
        '''
        if table not in DbQuery.supported_tables:
            raise WrongTableError(table)
        for version in versions:
            self.fetch_record(version)

    def fetch_record (self, nucleotide_accession):
        '''
        @param nucleotide_accession (str)
        @return Record (ncbi/db/[genbank/embl])
        '''
        self._add_record(nucleotide_accession)
        return self.record_repository[nucleotide_accession]

    def fetch_existing_record (self, nucleotide_accession):
        '''
        @param nucleotide_accession (str)
        @return UnityRecord (ncbi/db/[genbank/embl])
        '''
        return self.record_repository.get(nucleotide_accession)

    def fetch_all_records (self, format=iter):
        '''
        Fetches all loaded records in a specified format
        @param: format (iter, list, set)
        @return format(records)
        '''
        assert (format in [iter, list, set])
        return format(self.record_repository.values())

    def _add_record (self, record_id):
        ''' Adds the record from database if not already present
	   If unable to find entry in database, stores None instead.
        '''
        try:
            getattr(self, 'db_query')
        except AttributeError:
            raise AttributeError("RecordContainer has not attribute 'db_query'. Did you forget to envoke set_db_access()?")
        if not self.record_repository.has_key(record_id):
            record = self.db_query.get_record(record_id) # What is type of this object?
            try :
                getattr(record, 'version')
                self.record_repository[record_id] = record
            except AttributeError:
                self.log.info("No record with ID %s", str(record_id))
                self.record_repository[record_id] = None
                self.num_missing_records += 1
