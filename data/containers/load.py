from data.containers.read   import ReadContainer
from data.containers.record import RecordContainer
from data.containers.cdsaln import CdsAlnContainer
from ncbi.db.mock_db_access import MockDbQuery

def fill_containers (alignment_file, db_access):
    '''
    Populates read, record and CDS alignment container.
    @return tuple(ReadContainer, RecordContainer, CdsAlnContainer)
    '''

    read_cont   = ReadContainer()
    record_cont = RecordContainer()
    record_cont.set_db_access(db_access)
    cdsaln_cont = CdsAlnContainer()

#   1. Load all the information available in the alignment file
    read_cont.populate_from_aln_file(alignment_file)
#   2. Fetch all the records reported in the alignment file from the database
    record_cont.populate(read_cont.fetch_all_reads_versions())
#   3. Find to which coding sequences reads map
    read_cont.populate_cdss(record_cont)
#   4. Populate Cds Alignment container
    cdsaln_cont.populate(read_cont.fetch_all_reads())

    return (read_cont, record_cont, cdsaln_cont)
