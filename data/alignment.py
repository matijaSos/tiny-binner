from data.containers.record import RecordContainer
from utils.location import Location
from utils.location import LoactionParsingException

class ReadAlnLocation (object):
    """ Contains information on alignment location on 
        an NT nucleotide string
    """
    
    def __init__ (self, read_id, nucleotide_accession, db_source, genome_index, score, 
                  location_span, complement, active=True):
        self.read_id                = read_id
        self.nucleotide_accession   = nucleotide_accession
        self.db_source              = db_source
        self.genome_index           = genome_index
        self.score                  = score
        self.location_span          = location_span
        self.complement             = complement
        self.active                 = active
        self.tax_id                 = None
        # self.determine_coding_seqs()
        # Sto je sa .aligned_cdss? Navesti to negdje u komentarima ako postoji!
    
    def set_active (self, active):
        '''
        Sets active status for the read alignment.
        Inactive reads do not go into CDS alignments.
        '''
        self.active = active

    def set_potential_host_status (self, potential_host):
        '''
        Set to true if organism is potential host [child of 
        animalia kingdom]
        @param potential_host (boolean) 
        '''
        self.potential_host = potential_host

    def is_potential_host (self):
        """ Returns true if organism is potential host 
        (child of animalia kingdom), false otherwise.
        @return (boolean)
        """
        return self.potential_host

    def determine_coding_seqs (self, record_container):
        ''' Determines which of the CDSs in the record aligned_regions
            aligned to the read.
            @return list of tuples (cds, intersecting_location) if such exist, 
            None if record is not available from the database
        '''
        self.aligned_cdss = []
        record = record_container.fetch_record (self.nucleotide_accession)

        # if not possible to fetch a record from the db, return None
        if not record:
            return None

        (start,stop) = self.location_span
        try:
            location = Location.from_location_str("%d..%d" % (start, stop))
        except LoactionParsingException, e:
            print "ReadAlignment/determine_coding_seqs:", e
            self.aligned_cdss = []
            return self.aligned_cdss
        
        for cds in record.cds:
            try:
                cds_location = Location.from_location_str(cds.location)
            except LoactionParsingException, e: 
                print "ReadAlignment/determine_coding_seqs:", e
                continue
            location_intersection = cds_location.find_intersection (location)
            if location_intersection is not None:
                self.aligned_cdss.append ((cds, location_intersection))
        
        return self.aligned_cdss

    # -------------------------------- Detrmine CDSs - Optimal ---------------------------- #

    def __overlap (self, cds_loc, aln_loc):
        ''' Returns True if given CDS and alignment overlap.
            Only (start, end) is taken in consideration (no sublocations)

            @param   (Location)  cds_loc CDS Location
            @param   (Location)  aln_loc Alignment Location
            @return  (Boolean)   True if CDS overlaps with alignment
        '''
        return not (cds_loc.end < aln_loc.start or cds_loc.start > aln_loc.end)

    # ---------- #

    def __get_cds_rel_pos (self, cdss, cds_id, aln_loc): 
        '''Returns position of CDS relative to alignment.

            @param  [Cds]       List of CDSs
            @param  (Location)  cds_id CDS we are looking at
            @param  (Location)  aln_loc Alignment Location

            @return (string)    LEFT_OF_ALN  - fully left of alignment      
                                RIGHT_OF_ALN - fully right or not first which overlaps 
                                FIRST        - first to overlap
        '''
        cds_loc = Location.from_location_str(cdss[cds_id].location)

        if (cds_loc.end < aln_loc.start):
            return "LEFT_OF_ALN"

        if (cds_loc.start > aln_loc.end):
            return "RIGHT_OF_ALN"

        # Overlap occured - check if it is first
            # If it is first CDS or the previous one does not overlap
        if (cds_id == 0):
            return "FIRST"

        cds_prev_loc = Location.from_location_str(cdss[cds_id - 1].location)
        if (not self.__overlap(cds_prev_loc, aln_loc)):
            return "FIRST"

        return "RIGHT_OF_ALN"

    # ---------- #

    def __find_first_overlapping_CDS_id (self, aln_location, cdss):
        ''' Find id of the first CDS which overlaps with the given alignment.
            Uses binary search algorithm.

            @param   (Location) aln_location    Alignment Location
            @param   [Cds]      cdss            List of CDSs, sorted by start
            @returns (int|None)                 Id in cdss of described CDS, None if no overlap
        '''
        lo = 0
        hi = len(cdss) - 1

        # If cdss is empty -> return None
        if (hi < 0):
            return None

        while (lo < hi):
            mid = lo + (hi - lo) // 2   # '//' for python 3 compatibility

            cds_rel_pos = self.__get_cds_rel_pos (cdss, mid, aln_location)

            if (cds_rel_pos == "LEFT_OF_ALN"):
                lo = mid + 1
                
            if (cds_rel_pos == "RIGHT_OF_ALN"):
                hi = mid - 1

            if (cds_rel_pos == "FIRST"):
                return mid

        # Check lo
        cds_location = Location.from_location_str(cdss[lo].location)
        if self. __overlap(cds_location, aln_location):
            return lo
        else:
            return None

    # ---------- #

    def determine_coding_seqs_optimal (self, record):
        ''' Determines which of the CDSs in the record aligned_regions
            aligned to the read.

            @param (UnityRecord) record Record that is used
            @return                  list of tuples (cds, intersecting_location) if such exist, 
                                     None if record is not available from the database
        '''
        
        self.aligned_cdss = []

        # If not possible to fetch a record from the db, return None
        if not record:
            return None

        # Acquire alignment Location
        (start, stop) = self.location_span
        try:
            aln_location = Location.from_location_str("%d..%d" % (start, stop))
        except LoactionParsingException, e:
            print "ReadAlignment/determine_coding_seqs:", e
            self.aligned_cdss = []
            return self.aligned_cdss

        # Determine first overlapping CDS - binary search
        first_ovp_id = self.__find_first_overlapping_CDS_id (aln_location, record.cds) 

        # No CDS from the list overlaps - return []
        if (first_ovp_id == None):
            return self.aligned_cdss

        # Determine following overlapping CDSs - loop while overlaps
        for i in range(first_ovp_id, len(record.cds)):
            cds = record.cds[i]
            cds_location = Location.from_location_str(cds.location)

            # If this one does not overlap, the others also won't because it's sorted
            if not self.__overlap(cds_location, aln_location):
                break

            location_intersection = cds_location.find_intersection (aln_location)
            if location_intersection is not None:
                self.aligned_cdss.append ((cds, location_intersection))

        return self.aligned_cdss

    # ---------------------------------------------------------------------------- #

    def set_type (self):
        """ Location can be coding or non-coding
        """
        pass


class CdsAlignment (object):
    ''' Contains all the alignment information for a single
        CDS, meaning: 
        - all the reads mapped to that CDS and
        - their corresponding sublocations
    '''
    
    def __init__ (self, cds):
        self.cds = cds              # CDS object (from Mladen)
        self.aligned_regions = {}   # dictionary of CdsAlnSublocation
        pass
    
    def __hash__ (self):
        return hash ((self.cds.record_id, self.cds.location))

    def __eq__ (self, other):
        if (other == None): return False
        return (self.cds.record_id, self.cds.location) == (other.cds.record_id, other.cds.location)

    def add_aligned_sublocation (self, read_id, aligned_location, score):
        ''' Adds an aligned region to the cds unless it comes 
            from the read already present in the aligned regions
            @param read_id read id 
            @param (Location) aligned_location intersection between the CDS and the 
                    alignment location
            @param score alignment score for this read
        '''
        
        # if the CDS has already been covered by the same read in the past,
        # discard this one. 
        if self.aligned_regions.has_key(read_id):
            return
        
        aligned_sublocation             = CdsAlnSublocation (read_id, aligned_location, score)
        self.aligned_regions[read_id]   = aligned_sublocation

    def is_active(self):
        '''
        Checks whether CDS alignment is active. 
        If all the CDSAlnSublocations are inactive, then the whole CdsAlignment
        is inactive.
        '''
        for cds_aln_subloc in self.aligned_regions.values():
            if cds_aln_subloc.active:
                return True
        return False

    def get_active_alignment_cnt (self):
        '''
        Counts the number of CdsAlnSublocations which
        have been marked as active
        @return (int) number of active alignments
        '''
        active_sublocations = 0
        for cds_aln_subloc in self.aligned_regions.values():
            if cds_aln_subloc.active:
                active_sublocations += 1
        return active_sublocations

    def get_key (self):
        pass
        
    def contains_read (self, read_id):
        ''' Determines whether this CDS alignment contains 
            a subalignment mapped to the specified read.
        '''
        return True if self.aligned_regions.has_key(read_id) else False
        
    def __str__(self):
        tab = " " * 2
        ret = "CdsAlignment\n"
        ret += tab + "cds: " + str(self.cds) + "\n"
        ret += tab + "aligned_regions:\n"
        for (key, aln_reg) in self.aligned_regions.items():
            ret += tab*2 + "(key) " + key + ":\n"
            ret += tab*3 + str(aln_reg).replace("\n", "\n"+(tab*3)) + "\n"
        return ret

    
class CdsAlnSublocation (object):
    ''' Represents the sublocation of a CDS covered
        by a single read.
    '''
        
    def __init__ (self, read_id, location, score, active=True):
        ''' @param read_id read ID 
            @param (Location) location intersection location 
            @param score alignment score
            @param (boolean) active If active than it maps to CDS that contains it.
        '''
        self.read_id    = read_id
        self.location   = location
        self.score      = score
        self.active     = active

    def __str__(self):
        tab = " "*2
        ret = "CdsAlnSublocation\n"
        ret += tab + "read_id:  " + str(self.read_id) + "\n"
        ret += tab + "location: " + str(self.location) + "\n"
        ret += tab + "score:    " + str(self.score) + "\n"
        ret += tab + "active:   " + str(self.active)
        return ret

