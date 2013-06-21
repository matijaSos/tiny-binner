
from data.alignment import ReadAlnLocation
import logging
from utils.autoslots import Autoslots

log = logging.getLogger(__name__)

class Read (Autoslots):
    """ Contains all the read-related information,
        such as the its identifier and the list of alignment
        locations
    """
    def __init__ (self, read_id, read_length, alignment_locations):
        self.id                     = read_id
        self.length                 = read_length
        self.alignment_locations    = alignment_locations  # Jel ovo [ReadAlnLocation]? Treba iskomentirati!!!
        self.potential_host         = None


    @staticmethod
    def from_read_str (read_str):
        """ Parses the description string and creates a new read from it with
            accompanying alignment locations
        """
        # Attributes for creating new read
        newRead_id       = None
        newRead_length   = None # Not available for now, should be in the future
        newRead_aln_locs = []

        # valList: ["read_id, num_align", "alignInfo1", "alignInfo2", ... "alignInfoN"]
        read_str = read_str.strip()
        if read_str.endswith(';'):
            read_str = read_str[0:-1]
        valList = read_str.split(';');

        # Get header info
        # headerList: [read_id, num_align]
        headerList = valList[0].split(',');
        newRead_id = headerList[0];
        if newRead_id.startswith('@'):
            newRead_id = newRead_id[1:]
        num_align = headerList[1];

        # Store every alignInfo
        for alignInfo in valList[1:]:
            data = alignInfo.split(','); # [nucl_acc, db_source, GI, score, start, stop, strand]

            nucl_acc    = data[0]
            db_source   = data[1]
            GI          = int   (data[2])
            score       = float (data[3])
            start       = int   (data[4])
            stop        = int   (data[5])
            strand      = data[6]
            assert (strand in ['+', '-'])
            complement  = False if strand=='+' else True

            # Create and store new ReadAlnLocation object
            try:
                newAlignInfo = ReadAlnLocation(newRead_id, nucl_acc, db_source, GI, score,
                                           (start, stop),   complement);
                newRead_aln_locs.append(newAlignInfo)
            except Exception as e:
                log.error("Location parsing error.", exc_info=True)

        return Read(newRead_id, newRead_length, newRead_aln_locs)


    def is_host(self):
        '''
        Returns host status of read

        :rtype boolean: True if read is potential host read,
        False if not, and None if still undecided
        '''
        return self.potential_host

    def get_alignments (self, format=list):
        '''
        Get read alignments for the read.
        @param: format (collection or iterator) format in which to
        acquire the alignments
        '''
        assert (format in [iter, list, set])
        return format(self.alignment_locations)

    def set_alignments (self, alignments):
        self.alignment_locations = alignments

    def has_alignments (self):
        return len(self.alignment_locations) > 0
