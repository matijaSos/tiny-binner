import pysam
import os, sys
sys.path.append(os.getcwd())

class SamParser(object):
    def convert_file(self, sam_input_fname, output_fname):

        #----- OPEN INPUT SAM FILE ----#
        sam_file = pysam.Samfile(sam_input_fname)
        #------  OPEN OUTPUT FILE -----#
        output_file = open(output_fname, 'w')

        #---- PARSE ----#
        last_read_id = None
        reads_to_process = []
        for readAlignment in sam_file.fetch():
            # handle initial case
            if last_read_id is None:
                reads_to_process.append(readAlignment)
                last_read_id = readAlignment.qname
                continue
            # different alignment of the same read
            if last_read_id == readAlignment.qname:
                reads_to_process.append(readAlignment)
                continue
            # next read
            else:
                self._process_reads(reads_to_process, sam_file, output_file)
                reads_to_process = [readAlignment]
                last_read_id = readAlignment.qname
        self._process_reads(reads_to_process)

        output_file.close()

    def _process_reads(self, reads_to_process, sam_file, output_file):
        read_id = freads_to_process[0].qname
        read_str = '{0},{1};'.format(read_id, len(reads_to_process))

        for read in reads_to_process:
            assert (read_id == read.qname)
            read_str += self._format_str(read, sam_file)
            read_str += ';'

        output_file.write('%s\n' % read_str)


    def _format_str(self, sam_aligned_read, sam_file):
        (x,gi,db_source,accession,x) = sam_file.getrname(sam_aligned_read.rname).split('|')
        score = sam_aligned_read.mapq
        start = sam_aligned_read.aend - sam_aligned_read.alen
        stop  = sam_aligned_read.aend
        strand = '-' if sam_aligned_read.is_reverse else '+'

        return '{0},{1},{2},{3},{4},{5},{6}'.format(accession, db_source, gi, score, start, stop, strand)


def main():
    if len(sys.argv) < 3:
        print 'Usage:\npython sam2input.py <INPUT SAM FILE> <BINNER INPUT FILE>'
        sys.exit(-1)
    samParser = SamParser()
    samParser.convert_file(sys.argv[1], sys.argv[2])


if __name__=='__main__':
    main()

