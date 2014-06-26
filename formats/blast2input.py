


class BLASTParser (object):
    ''' Enables BLAST to input format parsing
    '''
    def __init__ (self, output_format = None):
        ''' @param output_format: BLAST output format. 
            If none provided, standard output is taken.
            Standard output format is:
            qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        '''
        if not output_format:
            output_format = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        self._set_up_parser(output_format)
            
            
    
    ############ ############ ##############
    def _set_up_parser (self, output_format):
        ''' Takes output format specification and sets up
            the parser for the specific format parsing
        '''
        self.fmt_values = {}
        values = output_format.strip().split()
        i = 0
        for value in values:
            self.fmt_values[value] = i
            i += 1
         
    ############ ############ ##############
    def convert_file (self, blast_output_fname, output_fname):
        blast_output_file = open(blast_output_fname)
        output_file       = open(output_fname, 'w')
        readline = blast_output_file.readline
        
        read_id     = None
        alignment_list = []

        while (True):
            line = readline()
            line = line.strip()
            if not line:
                break
           
            if line.startswith('#'):
                continue 
            (new_read_id, aln_data) = self.parse_line(line)
            if not read_id:
                read_id = new_read_id

            if new_read_id != read_id:
                print >> output_file, self.get_input_line(read_id, alignment_list)
                read_id = new_read_id
                alignment_list = [aln_data]
            else:
                alignment_list.append(aln_data)

        print >> output_file, self.get_input_line(read_id, alignment_list)            
        blast_output_file.close()
    
    def get_input_line (self, read_id, alignment_data):
        output_line = "{0},{1};".format(read_id, len(alignment_data))
        for alignment in alignment_data:
            output_line += str(alignment) + ";"
        return output_line
    
    def parse_line (self, line): 
        aln_data_list   = line.split()
        subject_id      = aln_data_list[self.fmt_values['sseqid']]
        (gi, db_source, nucl_accession) = subject_id.split('|')[1:4]
        
        score       = aln_data_list[self.fmt_values['bitscore']]
        start       = aln_data_list[self.fmt_values['sstart']]
        stop        = aln_data_list[self.fmt_values['send']]
        query_id    = aln_data_list[self.fmt_values['qseqid']]
        
        aln_data    = AlignmentData(nucl_accession, db_source, gi, 
                                    score, start, stop)
        return (query_id, aln_data)
    

class AlignmentData (object):
    
    def __init__ (self, nucleotide_accession, db_source, gi, score, start, stop):
        self.nucleotide_accession   = nucleotide_accession
        self.db_source              = db_source
        self.gi                     = int (gi)
        self.score                  = float (score)
        self.start                  = int(start)
        self.stop                   = int(stop)
        if start > stop:
            self.strand = '+'
            tmp = self.start
            self.start = self.stop
            self.stop = tmp
        else:
            self.strand = '-'
            
    def __str__(self, *args, **kwargs):
        return "{0},{1},{2},{3},{4},{5},{6}".format (self.nucleotide_accession,
                                                     self.db_source, self.gi,
                                                     self.score, self.start, 
                                                     self.stop, self.strand)
            
        
def main():
    argparser = argparse.ArgumentParser() 
    argparser.add_argument('BlastAlnFile', help='Blast alignment file generated using outfmt 6-10', type=str)
    argparser.add_argument('BinnerAlnFile', help='Binner input alignment file. Format specified in Binner Docs.')
    argparser.add_argument('-f', '--format', help='Blast format specification', type=str,
                           default='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore')

    args = argparser.parse_args()
    parser = BLASTParser(args.format)
    parser.convert_file(args.BlastAlnFile, args.BinnerAlnFile)



if __name__ == '__main__':
    import sys, argparse
    main()
