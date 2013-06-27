
def convertLisaToInput (lisa_alignment_fname, output_fname):
    input_fhandle = open (lisa_alignment_fname, 'r')
    output_fhandle = open (output_fname, 'w')

    while (True):
        line = input_fhandle.readline()
        if not line:
            break
        line = line.strip()
        alns = line.split(';')

        # write header
        output_fhandle.write("%s" % alns[0])
        # if only header present
        if len(alns) == 1:
            output_fhandle.write("\n")
            continue

        for aln_str in alns[1:]:
            try:
                (nucl_data,score,start,stop,strand) = aln_str.split(',')
            except ValueError:
                continue
            nucl_data_list = nucl_data.split('|')
            gi          = nucl_data_list[1]
            db_source   = nucl_data_list[2]
            nucl_acc    = nucl_data_list[3]

            output_strand = '+' if strand == '0' else '-'
            output_aln_str = "{0},{1},{2},{3},{4},{5},{6}".format(nucl_acc, db_source, gi, score, start, stop, output_strand)

            output_fhandle.write(";%s" % output_aln_str)
        output_fhandle.write('\n')

    input_fhandle.close()
    output_fhandle.close()


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 3:
        print "Usage:\npython lisa2input.py <LISA_ALN_FILE> <OUTPUT_FILE_NAME>"
        sys.exit(-1)
    convertLisaToInput (sys.argv[1], sys.argv[2])



