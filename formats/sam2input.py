import pysam
import os, sys
sys.path.append(os.getcwd())

class SamParser(object):
    def convert_file(self, sam_input_file, output_file):
        samfile = pysam.Samfile(sam_input_file)
        for read in samfile.fetch():
            if read.is_paired:
                print read

def main():
    if len(sys.argv) < 3:
        print 'Usage:\npython sam2input.py <INPUT SAM FILE> <BINNER INPUT FILE>'
        sys.exit(-1)
    samParser = SamParser()
    samParser.convert_file(sys.argv[1], sys.argv[2])


if __name__=='__main__':
    main()

