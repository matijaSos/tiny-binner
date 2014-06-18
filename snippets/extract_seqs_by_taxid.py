
import argparse
import sys,os
import time
sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree

from Bio import SeqIO


class ArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(ArgParser, self).__init__('''\
                Loads file with taxids, \
                fasta file \
                and extracts all sequences that match to given taxids.
                .''')
        self.add_argument('taxids_file',
                help='path to the file with taxids',
                type=str)
        self.add_argument('fasta_in',
                help='path to the fasta file',
                type=str)
        self.add_argument('fasta_out',
                help='path where filtered fasta will be saved',
                type=str)


def main():

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()

    # Access database
    dataAccess = DataAccess(args)

    # ------------------ #

    print '1. Loading tax tree...'
    start = time.time()

    tax_tree = TaxTree()

    end = time.time()
    print("done: {0:.2f} sec".format(end - start))

    # ------------------ #

    # Get taxid file contents
    lines = [line.rstrip('\n') for line in open(args.taxids_file)]
    tax_ids = set()
    for line in lines:
        row_data = line.split()
        tax_id = row_data[1]
        # Add to set
        tax_ids.add(int(tax_id))

        print("Added: " + tax_id)

    # Extract records
    extracted_records = []
    with open(args.fasta_in, "rU") as fasta_in:
        for record in SeqIO.parse(fasta_in, "fasta"):

            # Get tax_id
            gi = record.id.split('|')[1]
            tax_id = dataAccess.get_taxids([gi], format=list)[0]

            if tax_id in tax_ids:
                extracted_records.append(record)
                print record.id

    SeqIO.write(extracted_records, args.fasta_out, "fasta")

if __name__ == '__main__':
    main()
