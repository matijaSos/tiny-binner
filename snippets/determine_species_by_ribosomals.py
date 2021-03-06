
import argparse
import sys,os
sys.path.append(os.getcwd())

from utils.argparser import DefaultBinnerArgParser
from ncbi.db.data_access import DataAccess
from ncbi.taxonomy.tree import TaxTree

# Chart stuff
from pylab import *

class ArgParser(DefaultBinnerArgParser):
    def __init__(self):
        super(ArgParser, self).__init__('''\
                Takes file with ribosomal CDSs data, \
                and estimates present organisms. \
                .''')
        self.add_argument('cds_data_file',
                help='file with ribosomal CDSs',
                type=str)
        self.add_argument('img_export_path',
                help='where to store created pie chart',
                type=str)

def main():

    # Input arguments
    argparser   = ArgParser()
    args        = argparser.parse_args()
    # Access database
    dataAccess = DataAccess(args)

    print '1. Loading tax tree...'
    tax_tree = TaxTree()
    print 'done.'

    # Get file contents
    lines = [line.rstrip('\n') for line in open(args.cds_data_file)]
    # Skip header
    # lines = lines[6:]
    
    # Dictionary with tax_ids
    species = {}
    for line in lines:
        row_data = line.split()

        tax_id = int(row_data[2])
        tax_id_species = tax_tree.get_parent_with_rank(tax_id, 'species')

        # Add to dict
        species[tax_id_species] = species.get(tax_id_species, 0) + 1

    # Remove species with not enough CDSs
    minCDSNum = 10 # Set as parameter
    keysToPop = []
    for key, val in species.iteritems():
        if val < minCDSNum: keysToPop.append(key)
    for key in keysToPop:
        species.pop(key)
        
    print
    print("Number of species: {0}".format(len(species)))
    print

    # Total number of CDSs
    total = sum(species.values())

    # Print species - sorted
    labels  = []
    fracs   = []
    explode = []
    for tax_id in sorted(species, key=species.get, reverse=True):
        num     = species[tax_id]
        frac    = float(num)/total * 100
        try:
            name = tax_tree.nodes[int(tax_id)].organism_name
        except:
            name = "unknown"
        print("{0:10} {1:80} {2:8d} {3:10.4f}%".format(tax_id, name, num, frac))

        # Chart data
        label = name.replace(" ", "\n")
        if num < 40:
            label = ''

        labels.append(label)
        fracs.append(frac)
        explode.append(0)

    # Create and save pie chart
    figure(1, figsize=(12,12))
    ax = axes([0.10, 0.10, 0.80, 0.80])
    explode[0] = 0.1

    pie(fracs, explode=explode, labels=labels, 
        autopct='%1.1f%%', shadow=False, startangle=90)
    title('Microbiome: ' + str(total) + " CDSs")

    savefig(args.img_export_path, format='png')


if __name__ == '__main__':
    main()
