import argparse
import os

class DefaultBinnerArgParser (argparse.ArgumentParser):
    def __init__(self, description):
        self.formatter_class = argparse.ArgumentDefaultsHelpFormatter
        self.description = description

        self._set_default_binner_values()

    def _set_default_binner_values(self):
        self.add_argument('input', help='Input alignment file', 
                           type=str)
        self.add_argument('descr', help='XML description schema', 
                               type=str)
        self.add_argument('output', help='Output XML file', 
                               type=str)
        self.add_argument('-l', '--log_configuration',
                               help='Logging configuration file', type=str,
                               default='config' + os.path.sep + 'logging.ini')
        mutexgroup_cds = self.add_mutually_exclusive_group()
        mutexgroup_cds.add_argument('--cds-db-connection', 
            default='mysql+mysqldb://root:root@localhost/unity',
            help='CDS database connection string')
        mutexgroup_cds.add_argument('--cds-fasta',
            help='CDS fasta file location')

        mutexgroup_tax = self.add_mutually_exclusive_group()
        mutexgroup_tax.add_argument('--ncbitax-db-connection', 
           default='mysql+mysqldb://root:root@localhost/ncbitax',
            help='NCBI Taxonomy database connection string')
        ncbi_tax_files = mutexgroup_tax.add_argument_group()
        ncbi_tax_files.add_argument('--gi2taxid',
            help='NCBI Taxonomy gi2taxid dump file')
        ncbi_tax_files.add_argument('--nodes',
            help='NCBI Taxonomy nodes dump')
        ncbi_tax_files.add_argument('--names',
            help='NCBI Taxonomy names dump')
        self.add_argument('-tt', '--tax-tree', 
           help='Taxonomy tree location', 
           default='./ncbi/taxonomy/.data/ncbi_tax_tree')        


def validate_args(args):
    error = False
    if not os.path.exists(os.path.expanduser(args.input)):
        print "Input alignment file %s doesn't exist" % args.input
        error = True
    if not os.path.exists(os.path.expanduser(args.descr)):
        print "XML description schema %s doesn't exist" % args.descr
        error = True
    if not os.path.exists(os.path.expanduser(args.log_configuration)):
        print "Log configuration file %s doesn't exist" % args.log_configuration
        error = True
    if args.cds_fasta is not None:
        if not os.path.exists(os.path.expanduser(args.cds_fasta)):
            print "CDS Fasta file %s doesn't exists" % args.descrargs.cds_fasta
    if args.gi2taxid is not None:
        if not os.path.exists(os.path.expanduser(args.gi2taxid)):
            print "gi_taxid_[nucl/prot] file %s doesn't exist." % args.gi2taxid
            error = True
    if args.nodes is not None:
        if not os.path.exists(os.path.expanduser(args.nodes)):
            print "nodes dump file %s doesn't exist." % args.nodes
            error = True
    if args.names is not None:
        if not os.path.exists(os.path.expanduser(args.names)):
            print "names dump file %s doesn't exist." % args.names
            error = True
    return error