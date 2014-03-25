from testing.solution import Solution
from data.containers.read import ReadContainer
from ncbi.taxonomy.tree import TaxTree

class LCABinner(object):
    ''' Lowest common ancestor implentation of binner.

    Each read is assigned to single taxon in taxonomy tree as the LCA
    of all its alignments.
    '''

    def __init__(self, tax_tree):
        ''' Constructor

        Args:
            tax_tree (TaxTree): Taxonomy tree
        '''
        self.tax_tree = tax_tree

    def bin_reads(self, read_container):
        '''Assigns every read to one or none taxon.

        Uses LCA algorithm to determine assignment.

        Args:
            read_container (ReadContainer): Holds all reads and their alignments
        Returns:
            (Solution): New Solution instance with read assignments.
        '''
        sol = Solution.create_empty()
        for read in read_container.fetch_all_reads():

            # Find LCA
            tax_ids = map(lambda aln: aln.tax_id, read.get_alignments())
            if len(tax_ids) == 0:
                continue
            lca = self.tax_tree.find_lca(tax_ids)

            # Store assignment
            sol.add_assignment(read.id, lca)

        return sol
