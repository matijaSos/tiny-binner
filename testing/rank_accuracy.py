
from testing.solution import Solution
from data.containers.read import ReadContainer
from ncbi.taxonomy.tree import TaxTree

class RankAccuracy(object):
    '''
    Compares correct with given solution.
    Validity of each read's assignment (some taxon) is
    evaluated through all ranks.
    Read can be to each rank assigned either correct or incorrect.
    '''
    def __str__(self):
        return str(len(self.data))

    def printData(self):
        # Header
        print "Rank\t\tTrue\tFalse"
        print "========================="

        # Data
        for rank in reversed(self.ranks):
            ts = self.data[rank]['True']
            fs = self.data[rank]['False']

            tab = "" if rank == "superkingdom" else "\t"
            print rank + tab + "\t" + str(ts) + "\t" + str(fs)

    # Ranks for which accuracy is checked
    ranks = ['species', 
             'genus', 
             'family', 
             'order', 
             'class', 
             'phylum', 
             'superkingdom']

    def __init_data(self):
        self.data = {}
        for rank in self.ranks:
            self.data[rank] = {'True': 0, 'False': 0}

    def __init__(self, tax_tree, correct_sol, read_container):
        '''
        @param tax_tree         (TaxTree)
        @param correct_sol      (Solution)
        @param read_container   (ReadContainer)
        '''
        # Initialize dictionary holding results
        self.__init_data()

        # Get those data
        for read in read_container.fetch_all_reads():
            
            # Get tax_id of LCA
            tax_ids = map(lambda aln: aln.tax_id, read.get_alignments())
            if len(tax_ids) == 0:
                continue
            lca = tax_tree.find_lca(tax_ids)
            
            # Now evaluate read assignment through all ranks
            correct = correct_sol.get_tax_id(read.id)
            matches = False
            for rank in self.ranks:
                if matches:
                    self.data[rank]['True'] += 1
                    continue
                    
                lca_parent = tax_tree.get_parent_with_rank(lca, rank)
                correct_parent = tax_tree.get_parent_with_rank(correct, rank)

                # If taxon is not specific enough
                if lca_parent == 0:
                    continue

                lca_name = tax_tree.nodes[lca_parent]
                correct_name = tax_tree.nodes[correct_parent]

                if lca_parent == correct_parent:
                    self.data[rank]['True'] += 1
                    matches = True
                else:
                    self.data[rank]['False'] += 1
                

