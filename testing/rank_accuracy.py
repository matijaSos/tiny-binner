
from testing.solution import Solution
from data.containers.read import ReadContainer
from ncbi.taxonomy.tree import TaxTree

class RankAccuracy(object):
    ''' Compares correct with given solution.

    Validity of each read's assignment (some taxon) is
    evaluated through all ranks.
    Read can be to each rank assigned either correct or incorrect.
    '''
    # Ranks for which accuracy is checked
    ranks = ['species', 
             'genus', 
             'family', 
             'order', 
             'class', 
             'phylum', 
             'superkingdom']

    @staticmethod
    def print_comparison(rank_accs):
        ''' Prints multiple RankAccuracy objects in such
            manner to be easily comparable.

        Args:
            ranks_accs ([RankAccuracy]): list of RankAcc objects
        '''
        # Header
        header = "Rank\t"
        for i in range(len(rank_accs)):
            header += "\tTrue\tFalse"
        header += "\n" + 20*"=" + 15*len(rank_accs)*"="
        print header

        for rank in reversed(RankAccuracy.ranks):
            tab = "" if rank == "superkingdom" else "\t"
            line = rank + tab

            for rank_acc in rank_accs:
                ts = rank_acc.data[rank]['True']
                fs = rank_acc.data[rank]['False']
                line += "\t" + str(ts) + "\t" + str(fs)
            print line

    def __str__(self):
        return str(len(self.data))

    def print_data(self):
        # Header
        print "Rank\t\tTrue\tFalse"
        print "==============================="

        # Data
        for rank in reversed(self.ranks):
            ts = self.data[rank]['True']
            fs = self.data[rank]['False']

            tab = "" if rank == "superkingdom" else "\t"
            print rank + tab + "\t" + str(ts) + "\t" + str(fs)


    def __init_data(self):
        self.data = {}
        for rank in self.ranks:
            self.data[rank] = {'True': 0, 'False': 0}

    def __init__(self, tax_tree, correct_sol, binner_sol):
        '''Constructor
        
        Evaluates accuracy of given read assignment to the correct one.

        Args:
            tax_tree (TaxTree): Taxonomy tree
            correct_sol (Solution): Correct read assignment
            binner_sol (Solution): Given read assignment to be evaluated
        '''
        # Initialize dictionary holding results
        self.__init_data()

        # Examine every read
        for read_id, tax_id in binner_sol.id2taxon.iteritems():
            
            # Now evaluate read assignment through all ranks
            correct = correct_sol.get_tax_id(read_id)
            matches = False
            for rank in self.ranks:
                if matches:
                    self.data[rank]['True'] += 1 # If correct on lower, is also on higher level
                    continue
                    
                try:
                    parent = tax_tree.get_parent_with_rank(tax_id, rank)
                except KeyError:
                    # TODO: Do some logging for stats
                    continue
                correct_parent = tax_tree.get_parent_with_rank(correct, rank)

                # If taxon is not specific enough
                if parent == 0:
                    continue

                if parent == correct_parent:
                    self.data[rank]['True'] += 1
                    matches = True
                else:
                    # Propagate False to lower levels
                    for rank2 in self.ranks:
                        self.data[rank2]['False'] += 1
                        if rank2 == rank:
                            break

                # Go up in tree
                tax_id = parent
                correct = correct_parent
