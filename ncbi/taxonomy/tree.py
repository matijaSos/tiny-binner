from  collections       import defaultdict
import os,sys
sys.path.append(os.getcwd())
from utils.progressbar import print_progress

class TaxTree ():
    ''' Loads the NCBI taxonomy tree, creates both
        parent-child and child-parent relations,
        enables parent-child relationship testing and
        finding the least common ancestor.
    '''

    def __init__ (self, parent2child_fname=None, tax_nodes_fname=None):
        ''' Locates the ncbi taxonomy file and sets the important
            taxonomy assignments (such as animalia, bacteria ecc)

            :param parent2child_fname location of the ncbi taxonomy tree file
            :param tax_nodes_fname location of the file containing taxid,
            organism name and organism rank for each taxid in the tree.
        '''

        if not parent2child_fname:
            parent2child_fname = self._h_find_taxnode_file('parent2child')
        self.load(parent2child_fname)
        if not tax_nodes_fname:
            tax_nodes_fname = self._h_find_taxnode_file('taxdata')
        self.load_taxonomy_data(tax_nodes_fname)

        #--------- RELEVANT TAXONOMY ASSIGNMENTS ----------#
        self._h_set_relevant_taxonomy_assignments()
        self._h_map_taxids_to_relevant_tax_nodes()

    def load (self, parent2child_fname):
        self.parent_nodes   = self._h_get_tax_nodes(parent2child_fname)
        self.child_nodes    = self._h_populate_child_nodes()

    def load_taxonomy_data(self, tax_nodes_fname):
        '''
        Uses data access object to find organism name and
        rank of each of the tax IDs.
        For each tax ID creates a node of type TaxNode
        After invoking this method, there is nodes parameter
        of type dict(key=tax_id:int, value=node:TaxNode)
        '''
        self.nodes = {}
        total = len(self.parent_nodes)
        current = 0
        tax_nodes_file = open(tax_nodes_fname, 'r')
        readline = tax_nodes_file.readline
        while (True):
            line = readline()
            if not line: break
            (taxid, org_name, rank) = line.strip().split('|')
            node = TaxNode(org_name, rank)
            self.nodes[int(taxid)] = node
        tax_nodes_file.close()

    def is_child (self, child_taxid, parent_taxid):
        ''' Test if child_taxid is child node of parent_taxid
            Node is not the child of itself
        '''
        # check boundary conditions
        if child_taxid == parent_taxid:
            return False
        if parent_taxid == self.root:
            return True

        tmp_parent_taxid = child_taxid
        while True:
            if not self.parent_nodes.has_key(tmp_parent_taxid):
                return False
            tmp_parent_taxid = self.parent_nodes[tmp_parent_taxid]
            if tmp_parent_taxid == self.root:
                return False
            if tmp_parent_taxid == parent_taxid:
                return True

    def find_lca (self, taxid_list):
        ''' Finds the lowest common ancestor of
            a list of nodes
        '''
        if len(taxid_list) == 0:
            raise Exception ("taxid_list is empty, cannot find LCA!")

        # each of the visited nodes remembers how many
        # child nodes traversed it
        self.num_visited        = defaultdict(int)

        current_taxids     = taxid_list
        num_of_nodes       = len(current_taxids)

        # first check if all nodes exist (and sum up blast scores)
        for i in range (0, len(taxid_list)):
            taxid       = taxid_list[i]
            if not self.parent_nodes.has_key(taxid):
                raise Exception ("Key error, no element with id %d." % taxid)

        # now find the lowest common ancestor
        while (True):

            parent_nodes = []
            for taxid in current_taxids:
                # root node must not add itself to parent list
                if   taxid != self.root:    parent_taxid = self.parent_nodes[taxid]
                else:                        parent_taxid = None
                # if parent exists, append him to parent list
                # duplicates ensure that every traversal will count
                if parent_taxid:            parent_nodes.append(parent_taxid)

                self.num_visited[taxid]     += 1
                if self.num_visited[taxid] == num_of_nodes:

                    self.lca_root = taxid
                    return taxid
            # refresh current nodes
            current_taxids = parent_nodes


    def get_relevant_taxid (self, tax_id):
        return self.tax2relevantTax.get(tax_id, -1)

    def get_lineage(self,tax_id):
        lineage = []
        while (True):
            if tax_id == self.root:
                break
            lineage.append(tax_id)
            tax_id = self.parent_nodes[tax_id]
        return reversed(lineage)

    def get_parent_with_rank(self, tax_id, rank):
        parent = 0
        while (True):
            if tax_id == self.root:
                return 0
            if self.nodes[tax_id].rank == rank:
                return tax_id
            tax_id = self.parent_nodes[tax_id]

    def _h_get_tax_nodes        (self, parent2child_fname):
        '''Loads the taxonomy nodes in a dictionary
           mapping the child to parent node.
        '''
        # file format per line: child_taxid parent_taxid
        with open(parent2child_fname) as fd:
            d = dict(self._h_from_parent_child_str (line) for line in fd)
        return d

    def _h_from_parent_child_str (self, line):
        '''Loads two integers (taxids) from a line
        '''
        key, sep, value = line.strip().partition(" ")
        if key == value: self.root = int(key)
        return int(key), int(value)


    def _h_find_taxnode_file(self, file_type):
        ''' Searches recursively through the current
            working directory to find the ncbi_tax_tree file.
        '''
        if file_type == 'parent2child':
            for root, dirs, files in os.walk (os.getcwd()):
                if 'ncbi_tax_tree' in files:
                    return root + ''.join(dirs) + '/ncbi_tax_tree'
        elif file_type == 'taxdata':
            for root, dirs, files in os.walk (os.getcwd()):
                if 'taxid2namerank' in files:
                    return root + ''.join(dirs) + '/taxid2namerank'
        else:
            raise ValueError('Internal error: looking for %s file type failed. No such file type supported.' % file_type)



    def _h_populate_child_nodes (self):
        ''' Populates child nodes from parent to child
            mapping dictionary
        '''
        child_nodes = defaultdict(list)
        for (child, parent) in self.parent_nodes.items():
            child_nodes[parent].append(child)
        return child_nodes

    def _h_set_relevant_taxonomy_assignments (self):
        ''' Sets some of the more important taxonomy
            assignments which can help in checking which kingdom
            an organism belongs to.
        '''
        import ncbi.taxonomy.organisms as orgs
        for organism_name in dir(orgs):
            if organism_name.startswith('__'):
                continue
            setattr(self, organism_name, getattr(orgs, organism_name))
        self.potential_hosts = [self.human,
                                self.mouse,
                                self.rats,
                                self.rodents,
                                self.primates,
                                self.animalia,
                                self.green_plants]

        self.microbes =        [self.archea,
                                self.bacteria,
                                self.viruses,
                                self.fungi,
                                self.euglenozoa,
                                self.alveolata,
                                self.amoebozoa,
                                self.fornicata,
                                self.parabasalia,
                                self.heterolobosea,
                                self.viroids,
                                self.stramenopiles,
                                self.cryptomycota,
                                self.entomophthoromycota,
                                self.microsporidia,
                                self.neocallimastigomycota]


    def _h_map_taxids_to_relevant_tax_nodes(self):
        host_nodes = list(self.potential_hosts)
        microbe_nodes = list(self.microbes)

        self.tax2relevantTax = {}
        for microbe_node in self.microbes:
            microbe_children = self._h_list_all_children(microbe_node)
            for child in microbe_children:
                self.tax2relevantTax[child] = microbe_node

        for host_node in self.potential_hosts:
            host_children = self._h_list_all_children(host_node)
            for child in host_children:
                self.tax2relevantTax[child] = host_node

        tagged_nodes    = self.tax2relevantTax.keys()
        all_nodes       = self.parent_nodes.keys()
        untagged_nodes  = set(all_nodes).difference(tagged_nodes)
        for node in untagged_nodes:
            self.tax2relevantTax[node] = -1

    def _h_list_all_children(self, tax_id):
        if not self.child_nodes.has_key(tax_id):
            return []
        one_step_children = self.child_nodes[tax_id]
        all_children = []
        while (True):
            if not one_step_children:
                break
            new_one_step_children = []
            all_children.extend(one_step_children)
            for child in one_step_children:
                if self.child_nodes.has_key(child):
                    new_one_step_children.extend(self.child_nodes[child])
            one_step_children = new_one_step_children
        return all_children



class TaxNode (object):
    '''
    Taxonomy nodes hold information on relevant
    taxonomy data, which are:
    * organism name (scientific)
    * taxonomy rank
    * score (arbitrary)
    '''


    def __init__(self, organism_name, rank, score=0.):
        self.organism_name = organism_name
        self.rank = rank
        self.score = score




