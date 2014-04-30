import xml.etree.ElementTree as ET
import sys,os
sys.path.append(os.getcwd())

from utils import enum
from utils.autoassign import autoassign

MISSING_TAXID = -1

class Organism (object):
    '''
    Placeholder for organism data from XML as described here:
    http://commondatastorage.googleapis.com/solverfiles/Organisms_XML_Schema_overview.pdf
    '''

    organismType = enum(GENUS='genus',
                        GENUS_SPECIES='genus_species',
                        GENUS_SPECIES_STRAIN='genus_species_strain',
                        ORGANISM_NAME='organism_name',
                        NEAREST_NEIGHBOR='nearest_neighbor',
                        MISSING_TAXID_ORG='missing_taxid_org'
                        )
    @autoassign
    def __init__(self, relativeAmount, count, taxon_id, taxonomy, type, organismName=None, genus=None,
                species=None, strain=None, nearestNeighbor=None, reads=None, genes=None):
        '''
        :param relativeAmount: (float)
        :param count: (int)
        :param taxon_id: (int)
        :param taxonomy: list of taxonomical node names (lineage)
        :param type:  (str) level of taxonomic assignment (check organismType)
        :param reads: list of read IDs or None for Host
        :param genes: list of gene objects mapped to this organism or None
                      if there are no reported genes
        '''

    def __hash__(self):
        return hash(self.taxon_id)

    def __eq__(self, other):
        '''Equality - by taxon id
        '''
        return self.taxon_id == other.taxon_id

    @classmethod
    def from_xml_organism_node(org, organism_node):
        '''
        Loads available organism data from organism XML node
        :param organism_node: xml.etree.ElementTree.Element instance
        '''
        # determine organism count
        relative_amount_node = organism_node.find('relativeAmount')
        relative_amount = eval(relative_amount_node.text)
        count = eval(relative_amount_node.attrib['count'])

        genes = Organism.get_genes(organism_node)
        reads = Organism.get_reads(organism_node)

        # determine taxonomy
        (taxon_id, taxonomy) = Organism.determine_taxonomy(organism_node)
        # determine organism type and available naming
        if taxon_id == MISSING_TAXID:
            organism = Organism(relative_amount, count, taxon_id, taxonomy, org.organismType.MISSING_TAXID_ORG, genes=genes, reads=reads)
            return organism

        (organism_type, names) = Organism.determine_org_type (organism_node)
        if taxon_id == MISSING_TAXID:
            organims_type = Organism.organismType.MISSING_TAXID_ORG

        nearest_neighbor = names['nearestNeighbor']
        org_name         = names['organismName']
        strain           = names['strain']
        species          = names['species']
        genus            = names['genus']
       
        organism = Organism(relative_amount, count, taxon_id, taxonomy, organism_type, org_name, genus,
                species, strain, nearest_neighbor, reads, genes)
        return organism


    @classmethod
    def determine_taxonomy(org, organism_node):
        '''
        Determines tax_id and taxonomical lineage from organism XML node.

        :param organism_node: xml.etree.ElementTree.Element intance
        :rtype: tuple(tax_id(int), taxonomy(list))
        '''
        taxonomy_node = organism_node.find('taxonomy')
        try:
            taxon_id = eval(taxonomy_node.attrib['taxon_id'])
        except KeyError:
            # so it is a nearest neighbor node
            log.info('No taxon ID info for organism with taxonomy %s' % taxonomy_node.text)
        taxonomy_str = taxonomy_node.text
        if taxonomy_str is None:
            return (MISSING_TAXID, [])
        if taxonomy_str.endswith('.'):
            taxonomy_str = taxonomy_str[0:-1]
        taxonomy = taxonomy_str.split('; ')
        return (taxon_id, taxonomy)

    @classmethod
    def get_genes(org, organism_node):
        '''
        Creates a list of genes (:class:`Gene`) from organism XML node.

        :param organism_node: :class:`xml.etree.ElementTree.Element` instance
        :rtype: list of :class:`Gene` objects or empty list if no
                genes were reported
        '''
        genes = []
        genes_node = organism_node.find('genes')
        if genes_node is None:
            return genes
        for gene_node in genes_node:
            genes.append(Gene.from_xml_gene_node(gene_node))
        return genes

    @classmethod
    def get_reads(org, organism_node):
        '''
        Creates a list of reads (str) from organism XML node.

        :param organism_node: xml.etree.ElementTree.Element intance
        :rtype: list of strings (read IDs)
        '''
        reads = []
        reads_node = organism_node.find('reads')
        if reads_node is None:
            return reads
        for read_node in reads_node:
            reads.append(read_node.text)
        return reads

    @classmethod
    def determine_org_type(org, organism_node):
        '''
        Determines taxonomical assignment type. As described here:
        http://commondatastorage.googleapis.com/solverfiles/Organisms_XML_Schema_overview.pdf
        there are four types of tax assignment:
        * nearest neighbor
        * only organism name (viruses, plasmids)
        * genus
        * genus, species
        * genus, species, strain

        :param organism_node: xml.etree.ElementTree.Element intance
        :rtype: tuple(organism_type, names). Organism type is value from
                organismType enum, and names is dictionary with
                string keys: 'nearestNeighbor', 'organismName', 'strain',
                'species', 'genus'
        '''
        neighbor_node = organism_node.find('nearestNeighbor')
        org_name_node = organism_node.find('organismName')
        strain_node   = organism_node.find('strain')
        species_node  = organism_node.find('species')
        genus_node    = organism_node.find('genus')
        nearest_neighbor = None
        organism_name    = None
        strain           = None
        species          = None
        genus            = None

        # ruzne sifrice
        if neighbor_node is not None:
            organism_type = org.organismType.NEAREST_NEIGHBOR
            nearest_neighbor = neighbor_node.text
        else:
            organism_name = org_name_node.text
            if strain_node is not None:
                organism_type = org.organismType.GENUS_SPECIES_STRAIN
                strain = strain_node.text
            else:
                if species_node is not None:
                    species = species_node.text
                    organism_type = org.organismType.GENUS_SPECIES
                else:
                    if genus_node is not None:
                        genus = genus_node.text
                        organism_type = org.organismType.GENUS
                    else:
                        organism_type = org.organismType.ORGANISM_NAME
        names = {'nearestNeighbor':nearest_neighbor, 'organismName':organism_name,
                'strain': strain, 'species':species, 'genus':genus}
        return (organism_type, names)






class Gene (object):
    '''
    Placeholder for gene data as described here:
    http://commondatastorage.googleapis.com/solverfiles/Organisms_XML_Schema_overview.pdf
    '''

    attributes = ['protein_id',
                  'locus_tag',
                  'product',
                  'ref_name',
                  'ref_start',
                  'ref_end']

    @classmethod
    def from_xml_gene_node (gene, gene_node):
        '''
        Generates Gene object from gene XML node.

        :param gene_node: :class:`xml.etree.ElementTree.Element` instance
        :rtype :class:`Gene`
        '''
        gene = Gene(gene_node.text)
        for key, value in gene_node.attrib.items():
            setattr(gene, key, value)
        for attribute in set(Gene.attributes).difference(set(gene_node.attrib.keys())):
            setattr(gene, attribute, None)
        return gene

    def __init__ (self, name):
        self.name = name

def loadOrganismData (xml_file):
    '''
    Loads all organisms from XML solution file.

    :param xml_file: path to XML solution file. File has to be of format as described
    here: http://commondatastorage.googleapis.com/solverfiles/Organisms_XML_Schema_overview.pdf
    :rtype: list of :class:`Organism` objects.
    '''
    tree = ET.parse(xml_file)
    root_node = tree.getroot()

    organisms = []
    organisms_node = root_node.find('organisms')
    for organism_node in organisms_node:
        organisms.append(Organism.from_xml_organism_node(organism_node))
    return organisms
