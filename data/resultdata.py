from utils import enum
from data.read import Read
from formats.xml_output import Gene, Organism, Read

class BinnedRead (object):
    def __init__(self, read_id, target_tax_id=None, binning_status=None, mapping_status=None):
        self.read_id        = read_id
        self.target_tax_id  = target_tax_id
        self.binning_status = binning_status
        self.mapping_status = mapping_status
    def set_target_organism(self, target_tax_id):
        self.target_tax_id = target_tax_id
    def set_binning_status(self, binning_status):
        self.binning_status = binning_status
    def set_mapping_status(self, mapping_status):
        self.mapping_status = mapping_status
    def to_xml_read(self):
        xml_read = Read(self.read_id)
        return xml_read

class IdentifiedCds (object):
    def __init__(self, cds, binned_reads=[]):
        self.cds = cds
        self.binned_reads = binned_reads
    def set_binned_reads(self, binned_reads):
        self.binned_reads = binned_reads
    def add_binned_read(self, binned_read):
        self.binned_reads.append(binned_read)
    def to_xml_gene(self):

        xml_gene = Gene(self.cds.protein_id, 
                        self.cds.locus_tag, 
                        self.cds.product, 
                        self.cds.protein_id, 
                        self.cds.name)
        return xml_gene


class Organism (object):
    def __init__(self, tax_id, name, rank):
        self.tax_id = tax_id
        self.name = name
        self.rank = rank
        self.identified_coding_regions = {}
        self.reads_aligned_to_noncoding_regions = []
        self.reads_aligned_to_coding_regions = []
        self.ambiguous_organism_binning_reads = []
        self.ambiguous_coding_region_mapping_reads = []

    def add_identified_coding_region(self, identified_cds):
        self.identified_coding_regions[identified_cds.cds] = identified_cds
        self.reads_aligned_to_coding_regions.extend(identified_cds.binned_reads)

    def add_read_aligned_to_noncoding_region(self, binned_read):
        self.reads_aligned_to_noncoding_regions.append(binned_read)

    def set_reads_aligned_to_noncoding_region(self, binned_reads):
        self.reads_aligned_to_noncoding_regions = binned_reads

    def add_ambiguous_organism_read(self, binned_read):
        self.ambiguous_organism_binning_reads.append(binned_read)

    def set_ambiguous_organism_reads(self, binned_reads):
        self.ambiguous_organism_binning_reads = binned_reads

    def add_ambiguous_coding_region_mapped_read(self, binned_read):
        self.ambiguous_coding_region_mapping_reads.append(binned_read)

    def set_ambiguous_coding_region_mapped_reads(self, binned_reads):
        self.ambiguous_coding_region_mapping_reads = binned_reads

    def contains_identified_coding_region(self, cds):
        if self.identified_coding_regions.has_key(cds):
            return True
        else:
            return False

    def get_reads(self):
            return self.reads_aligned_to_noncoding_regions + self.reads_aligned_to_coding_regions


    def to_xml_organism(self, tax_tree):
        xml_genes = []
        for identified_cds in self.identified_coding_regions.values():
            xml_genes.append(identified_cds.to_xml_gene())
        xml_reads = []
        for binned_read in self.reads_aligned_to_noncoding_regions+self.reads_aligned_to_coding_regions:
            xml_reads.apppend(binned_read.to_xml_read())
            #amount_count, amount_relative, taxon_id, taxonomy, name, genus, species, genes, variants, reads, is_host=False
        amount_count = len(self.reads_aligned_to_coding_regions + self.reads_aligned_to_noncoding_regions)
        lineage = tax_tree.get_lineage(self.tax_id)
        lineage_names = []
        for tax_id in lineage:
            lineage_names.append(tax_tree.nodes[tax_id])
        taxonomy = '; '.append(lineage_names)
        genus_taxid = tax_tree.get_parent_with_rank(self.tax_id, 'genus')
        species_taxid = tax_tree.get_parent_with_rank(self.tax_id, 'species')
        genus = tax_tree.nodes[genus_taxid] if genus_taxid else ''
        species = tax_tree.nodes[species_taxid] if species_taxid else ''

        xml_organism = Organism(amount_count, amount_count, self.tax_id,
                                taxonomy, self.name, genus, species,
                                xml_genes, None, xml_reads)
        return xml_organism


