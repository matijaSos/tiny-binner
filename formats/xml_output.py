import sys
from utils.autoassign import autoassign
from xml.dom import minidom

class Dataset(object):

    @autoassign
    def __init__(self, desc_xml):
        ''' Dataset init
            @param name input file name
            @param host_genus host genus name
            @param host_species host species name
            @param common_name common name
            @param taxon_id host taxon id
            @param taxonomy host taxonomy
            @param sample_source e.g. Whole Blood
            @param sample_type e.g. DNA
            @param seq_method e.g. single-end
            @param sequncer Roche 454, Illumina, PacBio or IonTorrent
        '''
        pass

class Gene(object):

    @autoassign
    def __init__(self, protein_id, locus_tag, product, ref_name, name):
        ''' Gene init
            @param protein_id protein id e.g. AAS63914.1
            @param locus_tag e.g. YP_3766
            @param product e.g. putative carbohydrate kinase
            @param ref_name ref name
            @param name gene name e.g. xylB3
        '''
        pass
        
class Variant(object) :

    @autoassign
    def __init__(self, ref_name, ref_start, ref_seq, name, offset, context):
        ''' Variant init
            @param ref_name e.g. CAB54900.1
            @param ref_start e.g. 31
            @param ref_seq e.g. -
            @param name variant name e.g. A
            @param offset e.g. 28
            @param context e.g. ACTGGGGAGAGAGGAGCTTTTATATATATATTATAAGGCCC  
        '''
        pass

class Organism(object):

    @autoassign
    def __init__(self, amount_count, amount_relative, taxon_id, taxonomy, name,
                 genus, species, genes, variants, reads, is_host=False):
        ''' Organism init
            @param amount_count amount of reads e.g. 2156
            @param amount_relative relative to all others e.g. 0.158
            @param taxon_id e.g. 9606 for human
            @param taxonomy organism taxonomy
            @param name organism name e.g. Yersinia
            @param species organism species
            @param ([Gene]) genes list of genes
            @param ([Variant]) variants list of variants
            @param ([Read]) reads list of reads
            @param is_host a bool defining whether the organism is host or not
        '''
        pass

class Read(object):

    @autoassign
    def __init__(self, sequence):
        pass

class XMLOutput(object):

    def __init__(self, dataset, organisms, output_file_name=None):
        ''' XMLOutput init
            @param (Dataset) dataset structure with info about the dataset
            @param ([Organism]) organisms a list of all organims and coresponding data
            @param output_file_name path to xml output relative to binner/ e.g. 'formats/myout.xml'
        '''
        self.dataset = dataset
        self.organisms = organisms
        self.output_file_name = output_file_name

    def _get_xml_text(self, nodes):

        for node in nodes:
            if node.nodeType == node.TEXT_NODE:
                return node.data

    def _dataset_details_output(self, level):

        tab = " " * level * 2

        xmldoc = minidom.parse(self.dataset.desc_xml)
        
        nodes = xmldoc.getElementsByTagName('datasetName')[0].childNodes
        name = self._get_xml_text(nodes)

        nodes = xmldoc.getElementsByTagName('hostGenus')[0].childNodes
        host_genus = self._get_xml_text(nodes)

        nodes = xmldoc.getElementsByTagName('hostSpecies')[0].childNodes
        host_species = self._get_xml_text(nodes)

        nodes = xmldoc.getElementsByTagName('commonName')[0].childNodes
        common_name = self._get_xml_text(nodes)

        nodes = xmldoc.getElementsByTagName('taxonomy')[0].childNodes
        taxonomy = self._get_xml_text(nodes)

        taxon_id = xmldoc.getElementsByTagName('taxonomy')[0].attributes['taxon_id'].value

        nodes = xmldoc.getElementsByTagName('sampleSource')[0].childNodes
        sample_source = self._get_xml_text(nodes)
  
        nodes = xmldoc.getElementsByTagName('sampleType')[0].childNodes
        sample_type = self._get_xml_text(nodes)

        nodes = xmldoc.getElementsByTagName('sequencer')[0].childNodes
        sequencer = self._get_xml_text(nodes)

        seq_method = xmldoc.getElementsByTagName('sequencer')[0].attributes['method'].value
   
        if (name): 
            print(tab + "<datasetName>" + str(name) + "</datasetName>")
        if (host_genus):
            print(tab + "<hostGenus>" + str(host_genus) + "</hostGenus>")
        if (host_species):
            print(tab + "<hostSpecies>" + str(host_species) + "</hostSpecies>")
        if (common_name):
            print(tab + "<commonName>" + str(common_name) + "</commonName>")
        if (taxon_id):
            print(tab + "<taxonomy taxon_id=\"" + str(taxon_id)  + "\">" + str(taxonomy) + "</taxonomy>")
        if (sample_source):
            print(tab + "<sampleSource>" + str(sample_source) + "</sampleSource>")
        if (sample_type):
            print(tab + "<sampleType>" + str(sample_type) + "</sampleType>")
        if (seq_method and sequencer):
            print(tab + "<sequencer method=\"" + str(seq_method) + "\">" + str(sequencer) + "</sequencer>")

    def _dataset_output(self, level):
        
        tab = " " * level * 2

        print(tab + "<dataset>")
        self._dataset_details_output(level+1)
        print(tab + "</dataset>")

    def _gene_output(self, level, gene):

        tab = " " * level * 2

        gene_attributes = "protein_id =\"" + str(gene.protein_id) + "\""
        if(gene.locus_tag):
            gene_attributes = gene_attributes + " locus_tag=\"" + str(gene.locus_tag) + "\""
            
        if(not gene.product):
            gene.product = gene.protein_id
            
        gene_attributes = gene_attributes + " product=\"" + str(gene.product) + "\""

        if(not gene.ref_name):
            gene.ref_name = gene.protein_id

        gene_attributes = gene_attributes + " ref_name=\"" + str(gene.ref_name) + "\""

        if(not gene.name):
            gene.name = gene.protein_id

        print(tab + "<gene " + gene_attributes + ">" + str(gene.name) + "</gene>" )

    def _variant_details(self, level, variant):

        tab = " " * level * 2

        print(tab + "<variant ref_name=\"" + str(variant.ref_name) + "\" ref_start=\"" + str(variant.ref_start) + "\" ref_seq=\"" + str(variant.ref_seq) + "\">" + str(variant.name) + "</variant>")
        print(tab + "<context offset=\"" + str(variant.offset) + "\">" + str(variant.context) + "</context>")

    def _variant_output(self, level, variant):

        tab = " " * level * 2

        print(tab + "<sequenceDifference>")
        self._variant_details(level+1, variant)
        print(tab + "</sequenceDifference>")

    def _sequence_output(self, level, read):

        tab = " " * level * 2

        print(tab + "<sequence>" + str(read.sequence) + "</sequence>")

    def _organism_details_output(self, level, organism):

        tab = " " * level * 2

        print(tab + "<relativeAmount count=\"" + str(organism.amount_count) + "\">" + str("{0:.3f}".format(organism.amount_relative*100)) + "</relativeAmount>")
       
        xmldoc = minidom.parse(self.dataset.desc_xml)
        
        if organism.is_host:
            organism.taxon_id = xmldoc.getElementsByTagName('taxonomy')[0].attributes['taxon_id'].value
            organism.taxonomy = self._get_xml_text(xmldoc.getElementsByTagName('taxonomy')[0].childNodes)
        
        if (organism.taxon_id):
            print(tab + "<taxonomy taxon_id=\"" + str(organism.taxon_id) + "\">" + str(organism.taxonomy) + "</taxonomy>")
        if (organism.name):
            print(tab + "<organismName>" + str(organism.name) + "</organismName>")

        if organism.is_host:
            # rest of data not needed in this case
            return

        if (organism.genus):
            print(tab + "<genus>" + str(organism.genus) + "</genus>")
        if (organism.species):
            print(tab + "<species>" + str(organism.species) + "</species>")

        print(tab + "<genes>")
        for gene in organism.genes:
            self._gene_output(level+1, gene)
        print(tab + "</genes>")

        #print(tab + "<variants>")
        #for variant in organism.variants:
        #    self._variant_output(level+1, variant)
        #print(tab + "</variants>")

        print(tab + "<reads>")
        for read in organism.reads:
            self._sequence_output(level + 1, read)
        print(tab + "</reads>")

    def _organism_output(self, level):

        tab = " " * level * 2

        for organism in self.organisms:
            print(tab + "<organism>")
            self._organism_details_output(level+1, organism)
            print(tab + "</organism>")

    def _organisms_output(self, level):

        tab = " " * level * 2

        print(tab + "<organisms>")
        self._organism_output(level+1)
        print(tab + "</organisms>")

    def xml_output(self, level=0):
        ''' Generates the xml from data already present in dataset and organisms
            to stdout
            @param level default 0, start offset for level zero xml tags
        '''

        if (self.output_file_name):
            orig_stdout = sys.stdout
            out_file = file(self.output_file_name, 'w')
            sys.stdout = out_file
        print("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>")
        print("<organismsReport>")
        self._dataset_output(level+1);
        self._organisms_output(level+1);
        print("</organismsReport>")
        if (self.output_file_name):
            sys.stdout = orig_stdout
            out_file.close()
