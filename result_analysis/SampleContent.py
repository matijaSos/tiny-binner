import xml.etree.ElementTree as ET

from result_analysis.solutiondata import Organism
from result_analysis.ValidationResult import ValidationResult

class SampleContent(object):
    '''
    Contains contents of a metagenomic sample.

    Essentialy, holds a list of organisms present in the sample and reads and genes
    assigned to each of them.

    This object can be either created as an output of a binning algorithm,
    or it can be loaded from an XML file. (e.g. when correct solution is available)

    This object follows DTRA style of binner output.
    '''

    def __init__(self, organisms):
        '''Constructor

        Args:
            organisms ([Organism]): List of organisms present in a sample
        '''
        self.organisms = organisms

    def validate_with(self, real):
        '''Compare estimated with correct content of a sample.

        Args:
            real (SampleContent): Real content
        Returns:
            (ValidationResult): Object holding the result of validation
        '''
        return ValidationResult(self, real)

    def get_orgs_num(self):
        '''Returns number of organisms in a sample.

        Args:
            None
        Returns:
            (int): number of organisms in a sample
        '''
        return len(self.organisms)

    # ------------------------ Factory methods ----------------------- #

    @classmethod
    def from_xml(cls, xml_file):
        '''Loads object from XML file - format defined by DTRA

        Args:
            xml_file (string): Path to the XML file holding sample contents.
        Returns:
            (SampleContent): New instance containing loaded data
        '''
        tree = ET.parse(xml_file)
        root_node = tree.getroot()

        organisms = []
        for organism_node in root_node.find('organisms'):
            organisms.append(Organism.from_xml_organism_node(organism_node))

        return cls(organisms)
            
    

