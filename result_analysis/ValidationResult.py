
from result_analysis.solutiondata import Organism

class ValidationResult(object):
    '''
    Holds data acquired from validating estimated against real sample content.

    General format: correctly_estimated/all wrongly_est total_estimated
    '''

    def __init__(self, estimated, real):
        '''Constructor

        Args:
            estimated (SampleContent): Estimated content
            real (SampleContent): Real content
        '''

        # Keep references to contents
        self.estimated  = estimated
        self.real       = real
        
        # ------ Organisms ------ #
        self.orgs_total_est = self.estimated.get_orgs_num()
        self.orgs_TP, self.orgs_FP = self._validate_organisms()

        self.orgs_total_real = self.real.get_orgs_num()

        # ------ Reads ------ #
        # For each organism, have a list of good and bad reads (and missing?)
        
        # ------ Genes ----- #
        # The same thing for the genes

    def _validate_organisms(self):
        '''Do validation on an organism level.

        Args:
            None
        Returns:
            ([Organism], [Organism]): (orgs_TP, orgs_FP)

        orgs_TP: organisms that are detected and really are present
        orgs_FP: organisms that are detected but are not present
        '''
        # Get lists from SampleContent objects
        est_orgs    = self.estimated.organisms
        real_orgs   = self.real.organisms

        orgs_TP = list(set(est_orgs) & set(real_orgs))
        orgs_FP = list(set(est_orgs) - set(real_orgs))

        return (orgs_TP, orgs_FP)

    def short_output(self):
        print("------------ Organisms ---------------")
        print("total estimated: {0} TP/total_real: {1}/{2} FP: {3}".format( 
                self.orgs_total_est, 
                len(self.orgs_TP), self.orgs_total_real,
                len(self.orgs_FP))
             )
        
        
        


            
        
        


