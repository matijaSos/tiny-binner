import argparse
import sys,os
sys.path.append(os.getcwd())

from result_analysis.SampleContent import SampleContent

def parse_input_parameters():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compares two xml output files '
        'file')
    argparser.add_argument('correct_xml', help='Correct contents', 
                            type=str)
    argparser.add_argument('estimated_xml', help="Estimated solution", 
                            type=str)
    args = argparser.parse_args()
        
    error = False
    if not os.path.exists(os.path.expanduser(args.correct_xml)):
        print "Correct xml file %s doesn't exist" % args.correct_xml
        error = True
    if not os.path.exists(os.path.expanduser(args.estimated_xml)):
        print "Estimated xml file %s doesn't exist" % args.estimated_xml
        error = True
    if error:
        exit(-1)

    return args
                                            
def main():
    args = parse_input_parameters()

    correct   = SampleContent.from_xml(args.correct_xml)
    estimated = SampleContent.from_xml(args.estimated_xml)
    print("Loaded successfully!")

    vr = estimated.validate_with(correct)
    vr.short_output()

if __name__ == '__main__':
    main()
