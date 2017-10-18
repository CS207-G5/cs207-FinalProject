# python rxn_from_file.py input.xml
#
# given an input .xml file called input.xml with structure specified by the
# example file rxns.sml, this script computes the reaction rates for the
# information given in input assuming all the reactions are elementary and
# irreversible using the compute method from rxn_from_file

import rxn_from_file
import sys

# the name of the .xml file will be given in the first command line argument
if len(sys.argv) >= 2:
    print(rxn_from_file.compute(sys.argv[1]))

