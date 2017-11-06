# given an input .xml file, collect the information and write it into a rxnset object

import xml.etree.ElementTree as ET
import numpy as np
import reaction_coeffs

class rxnset:
    def _init__(self):
        self.r_stoich=None
        self.p_stoich=None
        self.k_list=None
        self.specieslist=None
        self.conc_list=None
        self.nuij=None

    def compute(self, filename, input_T=100):
        root = ET.parse(filename).getroot()
        specieslist = root.find('phase').find('speciesArray').text.strip().split(' ')
        conc_list = list(map(lambda x: float(x), 
            root.find('phase').find('concentrations').text.strip().split(' ')))
        k_list = []
        r_stoich = []
        p_stoich = []
        for reaction in root.find('reactionData').findall('reaction'):
            # the stoichiometric coefficients will be stored in these arrays according
            # to the ordering given in the variable specieslist
            r_coeffs = [0] * len(specieslist)
            p_coeffs = [0] * len(specieslist)
            # these 2 lists are the information of specie and stoichiometric coefficient
            # gleaned from the .xml file
            r_list = reaction.find('reactants').text.strip().split(' ')
            p_list = reaction.find('products').text.strip().split(' ')
            for r in r_list:
                specie_coeff = r.split(':')
                # this indexing ensures that the same ordering is kept as in specieslist
                r_coeffs[specieslist.index(specie_coeff[0])] = float(specie_coeff[1])
            for p in p_list:
                specie_coeff = p.split(':')
                p_coeffs[specieslist.index(specie_coeff[0])] = float(specie_coeff[1])
            r_stoich.append(r_coeffs)
            p_stoich.append(p_coeffs)
            ratecoeff = reaction.find('rateCoeff')
            if ratecoeff.find('Arrhenius') is None:
                # constant rate coeff
                k_list.append(reaction_coeffs.const(float(ratecoeff.find('k').text)))
            else:
                # Arrhenius
                if ratecoeff.find('Arrhenius').find('A') is None:
                    raise ValueError("A is not found")
                if ratecoeff.find('Arrhenius').find('E') is None:
                    raise ValueError("E is not found")
                A = float(ratecoeff.find('Arrhenius').find('A').text)
                E = float(ratecoeff.find('Arrhenius').find('E').text)
                T = ratecoeff.find('Arrhenius').find('T')
                b = ratecoeff.find('Arrhenius').find('b')

                if b is None:
                    if T is None:
                        k_list.append(reaction_coeffs.arrh(A, E, input_T))
                    else:
                        k_list.append(reaction_coeffs.arrh(A, E, float(T.text)))
                else:
                    if T is None:
                        k_list.append(reaction_coeffs.mod_arrh(A, float(b.text), E,input_T))
                    else:
                        k_list.append(reaction_coeffs.mod_arrh(A, float(b.text), E,\
                            float(T.text)))

        # the transpose here is just to have the stoichiometric coefficients for each
        # specie lie along a column, with each column being another reaction
        r_stoich = np.array(r_stoich).transpose()
        p_stoich = np.array(p_stoich).transpose()
        self.r_stoich=r_stoich
        self.p_stoich=p_stoich
        self.k_list=k_list
        self.specieslist=specieslist
        self.conc_list=conc_list
        self.nuij=self.p_stoich-self.r_stoich

if __name__ == "__main__":
    x=rxnset()
    x.compute("rxns.xml",input_T=100)
    print("species: {}\nconcentration:{}\nr:{}\np:{}\nnuij:{}\nk:{}\n ".format(x.specieslist,x.conc_list,x.r_stoich,x.p_stoich,x.nuij, x.k_list))

