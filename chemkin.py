import numbers
import xml.etree.ElementTree as ET
import numpy as np
import doctest
import reaction_coeffs


class ElementaryRxn():
    '''
    member attributes:
        r_stoich: numeric list of lists or array
                  column length must equal length of x
                  number of columns indicates number of reactions
                  stoichiometric coefficients of reactants
        p_stoich: numeric list of lists or array
                  must be equal in shape to stoich_r
                  stoichiometric coefficients of products
        k:        a numeric value (or values, for multiple reactions)
                  reaction rate coefficient
    '''

    def __init__(self, filename):
        self.parse(filename)

    def parse(self, filename):
        root = ET.parse(filename).getroot()
        specieslist = root.find('phase').find('speciesArray').text.strip().split(' ')
        # commenting out the below because we won't be receiving the
        # concentrations from the input .xml file, however this could be useful
        # later if we might want to specify initial concentrations in the .xml
        #conc_list = list(map(lambda x: float(x), 
        #    root.find('phase').find('concentrations').text.strip().split(' ')))
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
                        k_list.append(reaction_coeffs.arrh(A, E))
                    else:
                        k_list.append(reaction_coeffs.arrh(A, E, float(T.text)))
                else:
                    if T is None:
                        k_list.append(reaction_coeffs.mod_arrh(A, float(b.text), E))
                    else:
                        k_list.append(reaction_coeffs.mod_arrh(A, float(b.text), E,\
                            float(T.text)))

        # the transpose here is just to have the stoichiometric coefficients for each
        # specie lie along a column, with each column being another reaction
        self.r_stoich = np.array(r_stoich).transpose()
        self.p_stoich = np.array(p_stoich).transpose()
        self.k = k_list
        # rxn = chemkin.ElementaryRxn(conc_list, r_stoich, k_list)
        # print(chemkin.rxn_rate(conc_list, r_stoich, k_list, p_stoich))
        # return np.array(rxn.rxn_rate(p_stoich))

    def prog_rate(self, x):
        '''
        Returns the progress rate for N species going through M
        irreversible, elementary reactions.

        INPUTS
        ======
        x:       numeric list
                 concentrations of A, B, C

        RETURNS
        =======
        omega:   the progress rate for the reaction, numeric

        EXAMPLES
        =======
        >>> ElementaryRxn('test/doctest1.xml').prog_rate([1, 2, 4])
        [20.0]

        >>> ElementaryRxn('test/doctest2.xml').prog_rate([1,2,1])
        [40.0, 10.0]
        '''
        k = self.k
        x = np.array(x)
        stoich = np.array(self.r_stoich)
        k = np.array(k)

        assert (np.issubdtype(x.dtype, np.number)), "Species concentrations must be numeric"
        assert (np.issubdtype(stoich.dtype, np.number)), "Stoichiometric coefficients must be numeric"
        assert (len(x) == stoich.shape[0]), "All species must have stoichiometric coefficients"

        assert (np.issubdtype(k.dtype, np.number)), "Reaction rate coefficients must be numeric"

#       if len(np.shape(stoich)) == 1: # Handle the one-reaction case
#           return list(k * np.product(x**stoich.T))
#       else:
        return list(k * np.product(x**stoich.T, axis=1))

    def rxn_rate(self, x):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        elementary, irreversible reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric

        EXAMPLES
        =======
        >>> ElementaryRxn('test/doctest3.xml').rxn_rate([1,2,1])
        array([-30., -60.,  20.])

        '''
        p_stoich = np.array(self.p_stoich)
        r_stoich = np.array(self.r_stoich)
#        assert(stoich_r.shape == stoich_p.shape), "Reactant and product stoich. coefficients must be same dimensions"

        # Get a list of progress rates for each reaction
        omega = self.prog_rate(x)

        # Multiply by difference in stoichiometric coefficients,
        # then sum across reactions.
        if np.shape(p_stoich)[1] == 1: # Handle the one-reaction case
            return np.sum(omega * (p_stoich - r_stoich))
        else:
            return np.sum(omega * (p_stoich - r_stoich), axis=1)

    def __str__(self):
        return("Stoichiometric coefficients of reactants: {}\n\
            Stoichiometric coefficients of reactants: {}\n\
            Reaction rate coefficient: {}".format(self.r_stoich,\
            self.p_stoich, self.k))

class ReversibleRxn(ElementaryRxn):

    def reversible_rxn_rate(x):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        reversible reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric
        '''
        raise NotImplementedError



class NonelRxn(ElementaryRxn):
    def nonel_rxn_rate(x):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        nonelementary reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric
        '''
        raise NotImplementedError

