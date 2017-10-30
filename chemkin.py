import reaction_coeffs
import numbers
import numpy as np
import doctest


class ElementaryRxn():

    def __init__(self, stoichr, k, stoichp=None):
 
        self.stoichr = stoichr
        self.stoichp=stoichp 
        self.k = k
        self.x=0
    def prog_rate(self,x):
        '''
        Returns the progress rate for N species going through M
        irreversible, elementary reactions.

        INPUTS
        ======
        x:       numeric list
                 concentrations of A, B, C
        stoichr:  numeric list of lists or array, equal column length to x
                 stoichiometric coefficients of reactants A, B, C
        k:       a numeric value
                 reaction rate coefficient

        RETURNS
        =======
        omega:   the progress rate for the reaction, numeric

        EXAMPLES
        =======
        >>> ElementaryRxn([1, 2, 4], 10).prog_rate([2, 1, 0])
        20

        >>> ElementaryRxn(stoichr = [[1,2],[2,0],[0,2]], k = 10).prog_rate([1,2,1])
        [40, 10]
        '''

        stoich = self.stoichr
        k = self.k
        self.x=x
        x = np.array(x)
        stoich = np.array(stoich)

        assert (np.issubdtype(x.dtype, np.number)), "Species concentrations must be numeric"
        assert (np.issubdtype(stoich.dtype, np.number)), "Stoichiometric coefficients must be numeric"
        assert (len(x) == stoich.shape[0]), "All species must have stoichiometric coefficients"

        if isinstance(k, int) and len(stoich.shape) > 1:

            k = [k]*stoich.shape[1]
        k = np.array(k)
        assert (np.issubdtype(k.dtype, np.number)), "Reaction rate coefficients must be numeric"

        if len(np.shape(stoich)) == 1:
            return k * np.product(x**stoich)
        else:
            omega = k * np.product(x**stoich.T, axis=1)
            return omega.tolist()

    def rxn_rate(self,  x):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        elementary, irreversible reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants
        stoich_r: numeric list of lists or array
                  column length must equal length of x
                  number of columns indicates number of reactions
                  stoichiometric coefficients of reactants
        stoich_p: numeric list of lists or array
                  must be equal in shape to stoich_r
                  stoichiometric coefficients of products
        k:        a numeric value (or values, for multiple reactions)
                  reaction rate coefficient

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric

        EXAMPLES
        =======
        >>> ElementaryRxn( stoichr = [[1,0],[2,0],[0,2]], k=10,stoichp = [[0,1],[0,2],[1,0]]).rxn_rate(x = [1,2,1])
        array([-30, -60,  20])

        >>> ElementaryRxn([2,1,0], 10).rxn_rate([1,20,4])
        -600.0
        '''
        stoich_r = self.stoichr
        stoich_p = self.stoichp
  
        k = self.k
        self.x=x
        if stoich_p is None:
            stoich_p = np.zeros(np.shape(stoich_r))
        else:
            stoich_p = np.array(stoich_p)

        stoich_r = np.array(stoich_r)
#        assert(stoich_r.shape == stoich_p.shape), "Reactant and product stoich. coefficients must be same dimensions"

        # Get a list of progress rates for each reaction
        omega = self.prog_rate(x)

        # Multiply by difference in stoichiometric coefficients,
        # then sum across reactions.
        if len(np.shape(stoich_p)) == 1: # Handle the one-reaction case
            return np.sum(omega * (stoich_p - stoich_r))
        else:
            return np.sum(omega * (stoich_p - stoich_r), axis=1)

    def __str__(self):
        return("Stoichiometric coefficients of reactants: {}\n\
            Concentrations: {}\n\
            Reaction rate coefficient(s): {}".format(self.stoich, self.x,\
                self.k))

class ReversibleRxn(ElementaryRxn):

    def reversible_rxn_rate(stoich_r, stoich_p, k ):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        reversible reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants
        stoich_r: numeric list of lists or array
                  column length must equal length of x
                  number of columns indicates number of reactions
                  stoichiometric coefficients of reactants
        stoich_p: numeric list of lists or array
                  must be equal in shape to stoich_r
                  stoichiometric coefficients of products
        k:        a numeric value (or values, for multiple reactions)
                  reaction rate coefficient

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric
        '''
        raise NotImplementedError



class NonelRxn(ElementaryRxn):
    def nonel_rxn_rate(stoich_r, stoich_p,k):
        '''
        Returns the reaction rate, f, for each specie (listed in x)
        through one or multiple (number of columns in stoich_r)
        nonelementary reactions.

        f = sum(omega_j*nu_ij) for i species in j reactions.

        INPUTS
        ======
        x:        numeric list or array
                  concentrations of reactants
        stoich_r: numeric list of lists or array
                  column length must equal length of x
                  number of columns indicates number of reactions
                  stoichiometric coefficients of reactants
        stoich_p: numeric list of lists or array
                  must be equal in shape to stoich_r
                  stoichiometric coefficients of products
        k:        a numeric value (or values, for multiple reactions)
                  reaction rate coefficient

        RETURNS
        =======
        f:        the reaction rate for each specie, numeric
        '''
        raise NotImplementedError

