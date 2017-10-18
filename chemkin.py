import numpy as np
import reacton_coeffs

def prog_rate(x, stoich, k):
    '''
    Returns the progress rate for N species going through M
    irreversible, elementary reactions.

    INPUTS
    ======
    x:       numeric list
             concentrations of A, B, C...
    stoich:  numeric list of lists or array, equal column length to x
             stoichiometric coefficients of reactants A, B, C...
    k:       a numeric value
             reaction rate coefficient

    RETURNS
    =======
    omega:   numpy array of progress rate(s) for the reaction, numeric

    EXAMPLES
    =======
    >>> prog_rate([1,2,4], [2,1,0], 10)
    20

    >>> prog_rate(stoich = [[1,2],[2,0],[0,2]], x = [1,2,1], k = 10)
    [40, 10]
    '''
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

def rxn_rate(x, stoich_r, k, stoich_p = None):
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
    f:        numpy array of reaction rate for each specie, numeric

    EXAMPLES
    =======
    >>> rxn_rate(x = [1,2,1], stoich_r = [[1,0],[2,0],[0,2]], stoich_p = [[0,1],[0,2],[1,0]], k=10)
    array([-30, -60,  20])

    >>> rxn_rate([1,20,4], [2,1,0], 10)
    -600.0
    '''
    if stoich_p is None:
        stoich_p = np.zeros(np.shape(stoich_r))
    else:
        stoich_p = np.array(stoich_p)

    stoich_r = np.array(stoich_r)
    assert(stoich_r.shape == stoich_p.shape), "Reactant and product stoich. coefficients must be same dimensions"

    # Get a list of progress rates for each reaction
    omega = prog_rate(x, stoich_r, k)

    # Multiply by difference in stoichiometric coefficients,
    # then sum across reactions.
    if len(np.shape(stoich_p)) == 1: # Handle the one-reaction case
        return np.sum(omega * (stoich_p - stoich_r))
    else:
        return np.sum(omega * (stoich_p - stoich_r), axis=1)


def reversible_rxn_rate(x, stoich_r, k, stoich_p = None):
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



def nonel_rxn_rate(x, stoich_r, k, stoich_p = None):
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
