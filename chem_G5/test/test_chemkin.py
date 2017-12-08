import sys
sys.path.append('..')  # Needs to be before chemkin import
from chem_G5 import chemkin
import cmath
import numpy as np

client_x = [2, 1, 0.5, 1, 1]
client_temps = [750, 1500, 2500]
# We don't seem to use these k values. Should we or are they no longer useful?
k1 = [chemkin.reaction_coeffs.mod_arrh(A=10**8, b=0.5, E=5*10**4, T=i) for i in client_temps]
k2 = [chemkin.reaction_coeffs.const(10**4)]*3
k3 = [chemkin.reaction_coeffs.arrh(A=10**7, E=10**4, T=i) for i in client_temps]
client_sr = [ [2,0,0], [1,0,1], [0,1,0], [0,1,0], [0,0,1] ]
client_sp = [ [1,0,0], [0,1,0], [2,0,1], [0,0,1], [0,1,0] ]

# Tests for basic ElementaryRxn()

def test_rxnrate():
    assert (all(chemkin.ElementaryRxn('chem_G5/test/rxnrate.xml').rxn_rate(client_x, 1200) == [-40., -45.,  85.,   5., -5.]))

def test_onerxn():
    assert (chemkin.ElementaryRxn('chem_G5/test/onerxn.xml').rxn_rate([1, 2, 4], 1200) == 60.0)

def test_zerotemp_mod():
    assert(chemkin.reaction_coeffs.mod_arrh(A=10**8, b=0.5, E=5*10**4, T=0) == 0.0)

def test_zerotemp_reg():
    assert(chemkin.reaction_coeffs.arrh(A=10**5, E=5*10**2, T=0) == 0.0)

def test_complex():
    try:
        chemkin.reaction_coeffs.mod_arrh(A=10**8, b=cmath.sqrt(-5), E=5*10**4, T=0) == 0.0
    except ValueError as err:
        assert(type(err) == ValueError)

def test_parse_nonelem():
    r = chemkin.ElementaryRxn('chem_G5/test/parse_nonelem.xml')
    assert (r.efficiencies == [[2.5, 12.0, 1.0, 1.0, 1.0], [2.0, 11.0, 1.0, 0.78, 1.0]])
    assert (r.rxnparams == [[4.577e+19, 104380.0, -1.4], [6.366e+20, 524.8, -1.72, 1475000000000.0, 0.0, 0.6, 0.8, 1e+30, 1e+30, 1e-30]])

# Tests for ReversibleRxn()

def test_creation():
    rever = chemkin.ReversibleRxn('chem_G5/test/rxnrate.xml')
    assert(rever.s == ['H', 'O', 'OH', 'H2', 'O2'])
    assert((rever.r == [[2,  0,  0],
                        [1,  0,  1],
                        [0,  1,  0],
                        [0,  1,  0],
                        [0,  0,  1]]).all())

    assert((rever.p == [[1,  0,  0],
                        [0,  1,  0],
                        [2,  0,  1],
                        [0,  0,  1],
                        [0,  1,  0]]).all())
    assert((rever.nuij == [[-1,  0,  0],
                           [-1,  1, -1],
                           [2, -1,  1],
                           [0, -1,  1],
                           [0,  1, -1]]).all())
    assert(rever.p0 == 1.0e+05)
    assert(rever.R == 8.3144598)
    assert((rever.gamma == [0, 0, 0]).all())

rever = chemkin.ReversibleRxn('chem_G5/test/rxnrate.xml')

def test_read_db():
    rever.read_SQL(100)
    assert((rever.nasa[0] == [2.50000001e+00, -2.30842973e-11,
                              1.61561948e-14, -4.73515235e-18,
                              4.98197357e-22,  2.54736599e+04,
                              -4.46682914e-01]).all())

def test_trange_lesser():
    rever.read_SQL(2500)
    assert((rever.nasa[0] == [2.50000000e+00, 7.05332819e-13,
                              -1.99591964e-15, 2.30081632e-18,
                              -9.27732332e-22, 2.54736599e+04,
                              -4.46682853e-01]).all())

def test_badT():
    try:
        rever.read_SQL('I am not a number!')
    except AssertionError as err:
        assert(type(err) == AssertionError)

def test_Cp():
    out = rever.Cp_over_R(100)
    # Note the use of np.isclose() to handle VERY tiny rounding diffs.
    assert(np.isclose(out, np.array([2.5, 2.90084906,
                                     3.79431834, 2.96702128,
                                     3.57189904])).all())

def test_H():
    out = rever.H_over_RT(2500)
    assert(np.isclose(out, np.array([12.68946396, 17.12513405, 7.5535748,
                    -7.51349287, 7.6502427])).all())

def test_S():
    out = rever.S_over_R(2000)
    assert(np.isclose(out, np.array([18.5555733, 24.97024176, 29.77900614,
                                     19.75486987,  33.26701985])).all())

def test_backwards_coeffs():
    rever.backward_coeffs(900)
    assert((rever.kb == [10, 10, 10]).all())

rev = chemkin.ReversibleRxn('chem_G5/test/rxns_rev.xml')
x=[10,10,10,10]

def test_prog_rate_reversible():
    np.testing.assert_allclose(rev.prog_rate(x,1200),[-1.22703266e+09])

def test_rxnrate_reversible():
    np.testing.assert_allclose(rev.rxn_rate(x,1200),[ 1.22703266e+09,   1.22703266e+09,  -1.22703266e+09,
        -1.22703266e+09])
