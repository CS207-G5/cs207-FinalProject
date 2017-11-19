import sys
sys.path.append('..')  # Needs to be before chemkin import
import chemkin
import cmath

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
    assert (all(chemkin.ElementaryRxn('test/rxnrate.xml').rxn_rate(client_x) == [-40., -45.,  85.,   5., -5.]))

def test_onerxn():
    assert (chemkin.ElementaryRxn('test/onerxn.xml').rxn_rate([1, 2, 4]) == 60.0)

def test_zerotemp_mod():
    assert(chemkin.reaction_coeffs.mod_arrh(A=10**8, b=0.5, E=5*10**4, T=0) == 0.0)

def test_zerotemp_reg():
    assert(chemkin.reaction_coeffs.arrh(A=10**5, E=5*10**2, T=0) == 0.0)

def test_complex():
    try:
        chemkin.reaction_coeffs.mod_arrh(A=10**8, b=cmath.sqrt(-5), E=5*10**4, T=0) == 0.0
    except ValueError as err:
        assert(type(err) == ValueError)

# Tests for ReversibleRxn()

def test_creation():
    rever = chemkin.ReversibleRxn('test/rxnrate.xml')
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
    assert((rever.kf == [10, 10, 10]).all())
    assert(rever.p0 == 1.0e+05)
    assert(rever.R == 8.3144598)
    assert((rever.gamma == [0, 0, 0]).all())

def test_read_db():
    rever = chemkin.ReversibleRxn('test/rxnrate.xml')
    rever.read_SQL(100)
    assert((rever.nasa[0] == [[2.50000001e+00, -2.30842973e-11,
                               1.61561948e-14, -4.73515235e-18,
                               4.98197357e-22,  2.54736599e+04,
                               -4.46682914e-01]]).all())

def test_trange_greater():
    pass

def test_trange_lesser():
    pass

def test_getlowcoeffs():
    pass

def test_gethighcoeffs():
    pass

def test_badT():
    rever = chemkin.ReversibleRxn('test/rxnrate.xml')
    try:
        rever.read_SQL('I am not a number!')
    except AssertionError as err:
        assert(type(err) == AssertionError)

def test_Cp():
    pass

def test_H():
    pass

def test_S():
    pass

def test_backwards_coeffs():
    pass

def test_prog_rate_reversible():
    pass

def test_rxnrate_reversible():
    pass
