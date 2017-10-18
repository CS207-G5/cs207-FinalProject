
import chemkin
import numpy as np
import cmath

client_x = [2, 1, 0.5, 1, 1]
client_temps = [750, 1500, 2500]
k1 = [chemkin.reaction_coeffs.mod_arrh(A=10**8, b=0.5, E=5*10**4, T=i) for i in client_temps]
k2 = [chemkin.reaction_coeffs.const(10**4)]*3
k3 = [chemkin.reaction_coeffs.arrh(A=10**7, E=10**4, T=i) for i in client_temps]
client_sr = [ [2,0,0], [1,0,1], [0,1,0], [0,1,0], [0,0,1] ]
client_sp = [ [1,0,0], [0,1,0], [2,0,1], [0,0,1], [0,1,0] ]
    
def test_rxnrate():
    er = chemkin.ElementaryRxn(client_x, client_sr, 10)
    assert (all(er.rxn_rate(stoich_p = client_sp) == [-40., -45.,  85.,   5., -5.]))
    
def test_no_stoichproducts():
    assert (all(chemkin.rxn_rate(x = client_x, k = 10,
                             stoich_r = client_sr) == [-80., -50.,  -5.,  -5., -10]))
    
def test_onerxn():
    assert (chemkin.rxn_rate(x = [1,2,4], stoich_r = [2,1,0], k =10) == -60)
    
def test_zerotemp_mod():
    assert(chemkin.reaction_coeffs.mod_arrh(A=10**8, b=0.5, E=5*10**4, T=0) == 0.0)

def test_zerotemp_reg():
    assert(chemkin.reaction_coeffs.arrh(A=10**5, E=5*10**2, T=0) == 0.0)

def test_complex():
    try: 
        chemkin.reaction_coeffs.mod_arrh(A=10**8, b=cmath.sqrt(-5), E=5*10**4, T=0) == 0.0
    except ValueError as err:
        assert(type(err) == ValueError)
