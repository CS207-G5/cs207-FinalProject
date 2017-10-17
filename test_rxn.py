import rxn_from_file as rx

import numpy as np
import xml.etree.ElementTree as ET
import reaction_coeffs
import chemkin
import sys

def test_normal():
    assert all(rx.compute("rxns.xml")==[2.05828207e-09,-2.05828207e-09,2.05828207e-09,-2.05828207e-09,-4.67050150e-28])
def test_A():
    try:
        rx.compute("rxns_a.xml")    
    except ValueError as err:
        assert(True)
def test_E():
    try:
        rx.compute("rxns_e.xml")    
    except ValueError as err:
        assert(True)

def test_B():
    assert all(rx.compute("rxns_b.xml")==[8.19416853e-15,-8.19416853e-15,8.19416853e-15,-8.19416853e-15,-1.17317693e-26])
    
def test_T():
    assert all(rx.compute("rxns_t.xml")==[5.38119498e+01,-5.38119498e+01,5.38145831e+01,-5.38132664e+01,-1.31667343e-03])