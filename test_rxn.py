#import rxn_from_file as rx

import numpy as np
import xml.etree.ElementTree as ET
import sys
import reaction_coeffs
import chemkin

def test_normal():
    assert list(chemkin.ElementaryRxn('rxns.xml').rxn_rate([1, 2, 3, 4, 5]))==[2.058282074119749e-09, -2.058282074119749e-09, 2.058282074119749e-09, -2.058282074119749e-09, -4.6705014964287609e-28]

def test_A():
    try:
        chemkin.ElementaryRxn('rxns_a.xml')
    except ValueError as err:
        assert(True)

def test_E():
    try:
        chemkin.ElementaryRxn('rxns_e.xml')
    except ValueError as err:
        assert(True)

def test_B():
    assert list(chemkin.ElementaryRxn('rxns_b.xml').rxn_rate([1, 2, 3, 4, 5]))==[8.194168527276232e-15, -8.194168527276232e-15, 8.1941685272996958e-15, -8.1941685272879639e-15, -1.1731769337224591e-26]

def test_T():
    assert list(chemkin.ElementaryRxn("rxns_t.xml").rxn_rate([1, 2, 3, 4, 5]))==[53.811949752854751, -53.811949752854751, 53.814583099716863, -53.813266426285807, -0.0013166734310555181]
