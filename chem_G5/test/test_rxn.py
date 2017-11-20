import sys
sys.path.append('..')
import numpy as np
import xml.etree.ElementTree as ET
from chem_G5 import reaction_coeffs
from chem_G5 import chemkin

def test_normal():
    assert list(chemkin.ElementaryRxn('chem_G5/test/rxns.xml').rxn_rate([1, 2, 3, 4, 5], 1200))==[5042089.9481361862, -5042089.9481361862, 7005917.5594420284, -6024003.7537891073, -981913.80565292062]

def test_A():
    try:
        chemkin.ElementaryRxn('chem_G5/test/rxns_a.xml')
    except ValueError as err:
        assert(True)

def test_E():
    try:
        chemkin.ElementaryRxn('chem_G5/test/rxns_e.xml')
    except ValueError as err:
        assert(True)

def test_B():
    assert list(chemkin.ElementaryRxn('chem_G5/test/rxns_b.xml').rxn_rate([1, 2, 3, 4, 5], 1200))==[-140443045.44219205, 140443045.44219205, 140443045.50068799, -0.029247979488840573, -140443045.47144002]

def test_T():
    assert list(chemkin.ElementaryRxn("chem_G5/test/rxns_t.xml").rxn_rate([1, 2, 3, 4, 5], 300))==[53.811949752854751, -53.811949752854751, 53.814583099716863, -53.813266426285807, -0.0013166734310555181]
