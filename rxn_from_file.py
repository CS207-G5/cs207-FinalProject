# python rxn_from_file.py input.xml

import xml.etree.ElementTree as ET
import numpy as np
import reaction_coeffs
import chemkin
import sys

filename = sys.argv[1]
root = ET.parse(filename).getroot()
# throw error if we can't get this far
specieslist = root.find('phase').find('speciesArray').text.strip().split(' ')
conc_list = map(lambda x: float(x), 
    root.find('phase').find('concentrations').text.strip().split(' '))
k_list = []
r_stoich = []
p_stoich = []
for reaction in root.find('reactionData').findall('reaction'):
    r_coeffs = [0] * len(specieslist)
    p_coeffs = [0] * len(specieslist)
    r_list = reaction.find('reactants').text.strip().split(' ')
    p_list = reaction.find('products').text.strip().split(' ')
    for r in r_list:
        specie_coeff = r.split(':')
        r_coeffs[specieslist.index(specie_coeff[0])] = float(specie_coeff[1])
    for p in p_list:
        specie_coeff = p.split(':')
        p_coeffs[specieslist.index(specie_coeff[0])] = float(specie_coeff[1])
    r_stoich.append(r_coeffs)
    p_stoich.append(p_coeffs)
    ratecoeff = reaction.find('rateCoeff')
    if ratecoeff.find('Arrhenius') is none:
        # constant rate coeff
        k_list.append(reaction_coeffs.const(float(ratecoeff.find('k').text)))
    else:
        # Arrhenius
        A = float(ratecoeff.find('Arrhenius').find('A').text)
        E = float(ratecoeff.find('Arrhenius').find('E').text)
        T = ratecoeff.find('Arrhenius').find('T')
        b = ratecoeff.find('Arrhenius').find('b')
        if b is none:
            if T is none:
                k_list.append(reaction_coeffs.arrh(A, E))
            else:
                k_list.append(reaction_coeffs.arrh(A, E, float(T.text)))
        else:
            if T is none:
                k_list.append(reaction_coeffs.mod_arrh(A, float(b.text), E))
            else:
                k_list.append(reaciton_coeffs.mod_arrh(A, float(b.text), E,\
                    float(T.text)))

r_stoich = np.array(r_stoich).transpose()
p_stoich = np.array(p_stoich).transpose()
chemkin.rxn_rate(conc_list, r_stoich, k_list, p_stoich)
