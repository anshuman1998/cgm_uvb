import astropy.table as tab
import numpy as np
import matplotlib.pyplot as plt

"""
                                                  Gas Phase Chemical Composition
        H :  0.0000  He: -1.0881  Li:-10.9508  Be:-10.6198  B : -9.3002  C : -3.5702  N : -4.1701  O : -3.3098  F : -7.4401
        Ne: -4.0701  Na: -5.7595  Mg: -4.4001  Al: -5.5498  Si: -4.4895  P : -6.5901  S : -4.8794  Cl: -6.5003  Ar: -5.6003
        K : -6.9706  Ca: -5.6596  Sc: -8.8508  Ti: -7.0501  V : -8.0701  Cr: -6.3595  Mn: -6.5702  Fe: -4.5003  Co: -7.0101
                                               Ni: -5.7799  Cu: -7.8097  Zn: -7.4401
"""

def get_gass10_abundances(element_name):
    abd = None
    if element_name == 'O':
        abd = 10 **(-3.3098)
    if element_name == 'C':
        abd = 10 **(-3.5702)
    if element_name == 'N':
        abd = 10 **(-4.1701)
    if element_name == 'Si':
        abd = 10 **(-4.4895)
    if element_name == 'Ne':
        abd = 10 **(-4.0701)
    if element_name == 'Mg':
        abd = 10 **(-4.4001)

    return abd

def get_the_factions(filename, ion, metallicity):

    element_name = ion.split('+')[-1]
    abd= get_gass10_abundances(element_name)

    data  = tab.Table.read(filename)
    tot_col = (data['H']+ data['H+'])*metallicity*abd
    frac_array = data[ion]/tot_col
    den_array  = data['hden']

    return den_array, frac_array



path = '/home/vikram/cloudy_run/metal_NH15_new'

uvb = ['KS18', 'HM12',  'P19', 'FG20']
uvb_Q = [14, 15, 16, 17, 18, 19, 20]

uvb_models =[]
the_Q_values = []
for background in uvb:
    if background == 'KS18':
        for q in uvb_Q:
            uvb_models.append(background)
            the_Q_values.append(q)
    else:
        q = 18
        uvb_models.append(background)
        the_Q_values.append(q)


logZ = -1
fname = (logZ+4)*100

for uvb_name, q_name in zip(uvb_models, the_Q_values):
    filename = 'try_{}_Q{}_Z{:.0f}.fits'.format(uvb_name, q_name, fname)
    file_with_path = path + '/' + filename
    find_max_o6_frac(filename= file_with_path)

logT = 5.5
print('for hot gas')
path = '/home/vikram/cloudy_run/hybrid_NH15'
fname = logT * 100
for uvb_name, q_name in zip(uvb_models, the_Q_values):

    filename = 'try_{}_Q{}_logT{:.0f}.fits'.format(uvb_name, q_name, fname)
    file_with_path = path + '/' + filename
    find_max_o6_frac(filename= file_with_path)
