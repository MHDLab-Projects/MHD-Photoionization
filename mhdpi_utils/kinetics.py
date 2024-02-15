"""
module for handling kinetics (recombination, ionization) data
"""

from pint import Quantity
import numpy as np
import xarray as xr

def get_kinetics(ds_cfd):

    rho_number = ds_cfd['rho_number'].pint.to('1/cm^3')
    T = ds_cfd['T'].pint.to('K')
    # combining k4 from 
    # 1. Axford, S.D.T., and Hayhurst, A.N. (1997). Mass spectrometric sampling of negative ions from flames of hydrogen and oxygen: the kinetics of electron attachment and detachment in hot mixtures of H2O, O2, OH and HO2. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences 452, 1007–1033. 10.1098/rspa.1996.0051.

    # Weighted average of collison partners for O2 rate
    k4_species = {
        'N2' : Quantity(1e-31, 'cm^6/s'),
        'H2O': Quantity(8e-30, 'cm^6/s'),
        'CO2': Quantity(1.3e-31, 'cm^6/s'), 
        'O2' : Quantity(2e-30, 'cm^6/s'),
    }

    weighted_avg = 0

    for species, k4 in k4_species.items():
        weighted_avg += k4*ds_cfd[species]

    rxn_rates = {
        'K+':  Quantity(4e-24, 'K*cm^6/s')*(1/T)*rho_number,
        'OH': Quantity(3e-31, 'cm^6/s')*rho_number,
        'O2_A': Quantity(5e-31, 'cm^6/s')*rho_number,
        'O2_B': weighted_avg,
        'O2_C': Quantity(6e-34, 'cm^6/s')*rho_number,
        'H2O': Quantity(1.6e-6, 'cm^3/s')*np.exp(-(Quantity(36060, 'K')/T)),
    }


    ds_kr = xr.Dataset(rxn_rates).pint.to('cm^3/s')

    return ds_kr