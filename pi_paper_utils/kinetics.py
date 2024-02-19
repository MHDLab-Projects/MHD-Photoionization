"""
module for handling kinetics (recombination, ionization) data
"""

from pint import Quantity
import numpy as np
import xarray as xr
import os
import pandas as pd

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



#Calculate potassium thermal kinetics

# # Old method: 1. Ashton, A.F., and Hayhurst, A.N. (1973). Kinetics of collisional ionization of alkali metal atoms and recombination of electrons with alkali metal ions in flames. Combustion and Flame 21, 69–75. 10.1016/0010-2180(73)90008-4.
# Direct calcualtion of kth and kr from Ashton and Hayhurst, for H2, O2, N2 flames. 

# def calc_kth(rho):
#     T = rho.coords['T'].values.item()
#     return 1e-8*(np.sqrt(T))*np.exp(-(4.34)/(8.617e-5*T))*rho

# def calc_kr(rho):
#     T = rho.coords['T'].values.item()
#     return 4e-24*(1/T)*rho

# kth = rho.groupby('T').apply(calc_kth)
# kr = rho.groupby('T').apply(calc_kr)


e = Quantity(1.602e-19, 'coulomb')
mp = Quantity(1.672e-27, 'kg')
kb = Quantity(8.617e-5, 'eV/K')
Eion = Quantity(4.34, 'eV')

# Ion cross sections 
# 
# Table 4 p 361
packagepath = os.path.dirname(os.path.realpath(__file__))
s = pd.read_csv(os.path.join(packagepath, 'kinetic_data', 'IonCrossSection.csv'), index_col=(0,1))['cs']
cs = xr.DataArray.from_series(s) 
cs = cs.pint.quantify('nm**2/particle').pint.to('cm**2/particle')
cs = cs.sel(metal = 'K').drop('metal') 

masses = {
    'K': 39.098,
    'Ar': 39.948,
    'CO': 28.01,
    'CO2': 44,
    'H2': 2.01,
    'H2O':18.015,
    'N2': 28.01
}

masses = {
    k: Quantity(v, 'amu') for k,v in masses.items()
}

redmass = { k: (masses['K']*v)/(masses['K']+v) for k,v in masses.items()}

def calc_A(T):
    A = np.sqrt(8*kb*T*np.pi)*np.exp(-Eion/(kb*T))
    # A = A*1e2 #m to cm (units are weird since mass has been pulled into B, but should end up like m/s per kg**1/2 or something)
    return A

def calc_B(ds_species, cs, redmass):
    b_tot = 0
    for species in cs.indexes['species']:
        x_i = ds_species[species]
        cs_i = cs.sel(species=species)
        redmass_i = redmass[species]
        b_i = x_i*cs_i/np.sqrt(redmass_i)
        b_tot = b_tot + b_i
    return b_tot

def calc_kth_kr_kelly(ds_TP_params, ds_TP_species, cs=cs, redmass=redmass):
    # # New method: 1. Kelly, R., and Padley, P.J. (1972). Measurement of collisional ionization cross-sections for metal atoms in flames. Proc. R. Soc. Lond. A 327, 345–366. 10.1098/rspa.1972.0050.
    # Calculate kth from Kelly and Padley, using ionization cross sections for different species, then calculate kr using saha Keq from Ashton and Hayhurst.

    #Define reaction as kf = A(T)*B(X_i)

    Ts = Quantity(ds_TP_species.indexes['T'].values, 'K')
    A = xr.DataArray(calc_A(Ts), coords={'T':Ts.magnitude}, dims = 'T')
    B = calc_B(ds_TP_species, cs, redmass)

    kth = ds_TP_params['rhocm3'].pint.quantify("particle/cm**3")*A*B
    kth.name = '$k_{th}$'
    kth = kth.pint.to('1/s')

    # Keq from 1. Ashton, A.F., and Hayhurst, A.N. (1973). Kinetics of collisional ionization of alkali metal atoms and recombination of electrons with alkali metal ions in flames. Combustion and Flame 21, 69–75. 10.1016/0010-2180(73)90008-4.
    # Manually adding temperature to units to match overall units of molecule/ml specified in paper
    Keq_prefactor = Quantity(2.42e15, 'molecule/ml/K**1.5')
    Keq = Keq_prefactor*(Ts**(3/2))*np.exp(-(Eion/(kb*Ts)))
    Keq = xr.DataArray(Keq, coords={'T':Ts.magnitude}, dims = 'T')

    kr = kth/Keq
    kr.name =  '$k_{r}$'
    kr = kr.pint.to('cm**3/molecule/s')

    return kth, kr