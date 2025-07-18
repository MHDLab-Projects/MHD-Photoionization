"""
module for handling kinetics (recombination, ionization) data

krt - termolecular recombination rate (ml^2/particle^2/s)
krb - bimolecular recombination rate  (ml/particle/s)
krm - monomolecular recombination rate (1/s)
"""

from pint import Quantity
import numpy as np
import xarray as xr
import os
import pandas as pd


def calc_krbO2_weighted(ds_species):
    # combining krt from 
    # 1. Axford, S.D.T., and Hayhurst, A.N. (1997). Mass spectrometric sampling of negative ions from flames of hydrogen and oxygen: the kinetics of electron attachment and detachment in hot mixtures of H2O, O2, OH and HO2. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences 452, 1007–1033. 10.1098/rspa.1996.0051.

    # Weighted average of collison partners for O2 rate
    k4_species = {
        'N2' : Quantity(1e-31, 'ml^2/particle^2/s'),
        'H2O': Quantity(8e-30, 'ml^2/particle^2/s'),
        'CO2': Quantity(1.3e-31, 'ml^2/particle^2/s'), 
        'O2' : Quantity(2e-30, 'ml^2/particle^2/s'),
    }

    weighted_sum = 0

    for species, krt in k4_species.items():
        weighted_sum += krt*ds_species[species].pint.to('particle/ml')

    return weighted_sum.pint.to('ml/particle/s')


def gen_ds_krb(da_Ts, da_rho_number):

    da_rho_number = da_rho_number.pint.to('particle/ml')
    da_Ts = da_Ts.pint.to('K')


    krb_dict = {
        'K+':  Quantity(4e-24, 'K*ml^2/particle^2/s')*(1/da_Ts)*da_rho_number,
        'OH': Quantity(3e-31, 'ml^2/particle^2/s')*da_rho_number,
        'O2_A': Quantity(5e-31, 'ml^2/particle^2/s')*da_rho_number,
        'O2_G': Quantity(6e-34, 'ml^2/particle^2/s')*da_rho_number,
        'O2_exp_eff': Quantity(3.53e-32, 'ml^2/particle^2/s')*da_rho_number, # This comes from effective termolecular rate calculated from experimental data (fig:SI_krb_eff_53x)
        'H2O': Quantity(1.6e-6, 'ml/particle/s')*np.exp(-(Quantity(36060, 'K')/da_Ts)),
    }

    ds_krb = xr.Dataset(krb_dict).pint.to('ml/particle/s')

    return ds_krb

def calc_krm(ds_krb, ds_species):
    """
    Calculate monomolecular recombination rate from bimolecular recombination rate and species concentrations
    """
    ds_krm = xr.Dataset()
    for species in ds_krb.data_vars:
        species_name = species.split('_')[0]
        ds_krm[species] = ds_krb[species].pint.to('ml/particle/s')*ds_species[species_name].pint.to('particle/ml')
    
    ds_krm = ds_krm.pint.to('1/s')

    return ds_krm


#Calculate potassium thermal kinetics

# # Old method: 1. Ashton, A.F., and Hayhurst, A.N. (1973). Kinetics of collisional ionization of alkali metal atoms and recombination of electrons with alkali metal ions in flames. Combustion and Flame 21, 69–75. 10.1016/0010-2180(73)90008-4.
# Direct calcualtion of kth and kr from Ashton and Hayhurst, for H2, O2, N2 flames. 

# def calc_kth(rho):
#     T = rho.coords['T'].values.item()
#     return 1e-8*(np.sqrt(T))*np.exp(-(4.34)/(8.617e-5*T))*rho

# def calc_krb(rho):
#     T = rho.coords['T'].values.item()
#     return 4e-24*(1/T)*rho

# kth = rho.groupby('T').apply(calc_kth)
# krb = rho.groupby('T').apply(calc_krb)


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

def calc_kth_krb_kelly(ds_TP_params, ds_TP_species, cs=cs, redmass=redmass):
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

    krb = kth/Keq
    krb.name =  '$k_{r,b}$'
    krb = krb.pint.to('cm**3/molecule/s')

    return kth, krb