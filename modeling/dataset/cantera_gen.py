"""
Generation of a dataset of Cantera outputs for a range of conditions for oxyfuel combustion with seeded kerosene.
Initializes a mixture of C12H26, and O2 at fixed enthalpy and pressure, then varies the temperature and pressure to generate a dataset of species and parameters.
Finally a nonequilibrium model is applied to the dataset to calculate the thermal and recombination kinetics of potassium.
"""


#!/usr/bin/env python
# coding: utf-8

# In[1]:

import matplotlib.pyplot as plt
import pandas as pd
import cantera as ct
import numpy as np
import os
from collections import OrderedDict
import re
import xarray as xr
import xyzpy

from pint import Quantity
import pint_pandas
import pint_xarray

import mhdpy.cantera_utils.util as ct_utils 

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 16})

physical_data_folder = 'physical_data'
physical_data_folder = os.path.abspath(physical_data_folder) # Cantera needs absolute path

output_path = 'output'
if not os.path.exists(output_path): os.mkdir(output_path)

# In[2]:


#Setup cantera database and define initial gas object

paramnames = ['rho','enthalpy','Cp','Cv','sigma','mobility','viscosity','thermalConduct']



gas_tr, electron_trans = ct_utils.gen_gas(
                    GasThermo_tr = os.path.join(physical_data_folder, 'drm19_SeededKero.yaml'),
                    GasThermoName_tr = 'gas',
                    CrossSectionFile = os.path.join(physical_data_folder, 'NETL_CF2017.csv'),
                    basis = 'molar'
                    )

# In[4]:

#Define fuel/o2/seed combos to be applied and run over combinations



combos_HP = {
    'phi':[0.7,0.8,0.9,1.0,1.1,1.2,1.3],
    'Kwt':np.logspace(-4,-1,7),
    'inlet_T': np.linspace(300,1000,3),
    'P_combustor': np.logspace(5,7,3),
}

r = xyzpy.Runner(ct_utils.calc_combustor,
    constants = 
        {'gas_tr':gas_tr,
        'electron_trans':electron_trans,
        'fueltype' :'C12H26'},
    var_names = ['T_adiabatic', *paramnames,*gas_tr.species_names])
ds = r.run_combos(combos_HP)

ds.coords['Kwt'].attrs = dict(long_name = 'K mass fraction')
ds.coords['phi'].attrs = dict(long_name = 'Equivalence Ratio')
ds.coords['inlet_T'].attrs = dict(units = 'K', long_name = 'Inlet Temperature')
ds.coords['P_combustor'].attrs = dict(units = 'Pa', long_name = 'Combustor Pressure')
ds.attrs = dict() # Gets rid of attributes added by xyzpy that cannot be written to cdf
ds_HP_params, ds_HP_species, ds_HP_species_rho = ct_utils.process_ds_speciesparams(ds, gas_tr)

# %%

#Define TP combos to be applied and run over combinations. fuel and Kwt values can be downselected here. 

nT = 41
nP = 41
combos_TP = {
    'T' : np.linspace(1000,4000,nT),
    'P' : np.logspace(3,7,nP),
    'phi':combos_HP['phi'],
    'Kwt': combos_HP['Kwt'],
    'inlet_T': [300],
    'P_combustor': [1e6],
}


r = xyzpy.Runner(ct_utils.calc_TP,
    constants = 
        {'gas_tr':gas_tr,
        'electron_trans':electron_trans,
        'ds_HP_species': ds_HP_species
        },
    var_names = [*paramnames,*gas_tr.species_names])
ds = r.run_combos(combos_TP)

ds.coords['Kwt'].attrs = dict(long_name = 'K mass fraction')
ds.coords['phi'].attrs = dict(long_name = 'Equivalence Ratio')
ds.coords['T'].attrs = dict(units = 'K', long_name = 'Temperature')
ds.coords['P'].attrs = dict(units = 'Pa', long_name = 'Pressure')
ds.attrs = dict()
ds = ds.squeeze()

ds_TP_params, ds_TP_species, ds_TP_species_rho = ct_utils.process_ds_speciesparams(ds, gas_tr)

# %%

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


# # New method: 1. Kelly, R., and Padley, P.J. (1972). Measurement of collisional ionization cross-sections for metal atoms in flames. Proc. R. Soc. Lond. A 327, 345–366. 10.1098/rspa.1972.0050.
# Calculate kth from Kelly and Padley, using ionization cross sections for different species, then calculate kr using saha Keq from Ashton and Hayhurst.

#Define reaction as kf = A(T)*B(X_i)

e = Quantity(1.602e-19, 'coulomb')
mp = Quantity(1.672e-27, 'kg')
kb = Quantity(8.617e-5, 'eV/K')
Eion = Quantity(4.34, 'eV')

# Ion cross sections 
# 
# Table 4 p 361
s = pd.read_csv(os.path.join(physical_data_folder,'IonCrossSection.csv'), index_col=(0,1))['cs']
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

ds_TP_params = ds_TP_params.assign({'kth':kth, 'kr':kr})
Gth = ds_TP_params['kth']*ds_TP_species_rho['K'].pint.quantify("molecule/cm**3")
Gth.attrs = dict(units = '$\#/cm^3 s$', long_name = 'Thermal Generation Rate')
ds_TP_params = ds_TP_params.assign({'Gth':Gth})

#%%

ds_TP_species = ds_TP_species.drop('C12H26')
ds_TP_species_rho = ds_TP_species_rho.drop('C12H26')

ds_HP_species = ds_HP_species.drop('C12H26')
ds_HP_species_rho = ds_HP_species_rho.drop('C12H26')

# In[18]:

# kth calculation above converts temperature coordinate attrs to be Unit('K') which can't be written to cdf.
# quick fix, not sure if this is the best way to do it
for coord in ['T', 'P']:
    ds_TP_params.coords[coord] = ds_TP_params.coords[coord].pint.dequantify()


ds_dict = {
    'ds_TP_params': ds_TP_params,
    'ds_TP_species': ds_TP_species,
    'ds_TP_species_rho': ds_TP_species_rho,
    'ds_HP_params': ds_HP_params,
    'ds_HP_species': ds_HP_species,
    'ds_HP_species_rho': ds_HP_species_rho

}

for key, ds in ds_dict.items():
    ds = ds.pint.dequantify()
    ds.to_netcdf(os.path.join(output_path, f'{key}.cdf'))
