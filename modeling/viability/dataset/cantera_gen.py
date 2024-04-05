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
    'phi':[0.8,0.9,1.0],
    'Kwt':np.logspace(-6,-1,6),
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

#%%

ds_TP_species = ds_TP_species.drop('C12H26')
ds_TP_species_rho = ds_TP_species_rho.drop('C12H26')

ds_HP_species = ds_HP_species.drop('C12H26')
ds_HP_species_rho = ds_HP_species_rho.drop('C12H26')

# In[18]:

# # kth calculation above converts temperature coordinate attrs to be Unit('K') which can't be written to cdf.
# # quick fix, not sure if this is the best way to do it
# for coord in ['T', 'P']:
#     ds_TP_params.coords[coord] = ds_TP_params.coords[coord].pint.dequantify()


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
