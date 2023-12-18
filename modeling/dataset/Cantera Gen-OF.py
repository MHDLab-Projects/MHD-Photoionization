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


# import cantera_utils.et
# import cantera_utils.util as ct_utils

import mhdpy.cantera.util as ct_utils 

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 16})

physical_data_folder = 'physical_data'
physical_data_folder = os.path.abspath(physical_data_folder) # Cantera needs absolute path

# from mhdpy.fp import gen_path
# output_path = gen_path('sharepoint', 'lee','Data','Cantera Gasses')

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

# def calc_kth(rho):
#     T = rho.coords['T'].values.item()
#     return 1e-8*(np.sqrt(T))*np.exp(-(4.34)/(8.617e-5*T))*rho

# def calc_kr(rho):
#     T = rho.coords['T'].values.item()
#     return 4e-24*(1/T)*rho

# kth = rho.groupby('T').apply(calc_kth)
# kr = rho.groupby('T').apply(calc_kr)


#Define reaction as kf = A(T)*B(X_i)


e = 1.602e-19
mp = 1.672e-27 
kb = 8.617e-5
Eion = 4.34

s = pd.read_csv(os.path.join(physical_data_folder,'IonCrossSection.csv'), index_col=(0,1))['cs']
cs = xr.DataArray.from_series(s) * 1e-14 #nm^2 to cm^2
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

redmass = cs.copy(deep=True).rename('redmass')

for species in redmass.indexes['species']:
    redmass.loc[{'species':species}] = ( (masses['K']*masses[species])/(masses['K'] + masses[species]) )*mp
    
def calc_A(T):
    A = np.sqrt(8*e*kb*T*np.pi)*np.exp(-Eion/(kb*T))
    A = A*1e2 #m to cm (units are weird since mass has been pulled into B, but should end up like m/s per kg**1/2 or something)
    return A

def calc_B(ds_species, cs, redmass):
    b_tot = 0
    for species in cs.indexes['species']:
        x_i = ds_species[species]
        cs_i = cs.sel(species=species)
        redmass_i = redmass.sel(species=species)
        b_i = x_i*cs_i/np.sqrt(redmass_i)
        b_tot = b_tot + b_i
    return b_tot

Ts = ds_TP_species.indexes['T'].values
A = xr.DataArray(calc_A(Ts), coords={'T':Ts}, dims = 'T')
B = calc_B(ds_TP_species, cs, redmass)

kth = ds_TP_params['rhocm3']*A*B
kth.name = '$k_{th}$'
kth = kth.assign_attrs(units = '$(1/s)$' )

Keq = 2.42e15*(Ts**(3/2))*np.exp(-(Eion/(kb*Ts)))
Keq = xr.DataArray(Keq, coords={'T':Ts}, dims = 'T')

kr = kth/Keq
kr.name =  '$k_{r}$'
kr = kr.assign_attrs(units = '$(cm^3/molecule)(1/s)$' )

ds_TP_params = ds_TP_params.assign({'kth':kth, 'kr':kr})
Gth = ds_TP_params['kth']*ds_TP_species_rho['K']
Gth.attrs = dict(units = '$\#/cm^3 s$', long_name = 'Thermal Generation Rate')
ds_TP_params = ds_TP_params.assign({'Gth':Gth})


# %% 

ds_TP_species = ds_TP_species.drop('C12H26')
ds_TP_species_rho = ds_TP_species_rho.drop('C12H26')

ds_HP_species = ds_HP_species.drop('C12H26')
ds_HP_species_rho = ds_HP_species_rho.drop('C12H26')


# In[18]:



# df_input.to_csv(os.path.join(path, 'inputs.csv'))

ds_TP_params.to_netcdf(os.path.join(output_path, 'ds_TP_params.cdf'))
ds_TP_species.to_netcdf(os.path.join(output_path, 'ds_TP_species.cdf'))
ds_TP_species_rho.to_netcdf(os.path.join(output_path, 'ds_TP_species_rho.cdf'))


ds_HP_params.to_netcdf(os.path.join(output_path, 'ds_HP_params.cdf'))
ds_HP_species.to_netcdf(os.path.join(output_path, 'ds_HP_species.cdf'))
ds_HP_species_rho.to_netcdf(os.path.join(output_path, 'ds_HP_species_rho.cdf'))


# from mhdpy.analysis.xr import writeunitsrowcsv

# ds_TP_params.to_dataframe().to_csv(os.path.join(output_path, 'ds_TP_params.csv'))
# writeunitsrowcsv(ds_TP_params,os.path.join(output_path, 'ds_TP_params.csv'))

# ds_TP_species.to_dataframe().to_csv(os.path.join(output_path, 'ds_TP_species.csv'))
# writeunitsrowcsv(ds_TP_species,os.path.join(output_path, 'ds_TP_species.csv'))

# ds_TP_species_rho.to_dataframe().to_csv(os.path.join(output_path, 'ds_TP_species_rho.csv'))
# writeunitsrowcsv(ds_TP_species_rho,os.path.join(output_path, 'ds_TP_species_rho.csv'))

# ds_HP_params.to_dataframe().to_csv(os.path.join(output_path, 'ds_HP_params.csv'))
# writeunitsrowcsv(ds_HP_params,os.path.join(output_path, 'ds_HP_params.csv'))

# ds_HP_species.to_dataframe().to_csv(os.path.join(output_path, 'ds_HP_species.csv'))
# writeunitsrowcsv(ds_HP_species,os.path.join(output_path, 'ds_HP_species.csv'))

# ds_HP_species_rho.to_dataframe().to_csv(os.path.join(output_path, 'ds_HP_species_rho.csv'))
# writeunitsrowcsv(ds_HP_species_rho,os.path.join(output_path, 'ds_HP_species_rho.csv'))
