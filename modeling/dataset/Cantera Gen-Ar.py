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

import mhdpy.cantera.et
import mhdpy.cantera.util as ct_utils

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.size': 16})

## Jupyter VS plugin Assumes VS code directory is MHDLab folder
# get_ipython().run_line_magic('matplotlib', 'inline')
# folder = os.path.join(os.getcwd(), "Analysis\\Photoionization\\Simulation\\Cantera Photoionization")

## Pure Python Script
folder = os.path.split(os.path.realpath(__file__))[0]  #pure python script

# Jupyter notebook
# folder = os.getcwd() 

import sys
if not folder in sys.path:
    sys.path.append(folder)


# In[2]:


#Setup cantera database and define initial gas object


outputtype = 'X'
paramnames = ['rho','enthalpy','Cp','Cv','sigma','mobility','viscosity','thermalConduct']

gas_tr, electron_trans = ct_utils.gen_gas(
                    GasThermo_tr = 'data/drm19_SeededKero.cti',
                    GasThermoName_tr = 'gas',
                    CrossSectionFile = 'data/NETL_CF2017.csv',
                    basis = 'molar' #why is this molar?
                    )

nTP = 9
combos_TP = {
    'T' : np.linspace(1000,3000,nTP),
    'P' : np.logspace(3,7,nTP),
    'Kwt': np.logspace(-9,-1,9),
}

setK = combos_TP['Kwt'][0]
gas_tr.TPY =300,1e5, {'Ar': 1- setK, 'K': setK}
gas_tr.equilibrate("HP") #not sure if this is nessecary
gas_tr.equilibrate("TP") #not sure if this is nessecary

def calc_TP(gas_tr, T, P, Kwt):
        
    setK = Kwt
    gas_tr.TPX = T, P, {'Ar': 1-setK, 'K': setK}
    gas_tr.equilibrate('TP')

    params, dat_species = mhdpy.cantera.util.extract_params(gas_tr, electron_trans, outputtype,)

    return tuple([*params, *dat_species])   


r = xyzpy.Runner(calc_TP,
    constants = {'gas_tr':gas_tr},
    var_names = [*paramnames,*gas_tr.species_names])
ds = r.run_combos(combos_TP)

ds.coords['Kwt'].attrs = dict(units = '%', long_name = 'K weight percent')
ds.coords['T'].attrs = dict(units = 'K', long_name = 'Temperature')
ds.coords['P'].attrs = dict(units = 'Pa', long_name = 'Pressure')
ds.attrs = dict()
ds = ds.squeeze()


ds_TP_params, ds_TP_species, ds_TP_species_rho = mhdpy.cantera.util.process_ds_speciesparams(ds, gas_tr)

#drop all empty species
da = ds_TP_species.to_array('species')
da = da.where(da>0).dropna('species','all')
ds_TP_species = da.to_dataset('species')

# %%

#Calculate potassium thermal kinetics


e = 1.602e-19
mp = 1.672e-27 
kb = 8.617e-5
Eion = 4.34

s = pd.read_csv(os.path.join(folder,'IonCrossSection.csv'), index_col=(0,1))['cs']
cs = xr.DataArray.from_series(s)
cs = cs * 1e-14 #nm^2 to cm^2

cs = cs.sel(metal = 'K').drop('metal')
cs = cs.rename(gas='species')

masses = {
    'K': 39.098,
    'Ar': 39.948,
    'CO': 28.01,
    'CO2': 44,
    'H2': 2.01,
    'H2O':18.015,
    'N2': 28.01
}

for species in cs.indexes['species']:
    if species not in ds_TP_species.data_vars:
        cs = cs.drop(species, 'species')
        del masses[species]

redmass = cs.copy(deep=True)

for species in redmass.indexes['species']:
    redmass.loc[{'species':species}] = (masses['K']*masses[species])/(masses['K'] + masses[species])
    
redmass = redmass*mp

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

rhocm3 = ds_TP_params['rhocm3']

Ts = ds_TP_species.indexes['T'].values
A = xr.DataArray(calc_A(Ts), coords={'T':Ts}, dims = 'T')
B = calc_B(ds_TP_species, cs, redmass)

Keq = 2.42e15*(Ts**(3/2))*np.exp(-(Eion/(kb*Ts)))
Keq = xr.DataArray(Keq, coords={'T':Ts}, dims = 'T')

kth = rhocm3*A*B
kr = kth/Keq

kr.name =  '$k_{r}$'
kth.name = '$k_{th}$'
kr = kr.assign_attrs(units = '$(cm^3/molecule)(1/s)$' )
kth = kth.assign_attrs(units = '$(1/s)$' )
 
ds_TP_params = ds_TP_params.assign({'kth':kth, 'kr':kr})

def calc_Gth(kth,m):
    return kth*m

Gth = calc_Gth(ds_TP_params['kth'],ds_TP_species_rho['K'])

Gth.attrs = dict(units = '$\#/cm^3 s$', long_name = 'Thermal Generation Rate')

ds_TP_params = ds_TP_params.assign({'Gth':Gth})



# In[18]:

pathout = mhdpy.io.gen_path('sharepoint', 'lee','Data','Cantera Gasses','Ar')

# df_input.to_csv(os.path.join(path, 'inputs.csv'))

ds_TP_params.to_netcdf(os.path.join(pathout, 'ds_TP_params.cdf'))
ds_TP_species.to_netcdf(os.path.join(pathout, 'ds_TP_species.cdf'))
ds_TP_species_rho.to_netcdf(os.path.join(pathout, 'ds_TP_species_rho.cdf'))



ds_TP_params.to_dataframe().to_csv(os.path.join(pathout, 'ds_TP_params.csv'))
mhdpy.analysis.xr.writeunitsrowcsv(ds_TP_params,os.path.join(pathout, 'ds_TP_params.csv'))

ds_TP_species.to_dataframe().to_csv(os.path.join(pathout, 'ds_TP_species.csv'))
mhdpy.analysis.xr.writeunitsrowcsv(ds_TP_species,os.path.join(pathout, 'ds_TP_species.csv'))

ds_TP_species_rho.to_dataframe().to_csv(os.path.join(pathout, 'ds_TP_species_rho.csv'))
mhdpy.analysis.xr.writeunitsrowcsv(ds_TP_species_rho,os.path.join(pathout, 'ds_TP_species_rho.csv'))
