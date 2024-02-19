
"""Calculate curves for alpha =1 with linear interpolation of each pressure vs temperature curve"""
# %%
import os
import numpy as np
import pandas as pd
import xarray as xr
import xyzpy
from pint import Quantity
import pint_xarray

from dotenv import load_dotenv
load_dotenv()

from pi_paper_utils import abscs, noneq


physical_data_folder = 'physical_data'
physical_data_folder = os.path.abspath(physical_data_folder) # Cantera needs absolute path

canterapath = os.path.join('output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(canterapath, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})


#%%

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

ds_TP_params = ds_TP_params.pint.dequantify()

#Perfect Ionization

B_const =5e-4
B_hall= np.sqrt(3600/(ds_TP_params['mobility'])**2)*1e-4
B_hall.name = 'Bmax ion slip'

combos = {
    'P_in' : np.array([0, *np.logspace(-10,10,11)]),
}

constants = {
    'G_th': ds_TP_params['Gth'],
    'kr': ds_TP_params['kr'],
    'mue_cant': ds_TP_params['mobility']*10000,
    'u': 1e5
}

constants['eta'] = 1
constants['B'] = B_const
ds_NE_perf_Bconst = xyzpy.Runner(noneq.calc_NE_all, constants = constants, var_names=None).run_combos(combos).assign_coords(analysis='perf_Bconst')

constants['B'] = B_hall
ds_NE_perf_Bhall = xyzpy.Runner(noneq.calc_NE_all, constants = constants, var_names=None).run_combos(combos).assign_coords(analysis='perf_Bhall')

# #Photoionization
ds_cs = abscs.calc_ds_cs(ds_TP_species_rho.coords['T'].values, wls= [248]).squeeze()
gas_lam = noneq.calc_atten_lengths(ds_cs, ds_TP_species_rho)
FA = xr.DataArray([1], name='FA_1')
eta_PI = noneq.calc_eta_PI(gas_lam['tot'], gas_lam['K'], FA)
constants['eta'] = eta_PI

constants['B'] = B_const
ds_NE_PI_Bconst  = xyzpy.Runner(noneq.calc_NE_all, constants = constants, var_names=None).run_combos(combos).assign_coords(analysis='PI_Bconst')

constants['B'] = B_hall
ds_NE_PI_Bhall  = xyzpy.Runner(noneq.calc_NE_all, constants = constants, var_names=None).run_combos(combos).assign_coords(analysis='PI_Bhall')

ds_NE = xr.concat([ds_NE_perf_Bconst, ds_NE_perf_Bhall, ds_NE_PI_Bconst, ds_NE_PI_Bhall], 'analysis')
# %%
ds_NE.attrs = {}

ds_NE.to_netcdf('output/ds_NE.cdf')
# %%

sig_bk = ds_TP_params['sigma'].sel(T=3000, method='nearest')
# sig_bk = ds_NE['sigma'].sel(T=3000, method='nearest') #This was used previously, incorrectly

combos_l_bk = {'l_bk' :  [0, 0.5, 0.9, 0.99]}

r = xyzpy.Runner(noneq.calc_dsigma_tot, constants = {'sigma_bk' :sig_bk, 'sigma_bl': ds_NE['sigma'] }, var_names = ['sigma_tot'])
da_dsigma_tot = r.run_combos(combos_l_bk).squeeze()
da_dsigma_tot.name = 'enhancement factor'
da_dsigma_tot.attrs = {}

da_dsigma_tot.to_netcdf('output/da_dsigma_tot.cdf')

# alpha_bl = ds_NE['alpha']*da_dsigma_tot

# alpha_bl.name = 'alpha_bl'
# alpha_bl.to_netcdf('output/alpha_bl.cdf')



# %%
