
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

from pi_paper_utils import abscs, noneq, kinetics


physical_data_folder = 'physical_data'
physical_data_folder = os.path.abspath(physical_data_folder) # Cantera needs absolute path

canterapath = os.path.join('output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(canterapath, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})


#%%

# kth, krb = kinetics.calc_kth_krb_kelly(ds_TP_params, ds_TP_species)

from pi_paper_utils.kinetics import gen_ds_krb

e = Quantity(1.602e-19, 'coulomb')
mp = Quantity(1.672e-27, 'kg')
kb = Quantity(8.617e-5, 'eV/K')
Eion = Quantity(4.34, 'eV')

# 
Ts = xr.DataArray(ds_TP_params['T']).pint.quantify("K")

krb_all = gen_ds_krb(Ts, ds_TP_params['rhocm3'].pint.quantify("particle/ml"))

krb_O2 = krb_all['O2_A']
krm_O2 = krb_O2*ds_TP_species_rho['O2'].pint.quantify("particle/ml")
krm_O2 = krm_O2.pint.to('1/s')

krm_Kp = krb_all['K+']*ds_TP_species_rho['K+'].pint.quantify("particle/ml")
krm_Kp = krm_Kp.pint.to('1/s')

krb_Kp = krb_all['K+']

ne0 = ds_TP_species['e-']

#%%

# Keq from 1. Ashton, A.F., and Hayhurst, A.N. (1973). Kinetics of collisional ionization of alkali metal atoms and recombination of electrons with alkali metal ions in flames. Combustion and Flame 21, 69–75. 10.1016/0010-2180(73)90008-4.
# Manually adding temperature to units to match overall units of molecule/ml specified in paper
Keq_prefactor = Quantity(2.42e15, 'molecule/ml/K**1.5')
Keq = Keq_prefactor*(Ts**(3/2))*np.exp(-(Eion/(kb*Ts)))
Keq = xr.DataArray(Keq, coords={'T':Ts.pint.magnitude}, dims = 'T')

kth = Keq*krb_Kp
kth.name =  '$k_{th}$'
kth = kth.pint.to('1/s')

Gth = kth*ds_TP_species_rho['K'].pint.quantify("molecule/cm**3")
Gth.attrs = dict(units = '$\#/cm^3 s$', long_name = 'Thermal Generation Rate')

#Perfect Ionization

B_const =5e-4
B_hall= np.sqrt(3600/(ds_TP_params['mobility'])**2)*1e-4
B_hall.name = 'Bmax ion slip'

combos = {
    'P_in' : np.array([0, *np.logspace(-10,10,11)]),
}

constants = {
    'mue_cant': ds_TP_params['mobility']*10000,
    'u': 1e5,
    'eta': 1,
    'B': B_const
}

dss = []

constants_temp = constants.copy()
constants_temp['G_th'] = Gth.pint.dequantify()
constants_temp['krb'] = krb_Kp.pint.dequantify()
ds = xyzpy.Runner(noneq.calc_NE_all, constants = constants_temp, var_names=None).run_combos(combos)
dss.append(ds.assign_coords(rxn='Kp').assign_coords(eta='perf'))

constants_temp = constants.copy()
constants_temp['krm'] = krm_O2.pint.dequantify()
constants_temp['ne0'] = ne0.pint.dequantify()
ds = xyzpy.Runner(noneq.calc_NE_all_const_nx, constants = constants_temp, var_names=None).run_combos(combos)
dss.append(ds.assign_coords(rxn='O2').assign_coords(eta='perf'))

constants_temp = constants.copy()
constants_temp['krm'] = krm_O2.pint.dequantify() + 2*krm_Kp.pint.dequantify()
constants_temp['ne0'] = ne0.pint.dequantify()
ds = xyzpy.Runner(noneq.calc_NE_all_const_nx, constants = constants_temp, var_names=None).run_combos(combos)
dss.append(ds.assign_coords(rxn='mm_sum').assign_coords(eta='perf'))


# #Photoionization
ds_cs = abscs.calc_ds_cs(ds_TP_species_rho.coords['T'].values, wls= [248]).squeeze()
gas_lam = noneq.calc_atten_lengths(ds_cs, ds_TP_species_rho)

# FA = xr.DataArray([1], name='FA_1')
FA = gas_lam['KOH'].where(False).fillna(1.0)
FA.name = 'FA_1'

eta_PI = noneq.calc_eta_PI(gas_lam['tot'], gas_lam['KOH'], FA)

constants_temp = constants.copy()
constants_temp['eta'] = eta_PI
constants_temp['krm'] = krm_O2.pint.dequantify() + 2*krm_Kp.pint.dequantify()
constants_temp['ne0'] = ne0.pint.dequantify()
ds = xyzpy.Runner(noneq.calc_NE_all_const_nx, constants = constants_temp, var_names=None).run_combos(combos)
dss.append(ds.assign_coords(rxn='mm_sum').assign_coords(eta='eta_PI'))


constants_temp = constants.copy()
constants_temp['eta'] = eta_PI*0.3
constants_temp['krm'] = krm_O2.pint.dequantify() + 2*krm_Kp.pint.dequantify()
constants_temp['ne0'] = ne0.pint.dequantify()
ds = xyzpy.Runner(noneq.calc_NE_all_const_nx, constants = constants_temp, var_names=None).run_combos(combos)
dss.append(ds.assign_coords(rxn='mm_sum').assign_coords(eta='eta_PI_QY'))

eta_PI

#%%

dss_new = []

for ds in dss:
    ds_new = ds.expand_dims('rxn').expand_dims('eta')
    ds_new = ds_new.stack(temp=['eta', 'rxn'])
    dss_new.append(ds_new)

ds_NE = xr.concat(dss_new, 'temp').unstack('temp')

# ds_NE = xr.concat([ds_NE_Kp, ds_NE_O2, ds_NE_sum,], 'rxn')

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
