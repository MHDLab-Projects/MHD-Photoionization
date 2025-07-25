
"""Calculate curves for gamma =1 with linear interpolation of each pressure vs temperature curve"""
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

for species in ds_TP_species_rho.data_vars:
    ds_TP_species_rho[species] = ds_TP_species_rho[species].pint.quantify('particle/ml')

#%%
    
ds_TP_params['sigma'].sel(Kwt=0.01, phi=0.8, P=1e5, method='nearest').plot()

#%%

# kth, krb = kinetics.calc_kth_krb_kelly(ds_TP_params, ds_TP_species)

from pi_paper_utils.kinetics import gen_ds_krb, calc_krm

e = Quantity(1.602e-19, 'coulomb')
mp = Quantity(1.672e-27, 'kg')
kb = Quantity(8.617e-5, 'eV/K')
Eion = Quantity(4.34, 'eV')

# 
Ts = xr.DataArray(ds_TP_params['T']).pint.quantify("K")

krb_all = gen_ds_krb(Ts, ds_TP_params['rhocm3'].pint.quantify("particle/ml"))


krm = calc_krm(krb_all, ds_TP_species_rho)

krm_sum = krm['O2_exp_eff'] + 2*krm['K+'] + krm['H2O'] + krm['OH']
krm_sum = krm_sum.pint.to('1/s')
krm['mm_sum'] = krm_sum

#%%
ne0 = ds_TP_species_rho['e-']
krb_Kp = krb_all['K+']

# Keq from 1. Ashton, A.F., and Hayhurst, A.N. (1973). Kinetics of collisional ionization of alkali metal atoms and recombination of electrons with alkali metal ions in flames. Combustion and Flame 21, 69–75. 10.1016/0010-2180(73)90008-4.
# Manually adding temperature to units to match overall units of molecule/ml specified in paper
Keq_prefactor = Quantity(2.42e15, 'molecule/ml/K**1.5')
Keq = Keq_prefactor*(Ts**(3/2))*np.exp(-(Eion/(kb*Ts)))
Keq = xr.DataArray(Keq, coords={'T':Ts.pint.magnitude}, dims = 'T')

kth = Keq*krb_Kp
kth.name =  '$k_{th}$'
kth = kth.pint.to('1/s')

Gth = kth*ds_TP_species_rho['K']

#Perfect Ionization

# #Photoionization
ds_cs = abscs.calc_ds_cs(ds_TP_species_rho.coords['T'].values, wls= [248]).squeeze()
ds_cs = ds_cs.pint.quantify()
gas_lam = noneq.calc_atten_lengths(ds_cs, ds_TP_species_rho)

# FA = xr.DataArray([1], name='FA_1')
FA = gas_lam['KOH'].pint.dequantify().where(False).fillna(1.0).pint.quantify('dimensionless')
FA.name = 'FA_1'

eta_PI = noneq.calc_eta_PI(gas_lam['tot'], gas_lam['KOH'], FA).pint.to('dimensionless')

#TODO: is there a way combos can be passed in quantified?
combos = {
    # 'P_in' : np.array([0, *np.logspace(-10,10,11)]),
    # 'P_in' : [Quantity(0, 'W/cm**3').magnitude],
    'P_in' : [0],
}

constants = {
    'mue_cant': ds_TP_params['mobility'].pint.quantify(),
    'u': Quantity(1e5, 'cm/s'),
    'eta': Quantity(1, 'dimensionless'),
    'B': Quantity(5,'T')
}


# constants_temp = constants.copy()
# constants_temp['G_th'] = Gth.pint.dequantify()
# constants_temp['krb'] = krb_Kp.pint.dequantify()
# ds = xyzpy.Runner(noneq.calc_NE_all, constants = constants_temp, var_names=None).run_combos(combos)
# dss.append(ds.assign_coords(rxn='K+_bm').assign_coords(eta='perf'))

#TODO: incorportate these for loops into xyzpy runner, 
# also readd K+ bimolecular, but should be consistent wiht K+ monomolecular for P_in=0

dss = []

eta_dict = {'perf': 1, 'KOH': eta_PI}

for eta_str, eta in eta_dict.items():
    for krm_val in ['mm_sum', 'O2_A', 'O2_exp_eff', 'K+', 'H2O', 'OH']:
        constants_temp = constants.copy()
        constants_temp['eta'] = eta
        constants_temp['ne0'] = ne0
        constants_temp['krm'] = krm[krm_val]
        ds = xyzpy.Runner(noneq.calc_NE_all_const_nx, constants = constants_temp, var_names=None).run_combos(combos)
        dss.append(ds.assign_coords(rxn=krm_val).assign_coords(eta=eta_str))


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

# flow = bulk flow
# bl = boundary layer

sig_flow = ds_TP_params['sigma'].sel(T=3000, method='nearest').pint.quantify()
sig_bl = ds_NE['sigma'].pint.quantify()

combos_l_flow = {'l_flow' :  [0, 0.5, 0.9, 0.99]}

r = xyzpy.Runner(noneq.calc_dsigma_tot, constants = {'sigma_flow' :sig_flow, 'sigma_bl': sig_bl }, var_names = ['sigma_tot'])
da_dsigma_tot = r.run_combos(combos_l_flow).squeeze()
da_dsigma_tot.name = 'enhancement factor'
da_dsigma_tot.attrs = {}

da_dsigma_tot.pint.dequantify().to_netcdf('output/da_dsigma_tot.cdf')


