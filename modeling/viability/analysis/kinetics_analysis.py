# %%
import os
import numpy as np
import xarray as xr
import xyzpy

import matplotlib.pyplot as plt
import pint_xarray
from pint import Quantity

plt.rcParams.update({
    "savefig.facecolor": 'white',
    "font.size": 11, 
    'savefig.dpi': 300, 
    'font.sans-serif': 'arial', 
    # 'figure.figsize': (4.6, 3)
})

from dotenv import load_dotenv
load_dotenv()
REPO_DIR = os.getenv('REPO_DIR')

from mhdpy.fileio.path import chdir_if_nb_render; chdir_if_nb_render()
from pi_paper_utils import abscs, noneq

cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')




canterapath = os.path.join(REPO_DIR, 'modeling','viability', 'dataset', 'output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(canterapath, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

for species in ds_TP_species_rho.data_vars:
    ds_TP_species_rho[species] = ds_TP_species_rho[species].pint.quantify('particle/ml')


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
krm_O2 = krb_O2*ds_TP_species_rho['O2']
krm_Kp = krb_all['K+']*ds_TP_species_rho['K+']

krm_H20 = krb_all['H2O']*ds_TP_species_rho['H2O']
krm_OH = krb_all['OH']*ds_TP_species_rho['OH']

krm_sum = krm_O2 + 2*krm_Kp + krm_H20 + krm_OH
krm_sum = krm_sum.pint.to('1/s')


# %%


ne_max = krm_sum/krb_all['K+']

ne_max = ne_max.sel(phi=0.8, Kwt=0.01, P=1e5)

ne_max.plot()

plt.yscale('log')

plt.ylabel('Max $\Delta n_e$ [1/ml]')

