
# %%
import os
import numpy as np
import xarray as xr
import xyzpy

import matplotlib.pyplot as plt

plt.rcParams.update({
    "savefig.facecolor": 'white',
    "font.size": 11, 
    'savefig.dpi': 300, 
    # 'font.sans-serif': 'arial', 
    # 'figure.figsize': (4.6, 3)
})

from dotenv import load_dotenv
load_dotenv()
REPO_DIR = os.getenv('REPO_DIR')

from pi_paper_utils import abscs, noneq

cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_P_zero = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'P_zero.cdf'))

# %%

ds_NE = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE.cdf')).squeeze()
alpha = ds_NE['alpha']

da_dsigma_tot = xr.load_dataset(os.path.join(PI_modeling_dataset_dir, 'da_dsigma_tot.cdf'))['enhancement factor']
alpha = alpha*da_dsigma_tot



# %%

plt.figure()

combo_downsel = {
    'P_in' : 0,
    'Kwt': [0.001, 0.01, 0.1],
    'phi': 0.8,
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero.plot(col='l_bk', row='Kwt',hue='analysis', y='T', xscale='log')
# P_zero.plot(row='Kwt',hue='analysis', y='T', xscale='log')

plt.savefig('output/beta_boundary_layer.png')
# %%