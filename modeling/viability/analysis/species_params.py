#%%

# %%
import os
import numpy as np
import xarray as xr

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

cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')
# %%

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

# gamma = xr.open_dataset(os.path.join(PI_modeling_dataset_dir,'gamma_bl.cdf')).squeeze()['gamma_bl']
ds_NE_Kp = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE_Kp.cdf')).squeeze()
ds_NE_O2 = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE_O2.cdf')).squeeze()

#%%
# make a dictionary of 5 pressures and temperatures from ds_NE coordinates

# Generate 5 evenly spaced indices for 'P' and 'T'
indices_P = np.linspace(0, len(ds_NE_Kp.coords['P'].values) - 1, 5, dtype=int)
indices_T = np.linspace(0, len(ds_NE_Kp.coords['T'].values) - 1, 5, dtype=int)

# Use these indices to select the corresponding values from your coordinates
seldict = {
    'P': ds_NE_Kp.coords['P'].values[indices_P],
    'T': ds_NE_Kp.coords['T'].values[indices_T],
    'Kwt': 0.01
}

# %%
ds_TP_species['O2'].sel(Kwt=0.01).plot(col='phi', col_wrap=3)

# %%
ds_TP_species['O2'].sel(seldict).plot(col='phi', col_wrap=3, hue='P')

plt.yscale('log')
#%%


ds_NE_O2['krm'].sel(Kwt=0.01).plot(col='phi', col_wrap=3)
# %%

ds_NE_O2['krm'].sel(seldict).plot(col='phi', col_wrap=3, hue='P')

plt.yscale('log')
#%%
ds_NE_Kp['krm'].sel(seldict).plot(col='phi', col_wrap=3, hue='P')

plt.yscale('log')

# %%
