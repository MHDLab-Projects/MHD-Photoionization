# %%
import os
import numpy as np
import xarray as xr
import xyzpy

import matplotlib.pyplot as plt
import pint_xarray

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

cantera_data_dir = os.path.join(REPO_DIR, 'modeling','dataset','output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf')).sel({'phi': 0.8, 'Kwt': 0.01})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf')).sel({'phi': 0.8, 'Kwt': 0.01})
#%$ 

ds_TP_params['krb'].pint.quantify()


# %%
ds_TP_params['krb'].plot()

#%%

da_sel = ds_TP_params.sel(T=slice(1500,2000)).sel(P=1e5)['krb']

da_sel = da_sel.pint.quantify().pint.to('cm^3/(particle s)')


da_sel.plot()


# %%
