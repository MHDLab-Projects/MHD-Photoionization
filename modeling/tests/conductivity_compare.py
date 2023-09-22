"""
Check that the calculated conductivity at P_in=0 gives the same result as cantera conductivity
This is not really a test at this point because the values are slightly different, but this script prints them to confirm they are nearly identical

"""

# %%

import mhdpy

defaults = {
    'phi': 0.7,
    'Kwt': 0.001
}

from mhdpy.load import load_dataset_defaults

canterapath = mhdpy.fp.gen_path('sharepoint', 'lee', 'Data', 'Cantera Gasses')

ds_TP_params = load_dataset_defaults(os.path.join(canterapath, 'ds_TP_params.cdf'), defaults)
ds_TP_species_rho = load_dataset_defaults(os.path.join(canterapath, 'ds_TP_species_rho.cdf'), defaults)

sig_cant = ds_TP_params['sigma']/100 # convert to S/cm

# %%

import sys
sys.path.append('..')

from noneq_utils import noneq

G_NE = 0
G_th = ds_TP_params['Gth']
kr = ds_TP_params['kr']
mue_cant = ds_TP_params['mobility']*10000, #convert to cm2/Vs

ne = noneq.calc_ne(G_th, G_NE, kr)

sig_NE = noneq.calc_sig(ne, mue_cant)
sig_NE.name = 'sigma'

# %%

print(sig_cant.sel(P=1e5, T=3000, method='nearest'))
print(sig_NE.sel(P=1e5, T=3000, method='nearest'))
# %%
