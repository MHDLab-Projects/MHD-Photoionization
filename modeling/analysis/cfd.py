#%%

from mhdpy.analysis.standard_import import *
from mhdpy.pyvista_utils import CFDDatasetAccessor

import pint_pandas
import pint_xarray
from pint import Quantity

fp = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf')

ds = xr.load_dataset(fp)

# ds = ds.cfd.convert_all_rho_number()

ds = ds.sel(offset=0)
ds

#%%
ds['O2'].plot(hue='kwt', row='phi')

#%%


ds['KOH'].plot(row='kwt', hue='phi', sharey=False)

#%%

ds['rho'].plot(hue='kwt')

#%%

ds['T'].plot(hue='kwt')

#%%


#%%

ds['Yeq_KOH'].plot(hue='kwt')

plt.yscale('log')

plt.ylim(1e-4,)

#%%

ds[['Yeq_K', 'Yeq_K+']].to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

plt.ylim(1e-6,)


#%%

ds[['K', 'Yeq_K']].to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

#%%

ds[['Kp', 'Yeq_K+']].to_array('var').plot(row='kwt', hue='var')


plt.yscale('log')

#%%



#%%

species = [var for var in ds.data_vars if var not in ['rho', 'T', 'p']]

for species in species:
    sp_rho = ds[species]*ds['rho']
    sp_rho = sp_rho.pint.to('particle/ml')
    ds[species] = sp_rho

#%%

ds_sel = ds[['Yeq_K', 'Yeq_K+']]

g = ds_sel.to_array('var').plot(row='var', hue='kwt')

plt.yscale('log')

goldi_pos = ds['x'].min().item() + 0.18

for ax in g.axes.flatten():
    ax.axvline(goldi_pos, color='gray', linestyle='--')

#%%

ds_sel.sel(x=goldi_pos, method='nearest').to_dataframe()