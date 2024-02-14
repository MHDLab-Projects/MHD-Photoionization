#%%

from mhdpy.analysis.standard_import import *
from mhdpy.pyvista_utils import calc_rho

import pint_pandas
import pint_xarray
from pint import Quantity

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

ds = xr.load_dataset(fp)

ds['T'] = ds['T'].pint.quantify('K')
ds['p'] = ds['p'].pint.quantify('Pa')

ds['rho'] = calc_rho(ds['T'], ds['p'])  

#%%

ds['rho'].plot(hue='kwt')

#%%

ds['O2'].plot(hue='kwt')

#%%

ds['KOH'].plot(hue='kwt')

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
    sp_rho = sp_rho.pint.to('1/cm^3')
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