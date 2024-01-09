#%%

from mhdpy.analysis.standard_import import *

import pint_pandas
import pint_xarray
from pint import Quantity

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles.cdf')

ds = xr.load_dataset(fp)

#%%

ds['T'] = ds['T'].pint.quantify('K')
ds['p'] = ds['p'].pint.quantify('Pa')

#%%

R = Quantity(8.31446261815324, 'J/(mol*K)')
ds['rho'] = ds['p']/(ds['T']*R)

NA = Quantity(6.02214076e23, '1/mol')

ds['rho'] = ds['rho']*NA

ds['rho'].pint.to('1/cm^3')
#%%

ds['rho'].plot(hue='kwt')


#%%

species = [var for var in ds.data_vars if var not in ['rho', 'T', 'p']]

for species in species:
    sp_rho = ds[species]*ds['rho']
    sp_rho = sp_rho.pint.to('1/cm^3')
    ds[species] = sp_rho

#%%

ds_sel = ds[['K', 'Kp']]

g = ds_sel.to_array('var').plot(row='var', hue='kwt')

plt.yscale('log')

goldi_pos = ds['x'].min().item() + 0.18

for ax in g.axes.flatten():
    ax.axvline(goldi_pos, color='gray', linestyle='--')

#%%

ds_sel.sel(x=goldi_pos, method='nearest').to_dataframe()