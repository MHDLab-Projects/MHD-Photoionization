
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

from pi_paper_utils import abscs, noneq

cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf')).sel({'phi': 0.7})

seldict = {
    'P_combustor': ds_TP_params.coords['P_combustor'].item(),
    'inlet_T': ds_TP_params.coords['inlet_T'].item()
}

ds_HP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_HP_species.cdf')).sel(seldict)
ds_HP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_HP_params.cdf')).sel(seldict)
# %%


ds_cs = abscs.calc_ds_cs(ds_TP_species_rho.coords['T'].values)
ds_cs.coords['wl'] = ds_cs.coords['wl'].round(2)

# %%

fig, axes = plt.subplots(1,2, sharey = True, figsize = (7,4))
for var in ds_cs:
    ds_cs[var].sel(T = 2000, method='nearest').plot(ax = axes[0], yscale = 'log')
    ds_cs[var].sel(wl = 248, method = 'nearest').plot(ax = axes[1], label = var, yscale = 'log')

plt.legend(bbox_to_anchor = [0,0,1.4,1])

plt.savefig('output/abscs.png', bbox_inches = 'tight')

# %%

if len(ds_cs.coords['wl']) > 10:
    ds_cs = ds_cs.sel(wl=ds_cs.coords['wl'].values[::20])
    
gas_lam = noneq.calc_atten_lengths(ds_cs, ds_TP_species_rho)


das = []
for species in gas_lam:
    if species != 'tot':
        F = gas_lam['tot']/gas_lam[species]
        F.name = 'F_' + species
        F.attrs = dict(long_name = '$f_' + species + '$')
        das.append(F)
#         gas_lam = gas_lam.assign(temp=F).rename(temp=species)

f_species = xr.merge(das)

# %%

fig, axes = plt.subplots(3,1, figsize = (4,6), sharex = True)

ds_sel = gas_lam.sel(P = 1e5,  Kwt = 0.01, method = 'nearest')

lns = ds_sel['tot'].plot(hue = 'wl', yscale = 'log', ax  = axes[0])
# axes[0].get_legend().set_bbox_to_anchor((1.1, 1.05))
axes[0].get_legend().remove()
# axes[0].set_prop_cycle(None)
lns = ds_sel['KOH'].plot(hue = 'wl', yscale = 'log', ax  = axes[1])
# axes[1].get_legend().remove()

axes[1].get_legend().set_bbox_to_anchor((1.1, 1.05))

lns = f_species['F_KOH'].sel(P = 1e5,  Kwt = 0.01, method = 'nearest').plot(hue = 'wl', yscale = 'log', ax=axes[2])
axes[2].get_legend().remove()
axes[2].set_ylabel('$f_{KOH}$')
# lns[0].axes.get_legend().set_bbox_to_anchor((1.1, 1.05))
# fig.tight_layout()

for ax in axes:
    ax.set_title('')
plt.savefig('output/atten_length.png')
# %%
