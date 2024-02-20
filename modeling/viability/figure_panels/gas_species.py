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

cantera_data_dir = os.path.join(REPO_DIR, 'modeling' , 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))

seldict = {
    'P_combustor': ds_TP_params.coords['P_combustor'].item(),
    'inlet_T': ds_TP_params.coords['inlet_T'].item()
}

ds_HP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_HP_species.cdf')).sel(seldict)
ds_HP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_HP_params.cdf')).sel(seldict)

# %%

plt.figure(figsize=(6,3))

ds_sel = ds_TP_species.sel(P = 1e5, Kwt=0.01, phi = 0.8, method = 'nearest')

ds_sel = ds_sel.drop(['KO','KOH'])

for var in ds_sel.data_vars:
    if not np.isnan(ds_sel[var].where(ds_sel[var] > 1e-4)).all():
        if 'K' in var or 'e' in var:
            linestyle = '-'
            linewidth=3
        else:
            linestyle = '--'
            linewidth=2
        if 'e' in var:
            color = 'black'
            
        else:
            color = None
        
        ds_sel[var].plot(label = var, linestyle=linestyle, linewidth=linewidth, color=color)
lg = plt.legend()
lg.set_bbox_to_anchor([0,0,1.3,1])
plt.ylabel('Species Concentration')
plt.yscale('log')
# plt.xscale('log')
plt.ylim(1e-6,1)
plt.gca().set_title('')

plt.tight_layout()

plt.savefig('output/species_concentrations_T.png')

# %%
# K_species = [var for var in ds_TP_species_Kwt.data_vars if 'K' in var]
K_species = [var for var in ds_TP_species.data_vars if 'K' in var]

ds_TP_species_Kwt = ds_TP_species.copy()
ds_TP_species_Kwt.coords['Kwt'] = ds_TP_species_Kwt.coords['Kwt']*100
ds_TP_species_Kwt.coords['Kwt'].attrs = dict(long_name='K wt%')

total_K_frac = ds_TP_species_Kwt[K_species].to_array('species').sum('species')

ion_frac = ds_TP_species_Kwt['K+']/total_K_frac
# ion_frac.coords['Kwt']

# %%

plt.figure(figsize=(3,3))

ion_frac.sel(P=1e5, phi=0.8).sel(Kwt=[0.1,1,10]).plot(hue='Kwt', yscale='log', ylim=(1e-4, 1))

plt.grid()
plt.gca().set_title('')

plt.ylabel('Ionization Fraction')

plt.tight_layout()

plt.savefig('output/ionization_fraction_T.png')

# %%
ds_TP_params_Kwt = ds_TP_params.copy()
ds_TP_params_Kwt.coords['Kwt'] = ds_TP_params_Kwt.coords['Kwt']*100
ds_TP_params_Kwt.coords['Kwt'].attrs = dict(long_name='K wt%')

# %%

plt.figure(figsize=(3,3))

ds_TP_params_Kwt['sigma'].sel(P=1e5, phi=0.8).sel(Kwt=[0.1,1,10]).plot(hue='Kwt', yscale='log', ylim=(1e-4, 1e4))

plt.grid()
plt.gca().set_title('')
plt.gca().get_legend().remove()

plt.tight_layout()

plt.savefig('output/elecrical_conductivity_T.png')
# %%


# %%



