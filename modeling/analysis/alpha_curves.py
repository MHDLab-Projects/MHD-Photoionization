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

cantera_data_dir = os.path.join(REPO_DIR, 'modeling','dataset','output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_P_zero = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'P_zero.cdf'))

# %%


# alpha = xr.open_dataset(os.path.join(PI_modeling_dataset_dir,'alpha_bl.cdf')).squeeze()['alpha_bl']
ds_NE = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE.cdf')).squeeze()
alpha = ds_NE['alpha']


# Add enhancement factor
da_dsigma_tot = xr.load_dataset(os.path.join(PI_modeling_dataset_dir,'da_dsigma_tot.cdf'))['enhancement factor']
alpha = alpha*da_dsigma_tot


beta = alpha -1 



# %%


combo_sel = dict(l_bk=0, P_in=0, Kwt=0.01, phi=0.8, analysis='perf_Bconst')

cmap = plt.get_cmap('RdBu')
beta_sel = beta.sel(combo_sel)
g = beta_sel.plot(vmin=-1.2,vmax=1.2,xscale='log', cmap=cmap)

# Set the colorbar label
g.colorbar.set_label('$\\alpha - 1$')


ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=4, linestyle='--')

plt.gca().set_title('')

plt.ylabel('Temperature (K)')
plt.xlabel('Pressure (Pa)')



plt.savefig('output/alpha_curve_demo.png')

# %%


combo_downsel = {
    'P_in' : 0,
    'l_bk': 0,
    'Kwt': [0.001, 0.01, 0.1],
    'phi': [0.8,1,1.2],
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero.plot(col='phi', row='Kwt',hue='analysis', y='T', xscale='log')

# %%


combo_downsel = {
    'l_bk': [0],
    'phi': [0.8],
    'Kwt': [0.01],
    'analysis': ['perf_Bconst'],
    'P_in': [0,1e-6,1e-2,1e2,1e6]
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero.plot(hue='P_in', y='T', xscale='log')


# %%

combo_downsel = {
    'P_in' : 0,
    'phi': [0.8],
    'Kwt': [0.001,0.01,0.1]
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero.plot(col='l_bk', row='Kwt', hue='analysis', y='T', xscale='log')

# %%

combo_downsel = {
    'P_in' : 0,
    'l_bk': 0,
    'phi': [0.8],
    'Kwt': [0.01]
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero.plot(hue='analysis', y='T', xscale='log')

#%%

combo_downsel = {
    'P_in' : 0,
    'phi': [0.8],
    'Kwt': [0.01],
    # 'analysis': ['perf_Bconst']
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

g = P_zero.plot(hue='l_bk', col='analysis', y='T', xscale='log', col_wrap=2)

for ax in g.axes.flatten():
    ax.plot([1e5], [3000], marker='*', markersize=10)

plt.savefig('output/alpha_curve_analysis_lbk.png')