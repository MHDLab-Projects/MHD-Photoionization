
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

#%%

ds_P_zero = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'P_zero.cdf'))

ds_NE = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE.cdf')).squeeze()
gamma = ds_NE['gamma']

# Add enhancement factor
da_dsigma_tot = xr.load_dataset(os.path.join(PI_modeling_dataset_dir,'da_dsigma_tot.cdf'))['enhancement factor']
gamma = gamma*da_dsigma_tot

beta = gamma -1 

# %%

from matplotlib.colors import LogNorm

plt.figure(figsize=(5,3))

combo_sel = dict(l_bk=0,  Kwt=0.01, phi=0.8, eta='perf', rxn='mm_sum')

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)

norm = LogNorm(vmin=0.01, vmax=100)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm)

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')


# Set the colorbar label
g.colorbar.set_label('$\\gamma$')

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')

plt.gca().set_title('')

plt.ylabel('Temperature (K)')
plt.xlabel('Pressure (Pa)')

plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)


plt.savefig('output/gamma_curve_demo.png')
# %%


combo_downsel = {
    # 'P_in' : 0,
    # 'phi': [0.8],
    'l_bk': [0, 0.99],  
    'eta': ['perf', 'PI'],
    'Kwt': [0.01],
    'rxn': 'mm_sum'
    # 'analysis': ['perf_Bconst']
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

g = P_zero.plot(hue='phi', col='eta', row='l_bk', y='T', xscale='log', figsize=(5,3))

# for ax in g.axes.flatten():
#     ax.plot([1e5], [3000], marker='*', markersize=10)
# Get the legend and move it
legend = g.fig.legends[0]
legend.set_bbox_to_anchor((0.3, 0.35))  # coordinates are in figure units


plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)


g.axes[1,0].set_xlabel('Pressure (Pa)')
g.axes[1,1].set_xlabel('Pressure (Pa)')

plt.savefig('output/P_zero_l_bk_eta.png')
# %%

combo_downsel = {
    # 'P_in' : 0,
    'phi': [0.8, 1, 1.2],
    'l_bk': [0],  
    'eta': ['perf', 'PI'],
    'Kwt': [0.01],
    # 'rxn': 'mm_sum'
    # 'analysis': ['perf_Bconst']
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

g = P_zero.plot(hue='rxn', row='eta', col='phi', y='T', xscale='log', figsize=(5,3))

# for ax in g.axes.flatten():
#     ax.plot([1e5], [3000], marker='*', markersize=10)
# Get the legend and move it
legend = g.fig.legends[0]
legend.set_bbox_to_anchor((1.2, 0.5))  # coordinates are in figure units


plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)

plt.tight_layout()


plt.savefig('output/P_zero_rxn_component.png')