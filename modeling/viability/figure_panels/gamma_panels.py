
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

#%%
ds_P_zero['P_zero'].sel(eta='perf', l_b=0.99, rxn='mm_sum').plot(col='phi',hue='Kwt', y='T')

plt.xscale('log')

plt.savefig(os.path.join('output', 'viability_kwt.png'))
# %%

from matplotlib.colors import LogNorm

plt.figure(figsize=(5,3))

combo_sel = dict(l_b=0,  Kwt=0.01, phi=0.8, eta='perf', rxn='mm_sum')

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)

norm = LogNorm(vmin=0.01, vmax=100)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm)

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')


# Set the colorbar label
g.colorbar.set_label('$\\gamma$')

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')

plt.gca().set_title('')

plt.ylabel('Temperature [K]')
plt.xlabel('Pressure [Pa]')

plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)


plt.savefig('output/gamma_curve_demo.png')
# %%


combo_downsel = {
    # 'P_in' : 0,
    # 'phi': [0.8],
    'l_b': [0, 0.99],  
    'eta': ['perf', 'KOH'],
    'Kwt': [0.01],
    'rxn': 'mm_sum'
    # 'analysis': ['perf_Bconst']
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)


fig, axes = plt.subplots(2, 2, figsize=(4,3), sharex=True, sharey=True)


for i, eta in enumerate(['perf', 'KOH']):
    for j, l_b in enumerate([0, 0.99]):
        ax = axes[j,i]

        da_sel = P_zero.sel(eta=eta, l_b=l_b)

        lns = da_sel.plot(hue='phi', y='T', xscale='log', ax=ax)
        ax.set_title(f'$\eta$ = {eta}, $l_b$ = {l_b}')

        ax.get_legend().remove()
        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')


labs = [0.8, 0.9, 1.0]

leg = fig.legend(lns, labs, title='Equiv.\nRatio')
leg.set_bbox_to_anchor((1.19, 0.7))  # coordinates are in figure units

plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)


axes[0,0].set_title(r'$\eta$ = 1.0')
axes[0,1].set_title(r'$\eta_{KOH}$')

axes[1,0].set_xlabel('Pressure [Pa]')
axes[1,1].set_xlabel('Pressure [Pa]')

axes[0,0].set_ylabel('T [K]')
axes[1,0].set_ylabel('T [K]')

axes[0,1].text(1.05, 0.5, '$l_{bl} = 1$', transform=axes[0,1].transAxes, rotation=-90, va='center')
axes[1,1].text(1.05, 0.5, '$l_{bl} = 0.01$', transform=axes[1,1].transAxes, rotation=-90, va='center')

plt.tight_layout(pad=0.6)

plt.savefig('output/P_zero_l_b_eta.png')
# %%

combo_downsel = {
    # 'P_in' : 0,
    'phi': [0.8, 0.9, 1.0],
    'l_b': [0],  
    'eta': ['perf', 'KOH'],
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