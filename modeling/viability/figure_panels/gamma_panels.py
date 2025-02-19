# %%
import os
import numpy as np
import xarray as xr
import xyzpy

import matplotlib.pyplot as plt

plt.rcParams.update({
    "savefig.facecolor": 'white',
    # "font.size": 11, 
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
# %%
from matplotlib.colors import LogNorm

# Main text viability figure
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(3.5, 6))

gs = gridspec.GridSpec(2, 1, height_ratios=[0.67, 1])

# Top panel
ax_top = fig.add_subplot(gs[0])

# Bottom panel with four subplots
gs_bottom = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1], wspace=0.1, hspace=0.1)

ax_bottom_00 = fig.add_subplot(gs_bottom[0, 0])
ax_bottom_10 = fig.add_subplot(gs_bottom[1, 0])
ax_bottom_11 = fig.add_subplot(gs_bottom[1, 1])
ax_bottom_01 = fig.add_subplot(gs_bottom[0, 1])

ax_bottom_01.sharey(ax_bottom_00)
ax_bottom_11.sharey(ax_bottom_10)

ax_bottom_00.sharex(ax_bottom_10)
ax_bottom_01.sharex(ax_bottom_11)

axes_bottom = np.array([[ax_bottom_00, ax_bottom_01], [ax_bottom_10, ax_bottom_11]])

# Top panel of viability figure
combo_sel = dict(l_b=0, Kwt=0.01, phi=0.8, eta='perf', rxn='mm_sum')

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)

norm = LogNorm(vmin=0.01, vmax=100)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm, ax=ax_top)

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--', ax=ax_top)

# Set the colorbar label
g.colorbar.set_label('$\\gamma$')

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--', ax=ax_top)

ax_top.set_title('')

ax_top.set_ylabel('Temperature [K]')
ax_top.set_xlabel('Pressure [Pa]')

ax_top.set_xlim(0.6e4, 1.4e6)
ax_top.set_ylim(1200, 3500)

# Bottom panel of main text viability figure
combo_downsel = {
    'l_b': [0, 0.99],
    'eta': ['perf', 'KOH'],
    'Kwt': [0.01],
    'rxn': 'mm_sum'
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

for i, eta in enumerate(['perf', 'KOH']):
    for j, l_b in enumerate([0, 0.99]):
        ax = axes_bottom[j, i]

        da_sel = P_zero.sel(eta=eta, l_b=l_b)

        for phi in da_sel['phi'].values:
            da_sel.sel(phi=phi).plot(y='T', xscale='log', ax=ax, label=f'$\phi$ = {phi}')
        ax.set_title(f'$\eta$ = {eta}, $l_b$ = {l_b}')

        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')

labs = [0.8, 0.9, 1.0]

# Create a single legend for the right four axes and modify its position
handles, labels = axes_bottom[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, title='Equivalence\nRatio', loc='center right', bbox_to_anchor=(1.22, 0.56))


for ax in axes_bottom.flatten():
    ax.set_xlim(0.6e4, 1.4e6)
    ax.set_ylim(1200, 3500)

axes_bottom[0, 0].set_title(r'$\eta$ = 1.0')
axes_bottom[0, 1].set_title(r'$\eta_{KOH}$')

axes_bottom[1, 0].set_xlabel('Pressure [Pa]')
axes_bottom[1, 1].set_xlabel('Pressure [Pa]')

axes_bottom[0, 0].set_ylabel('T [K]')
axes_bottom[1, 0].set_ylabel('T [K]')

axes_bottom[0, 1].text(1.05, 0.5, '$l_{bl} = 1$', transform=axes_bottom[0, 1].transAxes, rotation=-90, va='center')
axes_bottom[1, 1].text(1.05, 0.5, '$l_{bl} = 0.01$', transform=axes_bottom[1, 1].transAxes, rotation=-90, va='center')

# Remove y tick labels for shared y axes
plt.setp(ax_bottom_01.get_yticklabels(), visible=False)
plt.setp(ax_bottom_11.get_yticklabels(), visible=False)

# Remove x tick labels for shared x axes
plt.setp(ax_bottom_00.get_xticklabels(), visible=False)
plt.setp(ax_bottom_01.get_xticklabels(), visible=False)

# fig.text(0.005, 0.98, 'A)', transform=fig.transFigure, verticalalignment='top')
# fig.text(0.005, 0.48, 'B)', transform=fig.transFigure, verticalalignment='top')

fig.tight_layout(h_pad=2)


# Add labels to the subplots
axes = [ax_top, ax_bottom_00]
labels_ab = ['A)', 'B)']

for ax, label in zip(axes, labels_ab):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - 0.2, Y , label)
 

fig.savefig(os.path.join(REPO_DIR, 'final','figures', 'Fig8_viability.svg'))

# %%
