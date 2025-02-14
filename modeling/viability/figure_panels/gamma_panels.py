
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

ds_NE = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE.cdf'), engine = 'scipy').squeeze()
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

fig = plt.figure(figsize=(8, 4))

gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

# Left panel
ax_left = fig.add_subplot(gs[0])

# Right panel with four subplots
gs_right = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[1])

ax_right_00 = fig.add_subplot(gs_right[0, 0])
ax_right_01 = fig.add_subplot(gs_right[0, 1], sharex=ax_right_00, sharey=ax_right_00)
ax_right_10 = fig.add_subplot(gs_right[1, 0], sharey=ax_right_00)
ax_right_11 = fig.add_subplot(gs_right[1, 1], sharex=ax_right_10, sharey=ax_right_00)

axes_right = np.array([[ax_right_00, ax_right_01], [ax_right_10, ax_right_11]])

# Left panel of viability figure
combo_sel = dict(l_b=0, Kwt=0.01, phi=0.8, eta='perf', rxn='mm_sum')

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)

norm = LogNorm(vmin=0.01, vmax=100)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm, ax=ax_left)

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--', ax=ax_left)

# Set the colorbar label
g.colorbar.set_label('$\\gamma$')

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--', ax=ax_left)

ax_left.set_title('')

ax_left.set_ylabel('Temperature [K]')
ax_left.set_xlabel('Pressure [Pa]')

ax_left.set_xlim(0.8e4, 1.2e6)
ax_left.set_ylim(1200, 3500)

# Right panel of main text viability figure
combo_downsel = {
    'l_b': [0, 0.99],
    'eta': ['perf', 'KOH'],
    'Kwt': [0.01],
    'rxn': 'mm_sum'
}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

for i, eta in enumerate(['perf', 'KOH']):
    for j, l_b in enumerate([0, 0.99]):
        ax = axes_right[j, i]

        da_sel = P_zero.sel(eta=eta, l_b=l_b)

        for phi in da_sel['phi'].values:
            da_sel.sel(phi=phi).plot(y='T', xscale='log', ax=ax, label=f'$\phi$ = {phi}')
        ax.set_title(f'$\eta$ = {eta}, $l_b$ = {l_b}')

        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('')

labs = [0.8, 0.9, 1.0]
# Create a single legend for the right four axes
handles, labels = axes_right[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, title='Equiv.\nRatio', loc='center right', bbox_to_anchor=(1.13, 0.5))


plt.xlim(0.8e4, 1.2e6)
plt.ylim(1200, 3500)

axes_right[0, 0].set_title(r'$\eta$ = 1.0')
axes_right[0, 1].set_title(r'$\eta_{KOH}$')

axes_right[1, 0].set_xlabel('Pressure [Pa]')
axes_right[1, 1].set_xlabel('Pressure [Pa]')

axes_right[0, 0].set_ylabel('T [K]')
axes_right[1, 0].set_ylabel('T [K]')

axes_right[0, 1].text(1.05, 0.5, '$l_{bl} = 1$', transform=axes_right[0, 1].transAxes, rotation=-90, va='center')
axes_right[1, 1].text(1.05, 0.5, '$l_{bl} = 0.01$', transform=axes_right[1, 1].transAxes, rotation=-90, va='center')

# Remove y tick labels for shared y axes
plt.setp(ax_right_01.get_yticklabels(), visible=False)
plt.setp(ax_right_11.get_yticklabels(), visible=False)

fig.tight_layout(w_pad = 2)
# Add labels to the subplots
fig.text(0.005, 0.98, 'A)', transform=fig.transFigure, verticalalignment='top')
fig.text(0.505, 0.98, 'B)', transform=fig.transFigure, verticalalignment='top', horizontalalignment='left')


plt.savefig('output/P_zero_l_b_eta.png')
