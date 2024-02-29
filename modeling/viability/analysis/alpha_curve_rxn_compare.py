# %%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils


cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_P_zero = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'P_zero.cdf'))

ds_NE = xr.open_dataset(pjoin(REPO_DIR, 'modeling', 'viability', 'dataset', 'output', 'ds_NE_rxn_comp.cdf')).squeeze()

#%%

alpha = ds_NE['alpha']

#%%

# alpha = ds_alpha.to_array('rxn').rename('alpha')
# alpha = ds_NE['alpha']

# Add enhancement factor
da_dsigma_tot = xr.open_dataset(pjoin(REPO_DIR, 'modeling', 'viability', 'dataset', 'output', 'da_dsigma_tot.cdf'))['enhancement factor']
alpha = alpha*da_dsigma_tot

# alpha.sel(Kwt=0.1, phi=1, analysis='PI_Bhall', P_in=0, T=1.6e3).plot()

#TODO: gamma used to be beta, search for other instances. 
gamma = alpha -1 

# %%

# %%
combo_downsel = {
    'P_in' : 0,
    'l_bk': 0,
    'Kwt': 0.01,
    'phi': [0.8,1,1.2],
    # 'analysis': 'perf_Bconst'

}

P_zero = ds_P_zero['P_zero'].sel(combo_downsel)

P_zero
# %%

P_zero.plot(col='rxn', hue='phi', y='T')

plt.xscale('log')

#%%

g = P_zero.plot(col='phi', hue='rxn', y='T', figsize=(4,2))


plt.xscale('log')
g.axes[0][0].set_ylabel('Temperature (K)')

for ax in g.axes.flatten():
    ax.set_xlabel('Pressure (Pa)')

plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)

plt.savefig(pjoin(DIR_FIG_OUT, 'P_zero_rxn_phi.png'))


# %%
combo_sel = dict(l_bk=0, P_in=0, Kwt=0.01, phi=0.8)

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)

#TODO: fig size with number of ticks
g = gamma_sel.plot(vmin=-1.2,vmax=1.2,xscale='log', col='rxn', cmap=cmap, figsize=(5,2.6))

# Set the colorbar label
# g.fig.colorbar().set_label('$\\alpha - 1$')


ax_O2 = g.axes[0][1]
line_O2 = ds_P_zero['P_zero'].sel(combo_sel).sel(rxn='O2')
line_O2.plot(y='T', color='green', linewidth=2, linestyle='--', ax=ax_O2)

ax_O2.set_title('O2')

ax_Kp = g.axes[0][0]
line_Kp = ds_P_zero['P_zero'].sel(combo_sel).sel(rxn='Kp')
line_Kp.plot(y='T', color='green', linewidth=2, linestyle='--', ax=ax_Kp)

ax_Kp.set_title('K+')

g.axes[0][0].set_ylabel('Temperature (K)')
g.axes[0][1].set_ylabel('')


x_ticks = np.logspace(3, 7, num=5)  # 5 ticks between 10^3 and 10^7
for ax in g.axes.flatten():
    ax.set_xlabel('Pressure (Pa)')
    ax.set_xticks(x_ticks)
    ax.set_xlim(0.8e4,1.2e6)
    ax.set_ylim(1200,3500)


# Get the colorbar
colorbar = g.axes[0, -1].collections[0].colorbar

# Set the colorbar label

colorbar.set_label('$1-\\beta$')

plt.savefig(pjoin(DIR_FIG_OUT, 'alpha_curve_demo_rxn.png'))
# %%
