#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()


data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()


#%%

fig, axes = plt.subplots(3, 1, figsize=(5,10), sharex=True)

var = 'nK_m3_barrel'
line_nK_barrel = axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='Barrel'
    )

var = 'nK_m3_mw_horns'
line_nK_mwhorns = axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='MW Horns'
    )


lineKOH = ds_species_cfd['Yeq_KOH'].pint.to('1/m**3').plot(ax=axes[0], label='CFD: KOH')
linenK = ds_species_cfd['Yeq_K'].pint.to('1/m**3').plot(ax=axes[0], label='CFD: K')


axes[0].set_ylabel("Species Concentration [#/m^3]")
axes[0].legend(
    [line_nK_barrel, line_nK_mwhorns, lineKOH[0], linenK[0]], 
    ['Expt. $n_K$ Barrel', 'Expt $n_K$ Goldi', 'CFD KOH (Goldi)', 'CFD K (Goldi)'],
    bbox_to_anchor=(0.85, 0.9), loc='upper left', framealpha=1
    )

axes[0].set_title('')

var = 'AS_max'
lineAS = axes[1].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='AS Maximum',
    )

axes[1].set_ylabel("AS Maximum")


ta = axes[1].twinx()
#No standard deviation for delta_pd1 as only one run
var = 'delta_pd1'
linePD = ta.plot(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    marker='x',
    label='Delta PD1',
    color='red'
    )

ta.set_ylabel("Delta PD1 [mV]")

axes[1].legend([lineAS, linePD[0]], ['AS Maximum', 'Delta PD1'])

var = 'mws_fit_decay'
axes[2].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit'
    )

for species in ds_tau.data_vars:
    ds_tau[species].plot(label="CFD: {}".format(species), ax=axes[2])

axes[2].legend(bbox_to_anchor=(0.85, 1), loc='upper left', framealpha=1) 
axes[2].set_ylabel("Time Constant [us]")
axes[2].set_title('')


for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')


axes[-1].set_xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params.png'), dpi=300, bbox_inches='tight')
# %%
