#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils #Sets matplotlib style


data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()
ds_mws_fit = xr.open_dataset(pjoin(data_directory, '53x_ds_fit_mws.cdf')).pint.quantify()

#%%

plt.figure(figsize=(4,2))

ds_plot = ds_mws_fit.sel(date='2023-05-12').sel(run_num=1)

# da_plot = ds_plot[['AS_fit','AS_all']].to_array('var')


ds_plot['AS_all'].plot(hue='kwt')


plt.gca().get_legend().set_title('K wt %')

for kwt in ds_plot.coords['kwt']:
    ds_plot['AS_fit'].sel(kwt=kwt).plot(color='black', linestyle='-', alpha=0.5)


plt.xlim(-1,40)

plt.yscale('log')
plt.ylim(1e-5, 0.1)

plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_exp.png'), dpi=300, bbox_inches='tight')


#%%

fig, axes = plt.subplots(2, 1, figsize=(5,8), sharex=True)

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


lineKOH = ds_species_cfd['Yeq_KOH'].pint.to('particle/m**3').plot(ax=axes[0], label='CFD: KOH')
linenK = ds_species_cfd['Yeq_K'].pint.to('particle/m**3').plot(ax=axes[0], label='CFD: K')
line_allK = ds_species_cfd['all_K_Yeq'].pint.to('particle/m**3').plot(ax=axes[0], label='CFD: All K')


axes[0].set_ylabel("Species Concentration [#/m^3]")
axes[0].legend(
    [line_nK_barrel, line_nK_mwhorns, lineKOH[0], linenK[0], line_allK[0]], 
    ['Expt. $n_K$ Barrel', 'Expt $n_K$ Goldi', 'CFD KOH (Goldi)', 'CFD K (Goldi)' , 'CFD All K (Goldi)'],
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

for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')

axes[-1].set_xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_ionization.png'), dpi=300, bbox_inches='tight')
# %%


fig, axes = plt.subplots(2, 1, figsize=(5,5), sharex=True)

ax = axes[0]

var = 'mws_fit_decay'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit'
    )


var = 'mws_fit_decay_exp'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit Expon.'
    )

for species in ds_tau.data_vars:
    ds_tau[species].plot(label="CFD: {}".format(species), ax=ax)

ax.legend(bbox_to_anchor=(1, 1), loc='upper left', framealpha=1) 
ax.set_ylabel("Time Constant [us]")
ax.set_title('')
ax.set_ylim(1e-1, 1e5)

ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_ylim(1e-2, 1e6)


ax = axes[1]

var = 'mws_fit_dne'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit'
    )

plt.yscale('log')

plt.ylabel("$\Delta n_e$ [um$^{-3}$]")
plt.xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_recomb.png'), dpi=300, bbox_inches='tight')


# %%
