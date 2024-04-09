#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu #Sets matplotlib style

import matplotlib.pyplot as plt
import matplotlib.cm as cm

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()

ds_mws_fit_exp = xr.open_dataset(pjoin(data_directory, '53x_ds_fit_mws_exp.cdf')).pint.quantify()
ds_mws_fit_dnedt = xr.open_dataset(pjoin(data_directory, '53x_ds_fit_mws_dnedt.cdf')).pint.quantify()


# Load in the raw lecroy data for confidence intervals. Doing this here to avoid having to save raw mnum acquisitions...Revisit. 
from mhdpy.fileio.ct import load_df_cuttimes
from mhdpy.coords.ct import downselect_acq_time
ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute')
# Downselect times to those in the CFD input kwt cases 
fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)


#TODO: the time grids are not the same for the two fitting methods. Why?
ds_mws_fit_exp = ds_mws_fit_exp.interp(time=ds_mws_fit_dnedt.time, method='linear')
#TODO: realized 'method' cannot be used as a coordinate as it is a reserved word in xarray(i.e. method='nearest'). Need to change this in the fitting code.
ds_mws_fit = xr.concat([
    ds_mws_fit_exp.assign_coords(fit_method=['exp']),
    ds_mws_fit_dnedt.assign_coords(fit_method=['dnedt'])
    ], 'fit_method')

# ds_mws_fit = ds_mws_fit.assign_coords(method=[str(s) for s in ds_mws_fit.coords['method'].values])

#%%

# ds_mws_fit_dnedt.coords['time']
ds_mws_fit_exp.coords['time']



#%%

plt.figure(figsize=(4,2))

ds_plot = ds_mws_fit.sel(date='2023-05-12').sel(run_num=1).sel(fit_method='exp')

# da_plot = ds_plot[['AS_fit','AS_all']].to_array('var')


ds_plot['AS_all'].plot(hue='kwt')


plt.gca().get_legend().set_title('K wt %')

for kwt in ds_plot.coords['kwt']:
    ds_plot['AS_fit'].sel(kwt=kwt).plot(color='black', linestyle='-', alpha=0.5)


plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5, 0.1)

plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_exp.png'), dpi=300, bbox_inches='tight')


#%%

ds_plot = ds_mws_fit.sel(kwt=1, method='nearest').sel(date='2023-05-12').sel(run_num=1)

g = ds_plot[['AS_fit','AS_all']].to_array('var').plot(hue='var', col='fit_method', sharey=False)

g.axes[0][1].set_ylabel('$\Delta AS$ Normalized [dimensionless]')

for ax in g.axes.flatten():
    ax.set_yscale('log')

plt.xlim(-1,40)
g.axes[0][0].set_ylim(1e-5, 2e-1)
g.axes[0][1].set_ylim(1e-4, 2)

g.axes[0][0].set_title('Method: Exponential')
g.axes[0][1].set_title('Method: Diff. Eq')


plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_exp_dnedt_individual_compare.png'), dpi=300, bbox_inches='tight')

#%%

plt.figure()

da_plot = ds_plot[['AS_fit','AS_all']].to_array('var').sel(fit_method='exp')

g = da_plot.plot(hue='var', figsize=(5,2))

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5, 0.1)

plt.title('')

# plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_exp.png'), dpi=300, bbox_inches='tight')
#%%
# %%


fig, axes = plt.subplots(2, 1, figsize=(5,5), sharex=True)

ax = axes[0]

var = 'mws_fit_decay'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit Diff. Eq.'
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
    label='MWS Fit Diff. Eq.'
    )

plt.yscale('log')

plt.ylabel("$\Delta n_e$ [um$^{-3}$]")
plt.xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_recomb_both.png'), dpi=300, bbox_inches='tight')


# %%

ds_2 = ds_lecroy.mws.calc_time_stats()

ds_2

#%%

ds_2['dpd1'].mean('run').mean('mnum').plot(hue='kwt')

plt.xlim(-2,3)


#%%


# ds_p_stats
# ds_species_cfd

fig, axes = plt.subplots(1, 1, figsize=(5,5), sharex=True)

ax = axes


var = 'dAS_abs_max'
ax.errorbar(
    ds_species_cfd['Yeq_K'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='MWS Fit Diff. Eq.'
    )

plt.yscale('log')

plt.ylabel("$\Delta AS$ Maximum")

plt.xlabel("K [um$^{-3}$]")
# plt.xscale('log')

#%%
# With a colorbar corresponding to K value

# Create a colormap
cmap = cm.get_cmap('viridis')

# Normalize the 'kwt' values to the range [0, 1] for the colormap
norm = plt.Normalize(ds_species_cfd['kwt'].min(), ds_species_cfd['kwt'].max())

fig, axes = plt.subplots(2, 1, figsize=(5,5), sharex=True)

vars = ['dAS_abs_max', 'dpd1_max']

xvar = 'Yeq_K'

for i, var in enumerate(vars):
    ax = axes[i]

    # Add connecting lines
    
    ax.plot(
        ds_species_cfd[xvar], 
        ds_p_stats['{}_mean'.format(var)], 
        color='black'
    )

    # No error bars for dpd1_max as only one run

    if var != 'dpd1_max':
        # Add error bars
        ax.errorbar(
            ds_species_cfd[xvar], 
            ds_p_stats['{}_mean'.format(var)], 
            yerr=ds_p_stats['{}_std'.format(var)], 
            fmt='', capsize=5, color='black', alpha=1.0,
            zorder=1
        )


    # Use scatter instead of errorbar to be able to specify a color for each point
    sc = ax.scatter(
        ds_species_cfd[xvar], 
        ds_p_stats['{}_mean'.format(var)], 
        c=ds_species_cfd['kwt'], cmap=cmap, norm=norm,
        marker='o',
        zorder=2
    )

    # Add a colorbar
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Nominal K wt %")


    ax.set_yscale('log')


    ax.set_xlabel("CFD K [cm$^{-3}$]")

axes[0].set_ylabel("$\Delta AS$ Maximum")
axes[1].set_ylabel("Delta PD1 [mV]")
axes[0].set_ylim(5e-3,1.1e-1)
axes[1].set_ylim(3.5e-1,2.5e0)

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_ionization_vsCFDK.png'), dpi=300, bbox_inches='tight')

# %%

# Create a colormap
cmap = cm.get_cmap('viridis')

# Normalize the 'kwt' values to the range [0, 1] for the colormap
norm = plt.Normalize(ds_species_cfd['kwt'].min(), ds_species_cfd['kwt'].max())

fig, axes = plt.subplots(2, 1, figsize=(5,5), sharex=True)

vars = ['dAS_abs_max', 'dpd1_max']

xvar = 'Yeq_KOH'

for i, var in enumerate(vars):
    ax = axes[i]

    # Add connecting lines
    
    ax.plot(
        ds_species_cfd[xvar], 
        ds_p_stats['{}_mean'.format(var)], 
        color='black'
    )

    # No error bars for dpd1_max as only one run

    if var != 'dpd1_max':
        # Add error bars
        ax.errorbar(
            ds_species_cfd[xvar], 
            ds_p_stats['{}_mean'.format(var)], 
            yerr=ds_p_stats['{}_std'.format(var)], 
            fmt='', capsize=5, color='black', alpha=1.0,
            zorder=1
        )


    # Use scatter instead of errorbar to be able to specify a color for each point
    sc = ax.scatter(
        ds_species_cfd[xvar], 
        ds_p_stats['{}_mean'.format(var)], 
        c=ds_species_cfd['kwt'], cmap=cmap, norm=norm,
        marker='o',
        zorder=2
    )

    # Add a colorbar
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Nominal K wt %")


    ax.set_yscale('log')


    ax.set_xlabel("CFD KOH [um$^{-3}$]")

axes[0].set_ylabel("$\Delta AS$ Maximum")
axes[1].set_ylabel("Delta PD1 [mV]")
axes[0].set_ylim(5e-3,1.1e-1)
axes[1].set_ylim(3.5e-1,2.5e0)



plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_ionization_vsCFDKOH.png'), dpi=300, bbox_inches='tight')