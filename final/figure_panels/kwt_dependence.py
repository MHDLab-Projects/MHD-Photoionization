#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu #Sets matplotlib style

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

#%%

run_sel = ('2023-05-12', 1)

plt.figure(figsize=(5,2))

da_plot = ds_lecroy['dAS_abs'].sel(run=run_sel).sel(kwt=1, method='nearest')
da_plot

# Calculate mean and standard deviation
mean = da_plot.mean('mnum')
std = da_plot.std('mnum')
count = da_plot.count('mnum')
std_err = std / np.sqrt(count)

print("Number of acquisitions: {}".format(count.isel(time=0).values))

# Get x values (assuming they are the same for mean and std)
x = mean.coords[mean.dims[0]].pint.magnitude
mean = mean.pint.magnitude
std_err = std_err.pint.magnitude
# Plot mean
plt.plot(x, mean, label='AS')

da_plot_fit = ds_mws_fit.sel(date=run_sel[0]).sel(run_num=run_sel[1]).sel(kwt=1, method='nearest').sel(fit_method='exp')['AS_fit']
da_plot_fit.plot(label='AS Fit')

# Plot confidence interval
plt.fill_between(x, mean - std_err, mean + std_err, color='b', alpha=0.1)


plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5,)
plt.ylabel('$\Delta AS$ [dimensionless]')

plt.title('')
plt.legend(['Data', 'Fit'])

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

da_KOH = ds_species_cfd[ppu.CFD_KOH_SPECIES_NAME].pint.to('particle/m**3').pint.dequantify()
lineKOH = da_KOH.plot(ax=axes[0], label='CFD: KOH')
da_nK = ds_species_cfd[ppu.CFD_K_SPECIES_NAME].pint.to('particle/m**3').pint.dequantify()
linenK = da_nK.plot(ax=axes[0], label='CFD: K')
da_allK = ds_species_cfd[ppu.CFD_allK_SPECIES_NAME].pint.to('particle/m**3').pint.dequantify()
line_allK = da_allK.plot(ax=axes[0], label='CFD: All K')


axes[0].set_ylabel("Species Concentration [$\#/m^3$]")
axes[0].legend(
    [line_nK_barrel, line_nK_mwhorns, lineKOH[0], linenK[0], line_allK[0]], 
    ['Expt. $n_K$ Barrel', 'Expt $n_K$ 180 mm', 'CFD KOH (180 mm)', 'CFD K (180 mm)' , 'CFD All K (180 mm)'],
    bbox_to_anchor=(0.75, 0.85), loc='upper left', framealpha=1
    )

axes[0].set_title('')
axes[0].set_xlabel('')

var = 'dAS_abs_max'
lineAS = axes[1].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='$\Delta AS$ Maximum',
    color='darkgreen'
    )

axes[1].set_ylabel("$\Delta AS$ Maximum")


ta = axes[1].twinx()
#No standard deviation for delta_pd1 as only one run
var = 'dpd1_max'
linePD = ta.plot(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    marker='x',
    label='Delta PD1',
    color='darkblue'
    )

axes[1].set_ylim(1e-2,2e-1)
axes[1].legend([lineAS, linePD[0]], ['AS Maximum', 'Delta PD1'])

ta.set_ylabel("Delta PD1 [mV]")

ta.set_yscale('log')
ta.set_ylim(3.5e-1,2.5e0)



for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')

axes[-1].set_xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_ionization.png'), dpi=300, bbox_inches='tight')
# %%

# %%


fig, ax = plt.subplots(figsize=(4.5,2.5))

var = 'mws_fit_decay_exp'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label='$\Delta AS$ Fit'
    )

for species in ds_tau.data_vars:
    ds_tau[species].plot(label="{}".format(species), ax=ax)

ax.legend(bbox_to_anchor=(1, 1), loc='upper left', framealpha=1) 
ax.set_ylabel("Time Constant [us]")
ax.set_title('')
ax.set_ylim(1e-1, 1e4)

ax.set_xscale('log')
ax.set_yscale('log')

plt.xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_recomb_exp.png'), dpi=300, bbox_inches='tight')
# %%
