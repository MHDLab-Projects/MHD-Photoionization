#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu #Sets matplotlib style
import matplotlib.ticker as mticker



# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()

ds_mwt_fit_exp = xr.open_dataset(pjoin(data_directory, '53x_ds_fit_mwt_exp.cdf')).pint.quantify()
ds_mwt_fit_dnedt = xr.open_dataset(pjoin(data_directory, '53x_ds_fit_mwt_dnedt.cdf')).pint.quantify()


# Load in the raw lecroy data for confidence intervals. Doing this here to avoid having to save raw mnum acquisitions...Revisit. 
from mhdlab.fileio.ct import load_df_cuttimes
from mhdlab.coords.ct import downselect_acq_time
ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute')
# Downselect times to those in the CFD input kwt cases 
fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)


#TODO: the time grids are not the same for the two fitting methods. Why?
ds_mwt_fit_exp = ds_mwt_fit_exp.interp(time=ds_mwt_fit_dnedt.time, method='linear')
#TODO: realized 'method' cannot be used as a coordinate as it is a reserved word in xarray(i.e. method='nearest'). Need to change this in the fitting code.
ds_mwt_fit = xr.concat([
    ds_mwt_fit_exp.assign_coords(fit_method=['exp']),
    ds_mwt_fit_dnedt.assign_coords(fit_method=['dnedt'])
    ], 'fit_method')

# ds_mwt_fit = ds_mwt_fit.assign_coords(method=[str(s) for s in ds_mwt_fit.coords['method'].values])

#%%

fig, axes = plt.subplots(2, 1, figsize=(3.5,6), sharex=True)

var = 'nK_m3_barrel'
data=ds_p_stats['{}_mean'.format(var)].pint.to('particle/cm^3').pint.dequantify()
yerr=ds_p_stats['{}_std'.format(var)].pint.to('particle/cm^3').pint.dequantify()
line_nK_barrel = axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    data,
    yerr=yerr,
    marker='o', capsize=5,
    label='Barrel'
    )

var = 'nK_m3_mw_horns'
data=ds_p_stats['{}_mean'.format(var)].pint.to('particle/cm^3').pint.dequantify()
yerr=ds_p_stats['{}_std'.format(var)].pint.to('particle/cm^3').pint.dequantify()
line_nK_mwhorns = axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    data,
    yerr=yerr,
    marker='o', capsize=5,
    label='MW Horns'
    )

da_KOH = ds_species_cfd[ppu.CFD_KOH_SPECIES_NAME].pint.to('particle/cm**3').pint.dequantify()
lineKOH = da_KOH.plot(ax=axes[0], label='CFD: KOH', ls = "--")
da_nK = ds_species_cfd[ppu.CFD_K_SPECIES_NAME].pint.to('particle/cm**3').pint.dequantify()
linenK = da_nK.plot(ax=axes[0], label='CFD: K')
# da_allK = ds_species_cfd[ppu.CFD_allK_SPECIES_NAME].pint.to('particle/cm**3').pint.dequantify()
# line_allK = da_allK.plot(ax=axes[0], label='CFD: All K')


axes[0].set_ylabel(r"Species Concentration [$\mathrm{cm^{-3}}$]")
axes[0].legend(
    [line_nK_barrel, line_nK_mwhorns, lineKOH[0], linenK[0]], 
    ['Expt. $n_K$ (Barrel)', 'Expt. $n_K$ (180 mm)', 'CFD KOH (180 mm)', 'CFD K (180 mm)'],
    bbox_to_anchor=(0.65, 0.75), loc='upper left', framealpha=1
    )

axes[0].set_title('')
axes[0].set_xlabel('')

var = 'dAS_abs_max'
lineAS = axes[1].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label=r'$\Delta AS_{max}$',
    color='darkgreen'
    )

axes[1].set_ylabel(r"$\Delta AS_{max}$")


ta = axes[1].twinx()
#No standard deviation for delta_pd1 as only one run
var = 'dpd1_max'
linePD = ta.plot(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    marker='x',
    label='$\\Delta PD_{max}$',
    color='darkblue'
    )

axes[1].set_ylim(1e-2,2e-1)
axes[1].legend([lineAS, linePD[0]], ['$\\Delta AS_{max}$', '$\\Delta PD_{max}$'])

ta.set_ylabel("$\\Delta PD_{max}$ [mV]")
ta.grid(False)

ta.set_yscale('log')


ta.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ta.yaxis.set_major_formatter(mticker.ScalarFormatter())
ta.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))

ta.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.2, 1.4, 1.6], minor=True)
ta.yaxis.set_minor_formatter(mticker.ScalarFormatter())
ta.yaxis.set_minor_formatter(mticker.FixedFormatter(['', '0.6', '', '0.8', '', '1.2', '1.4', '1.6']))

fig.subplots_adjust(hspace = 0)
for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')



axes[-1].set_xlabel("$K_{wt,nominal} [\\%]$")

labels = ['A)','B)','C)','D)']

for ax, label in zip(axes, labels):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - .25, Y - 0.015, label)
    # ax.tick_params(axis='y', which = "both", colors='white')



fig.savefig(pjoin(REPO_DIR, 'final','figures', 'Fig6_kwt_dependence.svg'), dpi=300, bbox_inches='tight')
# %%


fig, axes = plt.subplots(2,1, figsize=(3.5,6))


run_sel = ('2023-05-12', 1)



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
axes[0].plot(x, mean, label='AS')

da_plot_fit = ds_mwt_fit.sel(date=run_sel[0]).sel(run_num=run_sel[1]).sel(kwt=1, method='nearest').sel(fit_method='exp')['AS_fit']
da_plot_fit.plot(label='AS Fit', ax = axes[0])

# Plot confidence interval
axes[0].fill_between(x, mean - std_err, mean + std_err, color='b', alpha=0.1)


axes[0].set_yscale('log')
axes[0].set_xlim(-1,25)
axes[0].set_ylim(1e-5,)
axes[0].set_ylabel(r'$\Delta AS$')
axes[0].set_xlabel(r'Time [$\mathrm{\mu s}$]')

axes[0].set_title('')
axes[0].legend(['Data', 'Fit'])

# plt.savefig(pjoin(DIR_FIG_OUT, '53x_mwt_fit_exp.png'), dpi=300, bbox_inches='tight')




ax = axes[1]

var = 'mwt_fit_decay_exp'
ax.errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5,
    label=r'$\tau_{exp}$'
    )

ds_tau_plot = ds_tau.drop_vars('O2_exp_eff') # This is included for the viability analyis, is calculated from experimental data so is redundant with the lifetime that these data are eventually compared to. 

print(f"Removing H2O from 53x_params_recomb_exp plot for clarity. Average lifetime: {ds_tau_plot['H2O'].mean().values} us")
ds_tau_plot = ds_tau_plot.drop_vars('H2O') # Value is ~ 10^5, not including in plot for clarity.

label_map = {
    'O2_A': 'O_{2,A}',
    'O2_G': 'O_{2,G}',
    'O2_S': 'O_{2,S}',
    'K+': 'K^+',
}

for species in ds_tau_plot.data_vars:
    if species.startswith('O2'):
        marker = 'x'
    elif species.startswith('OH'):
        marker = '.'
    elif species.startswith('K'):
        marker = 's'

    if species in label_map:
        label = label_map[species]
    else:
        label = species

    label = r"$\tau_{pred}: " + label + "$"
    # label = r'${}$'.format(label)

    ds_tau_plot[species].plot(label=label, ax=axes[1], marker=marker)

axes[1].legend(bbox_to_anchor=(1, 1), loc='upper left', framealpha=1) 
axes[1].set_ylabel(r"Time Constant [$\mathrm{\mu s}$]")
axes[1].set_title('')
axes[1].set_ylim(5e-2, 2e3)

axes[1].set_xscale('log')
axes[1].set_yscale('log')

axes[1].set_xlabel("$K_{wt,nominal} [\\%]$")

# plt.savefig(pjoin(DIR_FIG_OUT, '53x_params_recomb_exp.png'), dpi=300, bbox_inches='tight')


for ax, label in zip(axes, labels):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - .2, Y - 0.015, label)
    # ax.tick_params(axis='y', which = "both", colors='white')


fig.savefig(os.path.join(REPO_DIR, 'final','figures', 'Fig7_kwt_dependence.svg'))

# %%