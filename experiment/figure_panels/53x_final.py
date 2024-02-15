#%%[markdown]

# # 53x (seedramp) Final figure panels

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.fileio.ct import load_df_cuttimes
from mhdpy.coords.ct import downselect_acq_time

tc = '53x'

#%%

# Load Absem Data 

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem = ds_absem.drop(0, 'kwt')

#%%

# # Load MWS Data

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

ds_lecroy = ds_lecroy.drop(0,'kwt')

#%%

# Load cfd line profiles

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)

ds_cfd = ds_cfd.interp(kwt=ds_lecroy.coords['kwt']).dropna('kwt', how='all')

ds_cfd = ds_cfd.cfd.convert_species_rho()

goldi_pos = ds_cfd['x'].min().item() + 0.18
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')


# Calculate kr

# combining k4 from 
# 1. Axford, S.D.T., and Hayhurst, A.N. (1997). Mass spectrometric sampling of negative ions from flames of hydrogen and oxygen: the kinetics of electron attachment and detachment in hot mixtures of H2O, O2, OH and HO2. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences 452, 1007–1033. 10.1098/rspa.1996.0051.

# Weighted average of collison partners for O2 rate
k4_species = {
    'N2' : Quantity(1e-32, 'cm^6/s'),
    'H2O': Quantity(8e-30, 'cm^6/s'),
    'CO2': Quantity(1e-31, 'cm^6/s'), 
    'O2' : Quantity(2e-30, 'cm^6/s'),
}

weighted_avg = 0

for species, k4 in k4_species.items():
    weighted_avg += k4*ds_cfd[species]

rxn_rates = {
    'K+':  Quantity(4e-24, 'K*cm^6/s')*(1/ds_cfd['T'])*ds_cfd['rho'],
    'OH': Quantity(3e-31, 'cm^6/s')*ds_cfd['rho'],
    'O2': weighted_avg,
    'H2O': Quantity(1.6e-6, 'cm^3/s')*np.exp(-(Quantity(36060, 'K')/ds_cfd['T'])),
}

ds_kr = xr.Dataset(rxn_rates).pint.to('cm^3/s')

# Calculate expected tau for recombination with those species

ds_species_cfd = ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O']]

ds_species_cfd = ds_species_cfd.pint.to('1/cm^3')

ds_species_cfd = ds_species_cfd.rename({'Yeq_K+': 'K+', 'Yeq_OH': 'OH'})
ds_tau = 1/(ds_species_cfd*ds_kr)

ds_tau = ds_tau.pint.to('us')

ds_tau


#%%

# Downselect times to those in the CFD input kwt cases 
fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)
ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs)

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit_absem = ds_absem
ds_fit_absem, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_2(ds_fit_absem)

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_fit = ds_lecroy.mean('mnum')
da_fit_lecroy = ds_fit.mws.calc_mag_phase_AS()['AS']
ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_exp(da_fit_lecroy, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))

#%%

da_stats = ds_lecroy['AS'].sel(time=slice(-1,1))
mws_max = da_stats.mean('mnum').max('time')
mws_std = da_stats.std('mnum').max('time')

ds_params = xr.merge([
    ds_p_absem['nK_m3'].sel(mp='barrel').drop('mp').rename('nK_m3_barrel'),
    ds_p_absem['nK_m3'].sel(mp='mw_horns').drop('mp').rename('nK_m3_mw_horns'),
    mws_max.rename('AS_max'),
    mws_std.rename('AS_std'),
    ds_p_mws['decay'].rename('mws_fit_decay')
    ])


# Perform averages over run
# TODO: combine with wma accessor

count = ds_params['nK_m3_barrel'].count('run')

ds_p_stats = xr.Dataset(coords=count.coords)

for var in ds_params.data_vars:
    mean = ds_params[var].mean('run')
    std = ds_params[var].std('run')
    stderr = std/np.sqrt(count)

    ds_p_stats["{}_mean".format(var)] = mean
    ds_p_stats["{}_std".format(var)] = std
    ds_p_stats["{}_stderr".format(var)] = stderr
# %%


vars = ['nK_m3_barrel', 'nK_m3_mw_horns', 'AS_max', 'mws_fit_decay']

fig, axes = plt.subplots(3, 1, figsize=(5,10), sharex=True)

var = 'nK_m3_barrel'
axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='Barrel'
    )

var = 'nK_m3_mw_horns'
axes[0].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_stderr'.format(var)], 
    marker='o', capsize=5,
    label='MW Horns'
    )


var = 'AS_max'
axes[1].errorbar(
    ds_p_stats.coords['kwt'], 
    ds_p_stats['{}_mean'.format(var)], 
    yerr=ds_p_stats['{}_std'.format(var)], 
    marker='o', capsize=5
    )


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

axes[2].legend()

axes[0].set_ylabel("nK_m3 [#/m^3]")
axes[1].set_ylabel("AS Maximum")
axes[2].set_ylabel("Time Constant [us]")


for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')


axes[-1].set_xlabel("K wt % nominal")

plt.savefig(pjoin(DIR_FIG_OUT, '53x_params.png'), dpi=300, bbox_inches='tight')