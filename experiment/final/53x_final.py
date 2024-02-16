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

ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

goldi_pos = ds_cfd['x'].min().item() + 0.18
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')


# Calculate kr

ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()

from mhdpi_utils.kinetics import get_kinetics

ds_kr = get_kinetics(ds_cfd)

# Calculate expected tau for recombination with those species

ds_species_cfd = ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O', 'Yeq_KOH', 'Yeq_K']]

ds_species_cfd = ds_species_cfd.pint.to('1/cm^3')

ds_species_cfd = ds_species_cfd.rename({'Yeq_K+': 'K+', 'Yeq_OH': 'OH'})

das = []
for var in ds_kr.data_vars:
    species_name = var.split('_')[0]
    da_tau = 1/(ds_species_cfd[species_name]*ds_kr[var])
    da_tau = da_tau.pint.to('us')
    da_tau = da_tau.rename('{}'.format(var))
    das.append(da_tau)

ds_tau = xr.merge(das)

ds_tau.pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, 'ds_tau.cdf'))

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

delta_pd1 = ds_lecroy['pd1'] - ds_lecroy['pd1'].sel(time=slice(-1,0)).mean('time')
delta_pd1 = delta_pd1.dropna('run', how='all')
delta_pd1 = delta_pd1.mean('mnum').max('time')
delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')


ds_params = xr.merge([
    ds_p_absem['nK_m3'].sel(mp='barrel').drop('mp').rename('nK_m3_barrel'),
    ds_p_absem['nK_m3'].sel(mp='mw_horns').drop('mp').rename('nK_m3_mw_horns'),
    mws_max.rename('AS_max'),
    mws_std.rename('AS_std'),
    ds_p_mws['decay'].rename('mws_fit_decay'),
    delta_pd1.rename('delta_pd1')
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
