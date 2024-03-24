#%%[markdown]

# # 53x (seedramp) Final Dataset


#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.fileio.ct import load_df_cuttimes
from mhdpy.coords.ct import downselect_acq_time


#%%

tc = '53x'
# Load Absem Data 
ds_absem = ppu.fileio.load_absem(tc)

# # Load MWS Data
# TODO: having to 'cancel' AS calulation then perform again to avoid memory errors when downselecting. Revist
ds_lecroy = ppu.fileio.load_lecroy(tc, AS_calc='absolute')

# Load cfd line profiles
ds_cfd = ppu.fileio.load_cfd_centerline(kwt_interp=ds_lecroy.coords['kwt'].values)

ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)

#%%

ds_lecroy
#%%


ds_lecroy['i'].mean('mnum').mean('run').plot(hue='kwt',marker='o')
plt.xlim(-1,1)
#%%

goldi_pos =  Quantity(178, 'mm')
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')

#%%
# Calculate kr

ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()

from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted

ds_krb = gen_ds_krb(ds_cfd['T'], ds_cfd['rho_number'])

ds_krb['O2_C'] = calc_krbO2_weighted(ds_cfd)

ds_krb


#%%

# Calculate expected tau for recombination with those species

ds_species_cfd = ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O', 'Yeq_KOH', 'Yeq_K', 'all_K_Yeq', 'T', 'rho_number', 'N2', 'CO2']]
ds_species_cfd = ds_species_cfd.rename({'Yeq_K+': 'K+', 'Yeq_OH': 'OH'})

from pi_paper_utils.kinetics import calc_krm

ds_krm = calc_krm(ds_krb, ds_species_cfd)

ds_tau = (1/ds_krm).pint.to('us')


#%%
ds_tau.pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, 'ds_tau.cdf'))

# ds_species_cfd = ds_species_cfd.pint.to('particle/ml')
ds_species_cfd.pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_species_cfd.cdf'))

#%%

# Downselect times to those in the CFD input kwt cases 
fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)
# mag_0 = ppu.fileio.load_mws_T0()

# ds_lecroy = ds_lecroy.mean('mnum')
# ds_lecroy = ds_lecroy.unstack('run').mws.calc_AS_abs(mag_0=mag_0).xr_utils.stack_run()

ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs)

ds_lecroy.mean('mnum').unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_lecroy.cdf'))
ds_absem.mean('mnum').unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_absem.cdf'))

#%%


#%%

ds_cfd_beam = ppu.fileio.load_cfd_beam(kwt_interp=ds_absem.coords['kwt'].values)

ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0).sel(motor=goldi_pos,method='nearest')

da_cfd_beam = ds_cfd_beam['Yeq_K']
da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

da_cfd_beam

#%%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_num_1

ds_fit = ds_absem.mean('mnum')
ds_fit = ds_fit.sel(kwt= da_cfd_beam.kwt.values) #TODO: downselecting as we don't have cfd for all kwt. Remove once we do
ds_fit_absem, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam})
ds_fit_absem['alpha'] = ds_fit['alpha']
#%%

da_plot = ds_fit_absem.mean('run').to_array('var')

da_plot.attrs['long_name'] = '$\\alpha$'

da_plot.plot.line(row='kwt', col='mp', hue='var')

plt.savefig(pjoin(DIR_FIG_OUT, '53x_absem_fit.png'))

# ds_p_absem['nK_m3'].sel(mp='barrel').mean('run').plot.line(x='kwt')

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3

ds_fit = ds_lecroy
da_fit_lecroy = ds_fit['dAS_abs'].mean('mnum')


# ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_exp(da_fit_lecroy, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))
ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_mws_3(da_fit_lecroy)
ds_p_mws['decay'] = 1/ds_p_mws['krm']

ds_fit_mws.unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_fit_mws_dnedt.cdf'))
ds_p_dnedt = ds_p_mws.copy()

#%%

da_fit_lecroy.sel(kwt=[0.05, 0.99], method='nearest').plot(hue='run_plot', row='kwt', x='time')

plt.xlim(-1,1)

# plt.gca().get_legend().set_bbox_to_anchor((1,0.5))



#%%

ds_fit_mws.mean('run')[['AS_sel','AS_fit','AS_all']].to_array('var').plot(col='kwt', col_wrap=2, hue='var', figsize=(5,10))

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-3, 1.1)

plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_dnedtv2.png'))

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_fit = ds_lecroy
da_fit_lecroy = ds_fit['dAS_abs'].mean('mnum')

ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_exp(da_fit_lecroy)

ds_p_exp = ds_p_mws.copy()

ds_fit_mws.unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_fit_mws_exp.cdf'))

#%%

#TODO: output a standard plot of fits for SI or save to file

da_sel = ds_lecroy['AS_abs'].mean('mnum').sel(kwt=0.05,method='nearest').dropna('run',how='all').isel(run=2)

da_mean = da_sel
da_std = da_sel

da_mean.plot()
plt.fill_between(da_mean.time, da_mean-da_std, da_mean+da_std, alpha=0.5)

plt.yscale('log')
plt.ylim(1e-3, 1e-2)
plt.xlim(-1,10)

#%%

ds_fit_mws.mean('run')[['AS_all','AS_fit','AS_all']].to_array('var').plot(col='kwt', col_wrap=2, hue='var', figsize=(5,10))

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-4, 1e-1)

plt.savefig(pjoin(DIR_FIG_OUT, '53x_mws_fit_exp.png'))

#%%



ds_stats_lecroy = ds_lecroy.mws.calc_time_stats()[['dAS_abs_max', 'mag_pp', 'mag_fluct', 'SFR_abs', 'dpd1_max']]
ds_stats_lecroy = ds_stats_lecroy.mean('mnum', keep_attrs=True)

#%%

ds_params = xr.merge([
    ds_p_absem['nK_m3'].sel(mp='barrel').drop('mp').rename('nK_m3_barrel'),
    ds_p_absem['nK_m3'].sel(mp='mw_horns').drop('mp').rename('nK_m3_mw_horns'),
    ds_p_dnedt['decay'].rename('mws_fit_decay'),
    ds_p_exp['decay'].rename('mws_fit_decay_exp'),
    ds_p_dnedt['dne'].rename('mws_fit_dne'),
    ds_stats_lecroy
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


ds_p_stats = ds_p_stats.pint.dequantify()

ds_p_stats.to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_p_stats.cdf'))

# %%
