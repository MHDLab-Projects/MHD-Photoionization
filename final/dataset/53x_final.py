#%%[markdown]
# # 53x (seedramp) Final Dataset

#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdlab.fileio.ct import load_df_cuttimes
from mhdlab.coords.ct import downselect_acq_time
from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted, calc_krm
from mhdlab.analysis.absem.fitting import pipe_fit_alpha_num_1
from mhdlab.analysis.mwt.fitting import pipe_fit_mwt_3, pipe_fit_exp

#%%
tc = '53x'
ds_absem = ppu.fileio.load_absem(tc)
ds_lecroy = ppu.fileio.load_lecroy(tc, AS_calc='absolute')
ds_cfd = ppu.fileio.load_cfd_centerline(kwt_interp=ds_lecroy.coords['kwt'].values)
ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)

#%%
ds_lecroy['i'].mean('mnum').mean('run').plot(hue='kwt', marker='o')
plt.xlim(-1, 1)

#%%
goldi_pos = Quantity(180, 'mm')
ds_cfd = ds_cfd.sel(x=goldi_pos, method='nearest')

#%%
ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()
ds_krb = gen_ds_krb(ds_cfd['T'], ds_cfd['rho_number'])


ds_krb['O2_S'] = calc_krbO2_weighted(ds_cfd)

#%%
ds_species_cfd = ds_cfd[[ppu.CFD_K_SPECIES_NAME, ppu.CFD_KOH_SPECIES_NAME, ppu.CFD_Kp_SPECIES_NAME, ppu.CFD_allK_SPECIES_NAME, ppu.CFD_OH_SPECIES_NAME, 'O2', 'H2O', 'T', 'rho_number', 'N2', 'CO2']]
ds_species_cfd = ds_species_cfd.rename({ppu.CFD_Kp_SPECIES_NAME: 'K+', ppu.CFD_OH_SPECIES_NAME: 'OH'})
ds_krm = calc_krm(ds_krb, ds_species_cfd)
ds_tau = (1 / ds_krm).pint.to('us')

#%%
ds_tau.pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, 'ds_tau.cdf'))
ds_species_cfd.pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_species_cfd.cdf'))

#%%
fp_ct_seedramp = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)
ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)
ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs)
ds_lecroy.mean('mnum').unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_lecroy.cdf'))
ds_absem.mean('mnum').unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_absem.cdf'))

#%%
ds_cfd_beam = ppu.fileio.load_cfd_beam(kwt_interp=ds_absem.coords['kwt'].values)
ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0).sel(motor=goldi_pos, method='nearest')
da_cfd_beam = ds_cfd_beam[ppu.CFD_K_SPECIES_NAME]
da_cfd_beam = da_cfd_beam / da_cfd_beam.max('dist')

#%%
ds_fit = ds_absem.mean('mnum')
ds_fit = ds_fit.sel(kwt=da_cfd_beam.kwt.values)
ds_fit_absem, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile': da_cfd_beam})
ds_fit_absem['alpha'] = ds_fit['alpha']

#%%
da_plot = ds_fit_absem.mean('run').to_array('var')
da_plot.attrs['long_name'] = '$\\alpha$'
da_plot.plot.line(row='kwt', col='mp', hue='var', figsize=(8, 12))
plt.savefig(pjoin(DIR_FIG_OUT, '53x_absem_fit.png'))

#%%
ds_fit = ds_lecroy
da_fit_lecroy = ds_fit['dAS_abs'].mean('mnum')
ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_mwt_3(da_fit_lecroy)
ds_p_mwt['decay'] = 1 / ds_p_mwt['krm']
ds_fit_mwt.unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_fit_mwt_dnedt.cdf'))
ds_p_dnedt = ds_p_mwt.copy()

#%%
da_fit_lecroy.sel(kwt=[0.05, 0.99], method='nearest').plot(hue='run_plot', row='kwt', x='time')
plt.xlim(-1, 1)

#%%
ds_fit_mwt.mean('run')[['AS_sel', 'AS_fit', 'AS_all']].to_array('var').plot(col='kwt', col_wrap=2, hue='var', figsize=(5, 10))
plt.yscale('log')
plt.xlim(-1, 40)
plt.ylim(1e-3, 1.1)
plt.savefig(pjoin(DIR_FIG_OUT, '53x_mwt_fit_dnedtv2.png'))

#%%
ds_fit = ds_lecroy
da_fit_lecroy = ds_fit['dAS_abs'].mean('mnum')
ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_exp(da_fit_lecroy)
ds_p_mwt['decay'] = ds_p_mwt['decay'].pint.quantify('us')
ds_p_exp = ds_p_mwt.copy()
ds_fit_mwt.unstack('run').pint.dequantify().to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_fit_mwt_exp.cdf'))

#%%
da_sel = ds_lecroy['AS_abs'].mean('mnum').sel(kwt=0.05, method='nearest').dropna('run', how='all').isel(run=2)
da_mean = da_sel
da_std = da_sel
da_mean.plot()
plt.fill_between(da_mean.time, da_mean - da_std, da_mean + da_std, alpha=0.5)
plt.yscale('log')
plt.ylim(1e-3, 1e-2)
plt.xlim(-1, 10)

#%%
ds_fit_mwt.mean('run')[['AS_all', 'AS_fit', 'AS_all']].to_array('var').plot(col='kwt', col_wrap=2, hue='var', figsize=(5, 10))
plt.yscale('log')
plt.xlim(-1, 40)
plt.ylim(1e-4, 1e-1)
plt.savefig(pjoin(DIR_FIG_OUT, '53x_mwt_fit_exp.png'))

#%%
ds_stats_lecroy = ds_lecroy.mwt.calc_time_stats()[['dAS_abs_max', 'mag_pp', 'mag_fluct', 'SFR_abs', 'dpd1_max']]
ds_stats_lecroy = ds_stats_lecroy.mean('mnum', keep_attrs=True)

#%%
ds_p_exp['krm'] = 1 / ds_p_exp['decay']
ds_params = xr.merge([
    ds_p_absem['nK_m3'].sel(mp='barrel').drop_vars('mp').rename('nK_m3_barrel'),
    ds_p_absem['nK_m3'].sel(mp='mw_horns').drop_vars('mp').rename('nK_m3_mw_horns'),
    ds_p_dnedt['decay'].rename('mwt_fit_decay'),
    ds_p_exp['decay'].rename('mwt_fit_decay_exp'),
    ds_p_exp['krm'].rename('mwt_fit_krm'),
    ds_p_dnedt['dne'].rename('mwt_fit_dne'),
    ds_stats_lecroy
])

count = ds_params['nK_m3_barrel'].count('run')
ds_p_stats = xr.Dataset(coords=count.coords)

for var in ds_params.data_vars:
    mean = ds_params[var].mean('run', keep_attrs=True)
    std = ds_params[var].std('run', keep_attrs=True)
    stderr = std / np.sqrt(count)
    stderr.attrs = std.attrs
    ds_p_stats["{}_mean".format(var)] = mean
    ds_p_stats["{}_std".format(var)] = std
    ds_p_stats["{}_stderr".format(var)] = stderr

ds_p_stats = ds_p_stats.pint.dequantify()
ds_p_stats.to_netcdf(pjoin(DIR_DATA_OUT, '53x_ds_p_stats.cdf'))

