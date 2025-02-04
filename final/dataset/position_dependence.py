#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu
from mhdlab.analysis import mwt, absem
from mhdlab.plot import dropna

stack_dims = ['date','run_num','motor','phi']

#%%
tc = '536_pos'
ds_absem_536_pos = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem_536_pos = ds_absem_536_pos.expand_dims('phi').stack(temp=stack_dims)

ds_lecroy_536_pos = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy_536_pos = ds_lecroy_536_pos.sortby('time').expand_dims('phi').stack(temp=stack_dims)

#%%
tc = '5x6_pos'
ds_lecroy_5x6 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy_5x6 = ds_lecroy_5x6.sel(phi=[0.79], method='nearest')
ds_lecroy_5x6 = ds_lecroy_5x6[['pd1','pd2']].stack(temp=stack_dims)

#%%
tc = '516_pos'
ds_absem_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem_516 = ds_absem_516.expand_dims('phi').stack(temp=stack_dims)

ds_lecroy_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy_516 = ds_lecroy_516.sortby('time').expand_dims('phi').stack(temp=stack_dims)

#%%
ds_absem = xr.concat([ds_absem_536_pos, ds_absem_516], dim='temp').unstack('temp')
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.absem.calc_alpha().sel(wavelength=slice(750,790))
ds_absem.mean('mnum').unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_absem.cdf'))

#%%
da_mwt_nothing = ppu.fileio.load_mwt_T0()
ds_lecroy = xr.concat([ds_lecroy_536_pos, ds_lecroy_5x6, ds_lecroy_516], dim='temp').unstack('temp')
ds_lecroy = ds_lecroy.mwt.calc_AS_abs(mag_0=da_mwt_nothing).mwt.calc_dPD()
ds_lecroy = ds_lecroy.xr_utils.stack_run()

#%%
ds_stats_lecroy = ds_lecroy.mwt.calc_time_stats()[['dAS_abs_max', 'mag_pp', 'mag_fluct', 'SFR_abs', 'dpd1_max']]
ds_stats_lecroy = ds_stats_lecroy.mean('mnum', keep_attrs=True)

# Now that statistics are calculated we average over mnum and unstack.
ds_lecroy = ds_lecroy.mean('mnum')
ds_lecroy.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_lecroy.cdf'))

#%%
ds_cfd_beam = ppu.fileio.load_cfd_beam()
ds_cfd_beam = ds_cfd_beam.sel(offset=0).sel(kwt=1)
ds_cfd_beam = ds_cfd_beam.assign_coords(phi=ds_absem.coords['phi'].values)

#TODO: motor coordinates not exactly the same, why?
ds_cfd_beam = ds_cfd_beam.sel(motor=ds_absem.coords['motor'].values, method='nearest')
ds_cfd_beam = ds_cfd_beam.assign_coords(motor=ds_absem.coords['motor'].values)
da_cfd_beam = ds_cfd_beam[ppu.CFD_K_SPECIES_NAME] / ds_cfd_beam[ppu.CFD_K_SPECIES_NAME].max('dist')

#%%

#%%

# Fit the position dependent mw horns data
# 

from mhdlab.analysis.absem.fitting import pipe_fit_alpha_num_1

ds_fit = ds_absem.mean('mnum').sel(mp='mw_horns').dropna('run', how='all')
nK_profile = da_cfd_beam.sel(mp='mw_horns')
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")


ds_alpha_fit.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_alpha_fit.cdf'))

nK_mw_horns = ds_p['nK_m3'].rename('nK_mw_horns')
# nK_barrel = ds_p['nK_m3'].sel(mp='barrel').rename('nK_barrel').drop_vars('mp')

ds_stats_out = xr.merge([nK_mw_horns, ds_stats_lecroy]).sortby('motor').dropna('run', how='all')
ds_stats_out.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_stats_pos.cdf'))


#%%
da_plot = ds.to_array('var').sel(phi=0.6, method='nearest').mean('run').dropna('motor',how='all')
motor_vals_06 = da_plot.coords['motor'].values
da_plot.attrs['long_name'] = '$\\alpha$'
da_plot.plot(hue='var', row='motor', figsize=(8,12))
plt.ylim(-0.1,1.1)
plt.savefig(pjoin(DIR_FIG_OUT, '516_pos_phi06_absem_fit.png'))

#%%
da_plot = ds.to_array('var').sel(phi=0.8, method='nearest').mean('run').dropna('motor',how='all')
da_plot = da_plot.sel(motor=motor_vals_06)
da_plot.attrs['long_name'] = '$\\alpha$'
da_plot.plot(hue='var', row='motor', figsize=(8,12))
plt.ylim(-0.1,1.1)
plt.savefig(pjoin(DIR_FIG_OUT, '536_pos_phi08_absem_fit.png'))

#%%

# Now fit the barrel data. Both datasets have it stored with 'motor' dimension.
# Data is duplicated for cfd data  (see load_cfd_beam fn)
# For absem, this corresponds to the times the stage was at this position...so we average over all motor


ds_fit = ds_absem.mean('mnum').sel(mp='barrel').dropna('run', how='all')
ds_fit = ds_fit.mean('motor', keep_attrs=True)

nK_profile = da_cfd_beam.sel(mp='barrel').isel(motor=0)

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds_p_cfdprofile = ds_p['nK_m3'].rename('nK_barrel_cfdprofile')

#%%

# was originally combining with tophat profile, but removed now. see position_dependence_extra_aes.py
ds_p_out = xr.merge([ds_p_cfdprofile]).dropna('run', how='all')

ds_p_out.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_barrel.cdf'))


