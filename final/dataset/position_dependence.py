#%%


from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

stack_dims = ['date','run_num','motor','phi']

# %%

tc = '536_pos'

ds_absem_536_pos = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))

ds_absem_536_pos = ds_absem_536_pos.expand_dims('phi').stack(temp = stack_dims)

ds_lecroy_536_pos = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy_536_pos = ds_lecroy_536_pos.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy_536_pos = ds_lecroy_536_pos.expand_dims('phi').stack(temp = stack_dims)
# ds_lecroy_536_pos = ds_lecroy_536_pos.mean('mnum')

#%%


#%%

# add 536 photodiode data from 5x6_pos run. Do not have photodiode data for 536_pos run

tc = '5x6_pos'

ds_lecroy_5x6 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy_5x6 = ds_lecroy_5x6.sel(phi=[0.79], method='nearest') # Cannot add 0.65 phi, because it conflicts with 516_pos. Need to change run_num for this measurement?

ds_lecroy_5x6 = ds_lecroy_5x6[['pd1','pd2']]
ds_lecroy_5x6 = ds_lecroy_5x6.stack(temp = stack_dims)
# ds_lecroy_5x6 = ds_lecroy_5x6.mean('mnum')

#%%

tc = '516_pos'

ds_absem_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem_516 = ds_absem_516.expand_dims('phi') # allows merging with 5x6_pos
ds_absem_516 = ds_absem_516.stack(temp = stack_dims)

ds_lecroy_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy_516 = ds_lecroy_516.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy_516 = ds_lecroy_516.expand_dims('phi') # allows merging with 5x6_pos
ds_lecroy_516 = ds_lecroy_516.stack(temp = stack_dims)  
# ds_lecroy_516 = ds_lecroy_516.mean('mnum')

#%%

ds_absem = xr.concat([ds_absem_536_pos, ds_absem_516], dim='temp').unstack('temp')
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem.mean('mnum').unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_absem.cdf'))

#%%

da_mws_nothing  = ppu.fileio.load_mws_T0()

ds_lecroy = xr.concat([ds_lecroy_536_pos, ds_lecroy_5x6, ds_lecroy_516], dim='temp').unstack('temp')

# # How to normalize without having to unstack...avoids memory error
# mag_0_stack = [da_mws_nothing.loc[date].item().magnitude for date in ds_lecroy.coords['date'].values]
# mag_0_stack = xr.DataArray(mag_0_stack, coords={'temp':ds_lecroy.coords['temp']}, dims='temp')
# mag_0_stack = mag_0_stack.pint.quantify('V')

ds_lecroy = ds_lecroy.mws.calc_AS_abs(mag_0=da_mws_nothing).mws.calc_dAS().mws.calc_dPD()
ds_lecroy = ds_lecroy.xr_utils.stack_run()

#%%

ds_stats_lecroy = ds_lecroy.mws.calc_time_stats()[['dAS_abs_max', 'mag_pp', 'mag_fluct', 'SFR_abs', 'dpd1_max']]
ds_stats_lecroy = ds_stats_lecroy.mean('mnum', keep_attrs=True)

#%%

# Now that statistics are calculated we average over mnum and unstack. 
ds_lecroy = ds_lecroy.mean('mnum')



#%%k

ds_lecroy.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_lecroy.cdf'))


#%%
ds_cfd_beam = ppu.fileio.load_cfd_beam()

ds_cfd_beam = ds_cfd_beam.sel(offset=0).sel(kwt=1)

ds_cfd_beam = ds_cfd_beam.assign_coords(phi=ds_absem.coords['phi'].values)

#TODO: motor coordinates not exactly the same, why? 
ds_cfd_beam = ds_cfd_beam.sel(motor=ds_absem.coords['motor'].values, method='nearest')
ds_cfd_beam = ds_cfd_beam.assign_coords(motor=ds_absem.coords['motor'].values) 

da_cfd_beam = ds_cfd_beam['Yeq_K']
da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

da_cfd_beam

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_num_1

ds_fit = ds_absem.mean('mnum')
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam})

#%%

ds = ds_alpha_fit[['alpha_red', 'alpha_fit']]
ds['alpha'] = ds_fit['alpha']
ds = ds.rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha':'raw'})
# ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds
#%%

da_plot = ds.to_array('var').sel(phi=0.6, method='nearest').mean('run').dropna('motor',how='all')

#going to downselect 0.8 phi motor to same in 0.6. Mainly to have something that fits on a page. 
motor_vals_06 = da_plot.coords['motor'].values

da_plot.attrs['long_name'] = '$\\alpha$'
da_plot.plot(hue='var', col='mp', row='motor', figsize=(8,12))

plt.ylim(-0.1,1.1)

plt.savefig(pjoin(DIR_FIG_OUT, '536_pos_phi06_absem_fit.png'))
#%%

da_plot = ds.to_array('var').sel(phi=0.8, method='nearest').mean('run').dropna('motor',how='all')
da_plot = da_plot.sel(motor=motor_vals_06)

da_plot.attrs['long_name'] = '$\\alpha$'
da_plot.plot(hue='var', col='mp', row='motor', figsize=(8,12))

plt.ylim(-0.1,1.1)

plt.savefig(pjoin(DIR_FIG_OUT, '536_pos_phi08_absem_fit.png'))

#%%


#%%

nK_mw_horns = ds_p['nK_m3'].sel(mp='mw_horns').rename('nK_mw_horns').drop('mp')
nK_barrel = ds_p['nK_m3'].sel(mp='barrel').rename('nK_barrel').drop('mp')

ds = xr.merge([
    nK_mw_horns,
    nK_barrel,
    ds_stats_lecroy,
    ]).sortby('motor').dropna('run', how='all')


ds

ds.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_stats_pos.cdf'))
# %%

#%%
ds_alpha_fit.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_alpha_fit.cdf'))
# %%
