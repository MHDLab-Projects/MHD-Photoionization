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

#%%

ds_lecroy_536_pos

#%%

# add 536 photodiode data from 5x6_pos run. Do not have photodiode data for 536_pos run

tc = '5x6_pos'

ds_lecroy_5x6 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy_5x6 = ds_lecroy_5x6.sel(phi=[0.79], method='nearest') # Cannot add 0.65 phi, because it conflicts with 516_pos. Need to change run_num for this measurement?

ds_lecroy_5x6 = ds_lecroy_5x6[['pd1','pd2']]

ds_lecroy_5x6 = ds_lecroy_5x6.stack(temp = stack_dims)

#%%

tc = '516_pos'

ds_absem_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem_516 = ds_absem_516.expand_dims('phi') # allows merging with 5x6_pos
ds_absem_516 = ds_absem_516.stack(temp = stack_dims)

ds_lecroy_516 = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy_516 = ds_lecroy_516.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy_516 = ds_lecroy_516.expand_dims('phi') # allows merging with 5x6_pos
ds_lecroy_516 = ds_lecroy_516.stack(temp = stack_dims)  

#%%

ds_absem = xr.concat([ds_absem_536_pos, ds_absem_516], dim='temp').unstack('temp')
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.absem.calc_alpha()

da_mws_nothing  = ppu.fileio.load_mws_T0()

ds_lecroy = xr.concat([ds_lecroy_536_pos, ds_lecroy_5x6, ds_lecroy_516], dim='temp').unstack('temp')
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS(mag_0=da_mws_nothing.pint.dequantify())
ds_lecroy = ds_lecroy.xr_utils.stack_run()


#%%

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
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, da_nK_profile=da_cfd_beam)


ds = ds_alpha_fit[['alpha_red', 'alpha_fit']]
ds = ds.rename({'alpha_red': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%

nK_mw_horns = ds_p['nK_m3'].sel(mp='mw_horns').rename('nK_mw_horns').drop('mp')
nK_barrel = ds_p['nK_m3'].sel(mp='barrel').rename('nK_barrel').drop('mp')

mws_max = ds_lecroy['AS_abs'].mean('mnum').max('time').rename('AS_max')
mws_pp = ds_lecroy['mag_pp'].mean('mnum').rename('mag_pp')
mws_pp_std = ds_lecroy['mag'].sel(time=slice(-50,-1)).std('time').mean('mnum').rename('mag_pp_std')
mws_max_std_ratio = mws_max/mws_pp_std
mws_max_std_ratio = mws_max_std_ratio.rename('AS_max_std_ratio')

delta_pd1 = ds_lecroy['pd1'] - ds_lecroy['pd1'].sel(time=slice(-1,0)).mean('time')
delta_pd1 = delta_pd1.dropna('run', how='all')
delta_pd1 = delta_pd1.mean('mnum').max('time')
# delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')

ds = xr.merge([
    mws_max, 
    mws_pp,
    mws_pp_std,
    mws_max_std_ratio,
    nK_mw_horns,
    nK_barrel,
    delta_pd1.rename('delta_pd1') 
    ]).sortby('motor').dropna('run', how='all')

ds['AS_max'].attrs = dict(long_name='AS Max')
ds['mag_pp'].attrs = dict(long_name='Mag. Pre Pulse (PP)', units='V')
ds['mag_pp_std'].attrs = dict(long_name='Pre Pulse (PP) Std Dev.', units='V')
ds['AS_max_std_ratio'].attrs = dict(long_name='AS Max / PP Std Dev.', units='1/V')

ds

ds.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_stats_pos.cdf'))
# %%

#%%
ds_alpha_fit['alpha'] = ds_fit['alpha']
ds_alpha_fit.pint.dequantify().unstack('run').to_netcdf(pjoin(DIR_DATA_OUT, 'ds_pos_alpha_fit.cdf'))