#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mws_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mws_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')

#%%


# %%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3, pipe_fit_exp

pipe_fit_mws_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(1, 'us'),Quantity(Quantity(20, 'us')))

dss_p = []
dss_fit_mws = []

ds_fit = ds_mws_53x
da_fit_lecroy = ds_fit.mws.calc_mag_phase_AS()['AS_abs']

da_fit_lecroy = da_fit_lecroy/da_fit_lecroy.mws._pulse_max()

#%%

da_fit_lecroy.plot(hue='run_plot', x='time', row='kwt')

plt.yscale('log')
plt.ylim(1e-2, 5.1)
plt.xlim(-1,1)

#%%

# ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_exp(da_fit_lecroy, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))
ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_mws_3(da_fit_lecroy)
ds_p_mws['decay'] = 1/ds_p_mws['krm']

dss_p.append(ds_p_mws.assign_coords(method='dnedt'))
dss_fit_mws.append(ds_fit_mws.assign_coords(method='dnedt'))

#%%

ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_exp(da_fit_lecroy)

dss_p.append(ds_p_mws.assign_coords(method='exp'))
dss_fit_mws.append(ds_fit_mws.assign_coords(method='exp'))

#%%

ds_p = xr.concat(dss_p, 'method')
ds_fit_mws = xr.concat(dss_fit_mws, 'method')

ds_p['decay'].mean('run').plot(hue='method')

plt.yscale('log')

plt.ylim(1,20)
# %%

g= ds_fit_mws.mean('run').to_array('var').plot(row='kwt', col='method', hue='var', figsize=(5,10))

plt.yscale('log')
plt.ylim(1e-4, 1.1)

plt.xlim(-1, )

dropna(g)

# %%

#TODO: why can't I use sel here?
ds_dnedt = ds_fit_mws.isel(method=0).isel(run=[-1])

g= ds_dnedt.mean('run').to_array('var').plot(row='kwt', hue='var', figsize=(5,10))
# g= ds_fit_mws.isel(method='dnedt').mean('run').to_array('var').plot(row='kwt', hue='var', figsize=(5,10)


g

plt.yscale('log')
plt.ylim(1e-2, 5.1)
plt.xlim(-1,5)

dropna(g)

#%%

g = ds_dnedt['AS_all'].plot(hue='run_plot', row='kwt', x='time')

dropna(g)

plt.yscale('log')
plt.ylim(1e-2, 5.1)
plt.xlim(-1,5)