#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mws_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mws_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')


# %%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3, pipe_fit_exp

dss_p = []
dss_fit_mws = []

ds_fit = ds_mws_53x
da_fit_lecroy = ds_fit.mws.calc_mag_phase_AS()['AS_abs']
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

dropna(g)
