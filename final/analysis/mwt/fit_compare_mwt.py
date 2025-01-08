#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mwt_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mwt_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')

#%%


from mhdpy.analysis.mwt.fitting import MWTFitPipeline3, MWTFitPipelineExp

def perform_fit_sequence(ds_fit):
    pipe_fit_exp = MWTFitPipelineExp()
    pipe_fit_mwt_3 = MWTFitPipeline3()

    pipe_fit_mwt_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(1, 'us'),Quantity(Quantity(20, 'us')))

    dss_p = []
    dss_fit_mwt = []

    ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_exp(ds_fit)
    dss_p.append(ds_p_mwt.assign_coords(fit_method='exp'))
    dss_fit_mwt.append(ds_fit_mwt.assign_coords(fit_method='exp'))

    ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_mwt_3(ds_fit)
    ds_p_mwt['decay'] = 1/ds_p_mwt['krm']
    dss_p.append(ds_p_mwt.assign_coords(fit_method='dnedt'))
    dss_fit_mwt.append(ds_fit_mwt.assign_coords(fit_method='dnedt'))

    ds_p = xr.concat(dss_p, 'fit_method')
    ds_fit_mwt = xr.concat(dss_fit_mwt, 'fit_method')

    return ds_fit_mwt, ds_p

# %%

ds_fit = ds_mwt_53x
da_fit_lecroy = ds_fit['dAS_abs']
da_fit_lecroy = da_fit_lecroy/da_fit_lecroy.mwt._pulse_max()

#%%

da_fit_lecroy

#%%

da_fit_lecroy.plot(hue='run_plot', x='time', row='kwt')

plt.yscale('log')
plt.ylim(1e-2, 5.1)
plt.xlim(-1,1)

#%%

# ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_exp(da_fit_lecroy, fit_method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))
# ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_mwt_3(da_fit_lecroy)

ds_fit_mwt, ds_p = perform_fit_sequence(da_fit_lecroy)

ds_p['decay'].mean('run').plot(hue='fit_method')

plt.yscale('log')

plt.ylim(1,20)
# %%

g= ds_fit_mwt.mean('run').to_array('var').plot(row='kwt', col='fit_method', hue='var', figsize=(5,10))

plt.yscale('log')
plt.ylim(1e-4, 1.1)

plt.xlim(-1, )

dropna(g)

# %%

#TODO: why can't I use sel here?
ds_dnedt = ds_fit_mwt.sel(fit_method='dnedt').isel(run=[-1])

g= ds_dnedt.mean('run').to_array('var').plot(row='kwt', hue='var', figsize=(5,10))
# g= ds_fit_mwt.isel(fit_method='dnedt').mean('run').to_array('var').plot(row='kwt', hue='var', figsize=(5,10)


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
# %%
