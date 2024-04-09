
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdpy.analysis.mws.fitting import pipe_fit_mws_3, pipe_fit_exp

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mws_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mws_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')
# %%

ds_fit = ds_mws_53x
da_fit_lecroy = ds_fit['dAS_abs']

da_fit_lecroy = da_fit_lecroy/da_fit_lecroy.mws._pulse_max()
da_fit_lecroy = da_fit_lecroy.unstack('run')

#%%
import xyzpy

def run_model(da_mws, pipe_fit, **kwargs):

    pipe_fit.perform_fit_kwargs['fit_timewindow'] = slice(
        Quantity(kwargs['fit_window_min'], 'us'),Quantity(Quantity(kwargs['fit_window_max'], 'us'))
        )

    ds_fit_mws, ds_p_mws, ds_p_stderr_mws = pipe_fit_mws_3(da_mws)
    ds_p_mws['decay'] = 1/ds_p_mws['krm']

    # Interpolate ds_fit_mws to the time grid of da_fit_lecroy
    # There is coardsening as well as fit timewindow, so we need to interpolate back to the original time grid
    ds_fit_mws = ds_fit_mws.interp(time=da_fit_lecroy.time)

    return ds_fit_mws['AS_fit'].values, ds_p_mws['decay'].values, ds_p_mws['dne'].values

combos = {
    'fit_window_max': [10,15,20,25],
    'fit_window_min': [0.1, 0.5, 1, 2, 5]
}

r = xyzpy.Runner(
    run_model,
    var_names  =['AS_fit', 'decay', 'dne'],
    constants=dict(da_mws=da_fit_lecroy, pipe_fit=pipe_fit_mws_3),
    var_dims = {
        'AS_fit': da_fit_lecroy.dims, 
        'decay': da_fit_lecroy.isel(time=0).dims,
        'dne': da_fit_lecroy.isel(time=0).dims
        },
    var_coords=da_fit_lecroy.coords,

)

ds_mws = r.run_combos(combos)

ds_mws = xr.merge([ds_mws, da_fit_lecroy])

ds_mws = ds_mws.dropna('time',how='all')

ds_mws = ds_mws.xr_utils.stack_run()
#
# %%

#%%

ds_mws['decay'].sel(fit_window_min=1).mean('run').plot(hue='fit_window_max')

plt.ylim(1,20)

#%%

ds_mws['dne'].sel(fit_window_min=1).mean('run').plot(hue='fit_window_max')

#%%


ds_mws['decay'].sel(fit_window_max=15).mean('run').plot(hue='fit_window_min')

plt.ylim(1,20)

#%%

ds_mws['dne'].sel(fit_window_max=15).mean('run').plot(hue='fit_window_min')
#%%

da_plot = ds_mws['AS_fit'].dropna('time','all').mean('run').sel(fit_window_max=15)


g = da_plot.plot(row='kwt',hue='fit_window_min')

plt.yscale('log')
# Iterate through the axes (i.e., the individual plots in the facet grid)
for i, ax in enumerate(g.axes.flat):
    # Get the kwt value for this plot
    kwt = da_plot.kwt.values[i]

    # Select the corresponding data from ds_mws['AS_abs']
    da_abs = ds_mws['dAS_abs'].sel(kwt=kwt).mean('run')

    # Add the data to the plot
    ax.plot(da_abs.time, da_abs, label='AS_abs', color='gray')

    # Set the yscale to log
    ax.set_yscale('log')

plt.xlim(-1, 20)
plt.ylim(1e-2, 1.1)

# %%
