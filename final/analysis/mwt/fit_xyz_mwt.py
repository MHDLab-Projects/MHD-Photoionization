
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdpy.analysis.mwt.fitting import pipe_fit_mwt_3, pipe_fit_exp

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mwt_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mwt_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')
# %%

ds_fit = ds_mwt_53x
da_fit_lecroy = ds_fit['dAS_abs']

da_fit_lecroy = da_fit_lecroy/da_fit_lecroy.mwt._pulse_max()
da_fit_lecroy = da_fit_lecroy.unstack('run')

#%%
import xyzpy

def run_model(da_mwt, pipe_fit, **kwargs):

    pipe_fit.perform_fit_kwargs['fit_timewindow'] = slice(
        Quantity(kwargs['fit_window_min'], 'us'),Quantity(Quantity(kwargs['fit_window_max'], 'us'))
        )

    ds_fit_mwt, ds_p_mwt, ds_p_stderr_mwt = pipe_fit_mwt_3(da_mwt)
    ds_p_mwt['decay'] = 1/ds_p_mwt['krm']

    # Interpolate ds_fit_mwt to the time grid of da_fit_lecroy
    # There is coardsening as well as fit timewindow, so we need to interpolate back to the original time grid
    ds_fit_mwt = ds_fit_mwt.interp(time=da_fit_lecroy.time)

    return ds_fit_mwt['AS_fit'].values, ds_p_mwt['decay'].values, ds_p_mwt['dne'].values

combos = {
    'fit_window_max': [10,15,20,25],
    'fit_window_min': [0.1, 0.5, 1, 2, 5]
}

r = xyzpy.Runner(
    run_model,
    var_names  =['AS_fit', 'decay', 'dne'],
    constants=dict(da_mwt=da_fit_lecroy, pipe_fit=pipe_fit_mwt_3),
    var_dims = {
        'AS_fit': da_fit_lecroy.dims, 
        'decay': da_fit_lecroy.isel(time=0).dims,
        'dne': da_fit_lecroy.isel(time=0).dims
        },
    var_coords=da_fit_lecroy.coords,

)

ds_mwt = r.run_combos(combos)

ds_mwt = xr.merge([ds_mwt, da_fit_lecroy])

ds_mwt = ds_mwt.dropna('time',how='all')

ds_mwt = ds_mwt.xr_utils.stack_run()
#
# %%

#%%

ds_mwt['decay'].sel(fit_window_min=1).mean('run').plot(hue='fit_window_max')

plt.ylim(1,20)

#%%

ds_mwt['dne'].sel(fit_window_min=1).mean('run').plot(hue='fit_window_max')

#%%


ds_mwt['decay'].sel(fit_window_max=15).mean('run').plot(hue='fit_window_min')

plt.ylim(1,20)

#%%

ds_mwt['dne'].sel(fit_window_max=15).mean('run').plot(hue='fit_window_min')
#%%

da_plot = ds_mwt['AS_fit'].dropna('time','all').mean('run').sel(fit_window_max=15)


g = da_plot.plot(row='kwt',hue='fit_window_min')

plt.yscale('log')
# Iterate through the axes (i.e., the individual plots in the facet grid)
for i, ax in enumerate(g.axes.flat):
    # Get the kwt value for this plot
    kwt = da_plot.kwt.values[i]

    # Select the corresponding data from ds_mwt['AS_abs']
    da_abs = ds_mwt['dAS_abs'].sel(kwt=kwt).mean('run')

    # Add the data to the plot
    ax.plot(da_abs.time, da_abs, label='AS_abs', color='gray')

    # Set the yscale to log
    ax.set_yscale('log')

plt.xlim(-1, 20)
plt.ylim(1e-2, 1.1)

# %%
