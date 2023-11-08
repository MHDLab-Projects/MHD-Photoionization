


#%%
from mhdpy.analysis.standard_import import *

data_folder = mhdpy.io.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-04-07')

figure_out_dir = pjoin(DIR_DATA_OUT, '1d_auto')

dsst = mhdpy.io.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%

tcs = [
    '536_pos_1',
    '536_power_1',
    '53x_2',
    '53x_1'
]

from mhdpy.mws_utils import calc_mag_phase_AS

dss = {}
dss_std = {}

for tc in tcs:

    fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

    ds_in = xr.load_dataset(fp_in)

    ds = calc_mag_phase_AS(ds_in)


    tc_dim = [dim for dim in ds.dims if dim != 'time'][0]
    mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

    ds_mean = ds.mean('mnum', keep_attrs=True)
    ds_std = ds.std('mnum', keep_attrs=True)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    if not os.path.exists(tc_fig_dir): os.makedirs(tc_fig_dir)

    for fn in os.listdir(tc_fig_dir):
        os.remove(os.path.join(tc_fig_dir, fn))

    plt.figure()

    mnum_counts.plot(marker='o')
    plt.ylabel('Acquisitions')


    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'acquisition_counts.png'))


    dss[tc] = ds_mean
    dss_std[tc] = ds_std

# %%

figsize = (8, 11)
for tc in dss:
    ds = dss[tc]

    tc_dim = [dim for dim in ds.dims if dim != 'time'][0]

    plot_names = np.array([ds.data_vars]).T
    fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

    for var in ds.data_vars:
        ax = axes[var]
        ds[var].plot(ax=ax, hue=tc_dim)

        if var != 'i': ax.get_legend().remove()
        if var != 'i': ax.set_title('')

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'raw_signals_hue.png'))

# %%

figsize = (8, 11)
for tc in dss:
    ds = dss[tc]

    plot_names = np.array([ds.data_vars]).T
    fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

    for var in ds.data_vars:
        ax = axes[var]
        ds[var].plot(ax=ax)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'raw_signals_color.png'))


#%%

for tc in dss:

    plt.figure()
    ds = dss[tc]

    tc_dim = [dim for dim in ds.dims if dim != 'time'][0]

    ds['AS'].plot(hue=tc_dim)

    plt.yscale('log')

    plt.xlim(-1,30)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'AS_log.png'))

#%%

from mhdpy.plot.common import xr_errorbar

for tc in dss:

    ds = dss[tc]
    ds_std = dss_std[tc]

    tc_dim = [dim for dim in ds.dims if dim != 'time'][0]

    plt.figure()

    tc_vals = ds.coords[tc_dim]
    fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)
    # ds['AS'].plot(hue=tc_dim)
    # dss_std[tc]['AS'].plot(hue=tc_dim, linestyle='--')

    ds = ds.drop_vars([var for var in ds.coords if var not in ds.dims])

    # xr_errorbar(ds['AS'], dss_std[tc]['AS'], huedim=tc_dim)

    for i, c in enumerate(ds.coords[tc_dim]):
        ds_sel_mean = ds.sel({tc_dim: c})['AS']
        xs = ds_sel_mean.coords['time'].values

        vals_mean = ds_sel_mean.values
        vals_std = ds_std.sel({tc_dim: c})['AS'].values

        plt.sca(axes[i])
        plt.plot(xs, vals_mean)

        plt.fill_between(x=xs,
                 y1=vals_mean - vals_std,
                 y2=vals_mean + vals_std,
                 alpha=0.5
                 )


    # plt.yscale('log')

    # plt.xlim(-1,30)
    # plt.ylim(-1e-2, 1e-2)


    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'AS_lin.png'))

#%%

import seaborn as sns
sns.set_theme(style="darkgrid")

# Load an example dataset with long-form data
fmri = sns.load_dataset("fmri")

# Plot the responses for different events and regions
sns.lineplot(x="timepoint", y="signal",
             hue="region", style="event",
             data=fmri)


# %%

max_window = slice(-1,1)
pre_pulse = slice(-50,-1)

for tc in dss:
    ds = dss[tc]

    plot_names = [['pulse_max'], ['pre_pulse_avg'], ['pre_pulse_std']]
    fig, axes = plt.subplot_mosaic(plot_names, figsize=figsize, sharex=True)

    pulse_max = ds['AS'].sel(time=max_window).max('time')
    pulse_max.name = 'pulse_max'
    pre_pulse_avg =ds['AS'].sel(time=pre_pulse).mean('time') 
    pre_pulse_avg.name = 'pre_pulse_avg'
    pre_pulse_std =ds['AS'].sel(time=pre_pulse).std('time') 
    pre_pulse_std.name = 'pre_pulse_std'

    ds_out = xr.merge([
        pulse_max,
        pre_pulse_avg,
        pre_pulse_std
    ])

    for var in ds_out.data_vars:
        ax = axes[var]
        ds_out[var].plot(ax=ax, marker='o')
    
    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'pulse_stats.png'))


# %%

from mhdpy.mws_utils.fitting import fit_fn 
from mhdpy.analysis.xr import fit_da_lmfit
from lmfit import Model

das_fit_input = {}
dss_fits = {}
dss_p = {}
dss_p_stderr ={}

for tc in dss:
    print(tc)
    ds = dss[tc]

    da_fit = ds['AS'].copy()
    pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

    tc_dim = [dim for dim in da_fit.dims if dim != 'time'][0]

    da_fit = da_fit.where(pulse_max > 5e-4) # Targeting low power...
    da_fit = da_fit.dropna(tc_dim,'all')
    pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

    da_fit = da_fit/pulse_max

    das_fit_input[tc] = da_fit

    mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0))

    params = mod.make_params()
    params['dne'].value = 1e14
    params['ne0'].value = 1e12
    params['dne'].min = 0
    params['ne0'].min = 0


    da_fit_region = da_fit.sel(time=slice(0,15))
    x_eval = da_fit.sel(time=slice(0,15)).coords['time'].values

    # da_fit_coarse = da_fit_region.coarsen(time=10).mean()
    fits, ds_p, ds_p_stderr = fit_da_lmfit(da_fit_region, mod, params, 'time', x_eval)


    dss_fits[tc] = fits
    dss_p[tc] = ds_p
    dss_p_stderr[tc] = ds_p_stderr

    # break
# %%
for tc in dss_p:

    plt.figure()
    ds_p = dss_p[tc]

    tc_dim = [dim for dim in ds_p.dims if dim != 'time'][0]

    ds_p.to_array('var').plot(row='var', yscale = 'log', sharey=False, marker='o')


    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'ds_p.png'))

#%%

for tc in dss_fits:
    plt.figure()

    fits = dss_fits[tc]
    da_fit = das_fit_input[tc]


    tc_dim = [dim for dim in fits.dims if dim != 'time'][0]
    tc_vals = da_fit.coords[tc_dim]

    fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)

    for i, tc_val in enumerate(tc_vals):

        ax = axes[i]

        da_fit.sel({tc_dim: tc_val}).plot(ax=ax)
        fits.sel({tc_dim: tc_val}).plot(ax=ax, linestyle='--', color='gray')

        ax.set_yscale('log')
        # plt.ylim(-0.1,1.1)
        ax.set_xlim(-1,30)

        if i != len(tc_vals) -1: ax.set_xlabel('')

        # if i != 0:
        #     ax.set_title('')

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'fits_all.png'))
