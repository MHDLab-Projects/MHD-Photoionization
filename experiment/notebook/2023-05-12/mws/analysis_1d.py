


#%%
from mhdpy.analysis.standard_import import *

datestr = '2023-05-12'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

figure_out_dir = pjoin(DIR_DATA_OUT, '1d_auto')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%

tcs = [
    '536_pos',
    '536_power',
    '53x',
]

from mhdpy.mws_utils import calc_mag_phase_AS


dss = {}

for tc in tcs:

    fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

    ds_in = xr.load_dataset(fp_in)
    ds_in = ds_in.sel(date=datestr).sel(run_num=1)

    ds = calc_mag_phase_AS(ds_in).drop('mag_pp')


    tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
    mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

    ds = ds.mean('mnum', keep_attrs=True)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    if not os.path.exists(tc_fig_dir): os.makedirs(tc_fig_dir)

    for fn in os.listdir(tc_fig_dir):
        os.remove(os.path.join(tc_fig_dir, fn))

    plt.figure()

    mnum_counts.plot(marker='o')
    plt.ylabel('Acquisitions')


    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'acquisition_counts.png'))


    dss[tc] = ds

# %%

figsize = (8, 11)
for tc in dss:
    ds = dss[tc]

    tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]

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

    tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]

    ds['AS'].plot(hue=tc_dim)

    plt.yscale('log')

    plt.xlim(-1,30)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'AS_log.png'))


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

    da_fit = da_fit.where(pulse_max > 5e-5) # Targeting low power...
    da_fit = da_fit.dropna(tc_dim, how='all')
    pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

    da_fit = da_fit/pulse_max

    das_fit_input[tc] = da_fit

    mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0))

    params = mod.make_params()
    params['dne'].value = 1e14
    params['ne0'].value = 1e12
    params['dne'].min = 1e11
    params['ne0'].min = 1e11


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

    g = ds_p.to_array('var').plot(row='var', yscale = 'log', sharey=False, marker='o')

    # g.axes[1].set_ylim()


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
