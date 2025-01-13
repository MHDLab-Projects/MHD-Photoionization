"""
Auto analysis script for mwt noise (taken from 2023-04-07/analysis_noise.py)
Keeping separate from mwt_1d.py to use medthods in absem_1d.py that avoid for loops, eventually want to use these methods in mwt_1d.py
"""
#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mwt_1d import tc_dict

import argparse
import shutil

tc = '536_pos'


# # Create the parser
# parser = argparse.ArgumentParser(description='Process some data.')

# # Add a positional argument for tc
# parser.add_argument('tc', type=str, help='a required tc argument')

# # Parse the arguments
# args = parser.parse_args()

# # Now you can use args.tc as your tc variable
# tc = args.tc

tc_dim_dict = {
    '53x': 'kwt',
    '536_pos': 'motor',
    '536_power': 'power'
}

tc_dim = tc_dim_dict[tc]

figure_out_dir = pjoin(DIR_DATA_OUT, 'mwt_1d_noise', tc)
if os.path.exists(figure_out_dir):
    shutil.rmtree(figure_out_dir)
os.makedirs(figure_out_dir)


mag_0 = ppu.fileio.load_mwt_T0()

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy = ds_lecroy.mwt.calc_AS_abs(mag_0=mag_0) # calculate before or after mnum mean?
ds_lecroy = ds_lecroy.xr_utils.stack_run()
ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected

ds_lecroy = ds_lecroy.drop_vars('acq_time').drop_vars('mag_0')
ds_lecroy = ds_lecroy.dropna('run', how='all')

# %%

ds_lecroy


#%%

def groupby_run_processor(ds_plot, mnum_limit=None):
    #TODO: assert that there is only one run coord, which should be true as we are grouping by it
    run_name = ds_plot.coords['run'].item()
    run_name = "{}_{}".format(run_name[0], run_name[1])

    ds_plot = ds_plot.dropna(tc_dim, how='all')

    if 'mnum' in ds_plot.dims:
        ds_plot = ds_plot.groupby(tc_dim).apply(lambda x: x.xr_utils.assign_mnum('mnum'))
        if mnum_limit is not None:
            ds_plot = ds_plot.isel(mnum=range(mnum_limit))

    return ds_plot, run_name

def gen_path_fig_out(run_name, plot_name, tc):
    fp_out = pjoin(figure_out_dir, plot_name, '{}.png'.format(run_name))

    # Make directory is not exists
    os.makedirs(os.path.dirname(fp_out), exist_ok=True)

    return fp_out


#%%

def plot_1(ds_plot, run_name, figsize=(8,11), savefig=True, plot_name='p1'):
    n_plots = len(ds_plot.coords[tc_dim])

    fig, axes = plt.subplots(n_plots, figsize=(4,20), sharex=True)

    for i, ax in enumerate(axes):

        da = ds_plot['AS_abs'].isel({tc_dim: i}).plot(ax=ax)

    if savefig:
        plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return ds_plot


ds_lecroy.groupby('run').apply(lambda x: plot_1(*groupby_run_processor(x, mnum_limit=150)))
# %%
new_dss = []

for mnum_max in range(150):
    mnum_range = slice(0,mnum_max)

    ds_sel=ds_lecroy.sel(mnum=mnum_range)

    ds_sel = ds_sel.mean('mnum')

    ds_sel = ds_sel.assign_coords(mnum_max=mnum_max)

    new_dss.append(ds_sel)


ds_lecroy_mm = xr.concat(new_dss, dim='mnum_max')

#%%

def plot_2(ds_plot, run_name, figsize=(8,11), savefig=True, plot_name='p2'):
    n_plots = len(ds_plot.coords[tc_dim])

    fig, axes = plt.subplots(n_plots, figsize=(4,20), sharex=True)

    for i, ax in enumerate(axes):

        ds_plot['AS_abs'].isel({tc_dim: i}).plot(ax=ax)

        if i < len(axes)-1:
            ax.set_xlabel('')

    if savefig:
        plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return ds_plot

ds_lecroy_mm.groupby('run').apply(lambda x: plot_2(*groupby_run_processor(x), plot_name='mm_all'))
ds_lecroy_mm.sel(time=slice(-50,-1)).groupby('run').apply(lambda x: plot_2(*groupby_run_processor(x), plot_name='mm_pp'))

#%%

def plot_3(ds_plot, run_name, figsize=(8,11), savefig=True, plot_name='p3'):
    plt.figure()
    ds_plot['AS_abs'].sel(time=slice(-50,-1)).std('time').plot(hue=tc_dim)

    plt.yscale('log')
    plt.ylabel("AS Pre Pulse Std Dev.")

    if savefig:
        plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return ds_plot


ds_lecroy_mm.groupby('run').apply(lambda x: plot_3(*groupby_run_processor(x), plot_name='mm_pp_std'))

#%%

da = ds_lecroy['mag'].pint.dequantify()

da

#%%

from xhistogram.xarray import histogram

def hist_plot_1(da_mag, run_name, figsize=(8,11), savefig=True, plot_name='hist1'):

    # bins = np.linspace(-.2,.2)
    bins = np.linspace(da_mag.min().item(),da_mag.max().item())
    h = histogram(da_mag, bins=bins, dim=['mnum'])

    h.plot(col=tc_dim, col_wrap=3)

    if savefig:
        # plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return da_mag

da.groupby('run').apply(lambda x: hist_plot_1(*groupby_run_processor(x), plot_name='hist1'))

#%%

# def hist_plot_2(da_mag, run_name, figsize=(8,11), savefig=True, plot_name='hist2'):


#     bins = np.linspace(da_mag.min(),da_mag.max())
#     h = histogram(da_mag, bins=bins, dim=['mnum'])

#     tc_vals = h.coords[tc_dim]
#     fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)

#     for i, c in enumerate(h.coords[tc_dim]):
#         h_sel = h.sel({tc_dim: c})
#         h_sel.plot(ax=axes[i])

#     if savefig:
#         # plt.tight_layout()
#         fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
#         plt.savefig(fp_fig_out)

#     return da_mag


# da.groupby('run').apply(lambda x: hist_plot_2(*groupby_run_processor(x), plot_name='hist2'))



#%%

from matplotlib.ticker import FuncFormatter

# Define a function to format the tick values
def format_tick(x, pos):
    return f'{x:.3f}'

formatter = FuncFormatter(format_tick)

def hist_plot_2(da_mag, run_name, figsize=(8,11), savefig=True, plot_name='hist2'):

    col_wrap = 3

    # da_mag = da.sel(run=('2023-05-18', 1)).dropna(tc_dim, how='all')

    bins = np.linspace(da_mag.min().item(),da_mag.max().item())
    h = histogram(da_mag, bins=bins, dim=['mnum'])
    h.attrs['long_name'] = 'Bin Count'

    # Calculate the number of rows needed
    nrows = np.ceil(len(h[tc_dim]) / col_wrap).astype(int)

    # Create a 2D array of Axes objects
    fig, axes = plt.subplots(nrows=nrows, ncols=col_wrap, figsize=(3*col_wrap, 3*nrows))

    # Flatten the array of Axes objects for easy iteration
    axes_flat = axes.flatten()

    # Loop over the Axes objects and plot the data
    for ax, (label, da_plot) in zip(axes_flat, h.groupby(tc_dim)):
        da_plot.plot(ax=ax)
        ax.set_title(label)
        ax.xaxis.set_major_formatter(formatter)
        ax.tick_params(axis='x', labelcolor='black', color='white', rotation=45)  # Set the tick labels to black and ticks to white

    # remove x labels from each row except the bottom
    # remove y labels from each column except the left

    if nrows > 1:
        for ax in axes[:, 1:].flatten():
            ax.set_ylabel('')

    
    for ax in axes_flat[:-col_wrap]:
        ax.set_xlabel('')

    plt.tight_layout()

    if savefig:
        # plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return da_mag



da.groupby('run').apply(lambda x: hist_plot_2(*groupby_run_processor(x), plot_name='hist2'))

#%%


def hist_plot_3(da_mag, run_name, figsize=(8,11), savefig=True, plot_name='hist3'):

    # da_mag = da.sel(run=('2023-05-18', 1)).dropna(tc_dim, how='all')
    da_mag = da_mag.sel(time=slice(-50,1))

    bins = np.linspace(da_mag.min().item(),da_mag.max().item())
    h = histogram(da_mag, bins=bins, dim=['mnum'])

    h.attrs['long_name'] = 'Bin Count'

    h.mean('time').plot(col=tc_dim, col_wrap=3, marker='o')

    if savefig:
        # plt.tight_layout()
        fp_fig_out = gen_path_fig_out(run_name, plot_name, tc)
        plt.savefig(fp_fig_out)

    return da_mag


da.groupby('run').apply(lambda x: hist_plot_3(*groupby_run_processor(x), plot_name='hist3'))

#%%



da_mag = da.sel(time=slice(-50,1))

bins = np.linspace(da_mag.min().item(),da_mag.max().item())
h = histogram(da_mag, bins=bins, dim=['mnum'])

h.attrs['long_name'] = 'Bin Count'
h = h.mean('time', keep_attrs=True)

h = h.unstack('run').xr_utils.stack_run()

h.plot(col=tc_dim, hue='run_plot', x='mag_bin', col_wrap=3, marker='o')


plt.yscale('log')

plt.ylim(1e-3,1e3)

gen_path_fig_out('all', 'hist4', tc)
plt.savefig(pjoin(figure_out_dir, 'hist4', 'hist4_{}.png'.format(tc)))

#%%

