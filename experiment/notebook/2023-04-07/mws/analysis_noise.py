"""
Analysis of the pre pulse noise vs motor position
"""

#%%
from mhdpy.analysis.standard_import import *

datestr = '2023-04-07'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')


dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%


from mhdpy.analysis import mws

tc = '536_pos'

fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)
ds_in = ds_in.sel(date=datestr).sel(run_num=1)

ds = ds_in.mws.calc_mag_phase_AS().drop('mag_pp')


tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

# ds = ds.mean('mnum', keep_attrs=True)


ds = ds.drop('power').drop('kwt')


# %%
ds['AS'].mean('mnum').plot(hue='motor')

#%%

mnum_counts

#%%

ds_2 = ds.sel(mnum=slice(0,150))
# ds_2 = ds_2.sel(time=slice(-50,-1))

ds_2

#%%

n_plots = len(ds_2.coords['motor'])

fig, axes = plt.subplots(n_plots, figsize=(4,20), sharex=True)

for i, ax in enumerate(axes):

    ds_2['AS'].isel(motor=i).plot(ax=ax)

# %%

new_dss = []

for mnum_max in range(150):
    mnum_range = slice(0,mnum_max)

    ds_sel=ds_2.sel(mnum=mnum_range)

    ds_sel = ds_sel.mean('mnum')

    ds_sel = ds_sel.assign_coords(mnum_max=mnum_max)

    new_dss.append(ds_sel)


# %%
ds_3 = xr.concat(new_dss, 'mnum_max')


# %%

n_plots = len(ds_3.coords['motor'])

fig, axes = plt.subplots(n_plots, figsize=(4,20), sharex=True)

for i, ax in enumerate(axes):

    ds_3['AS'].isel(motor=i).plot(ax=ax)

    if i < len(axes)-1:
        ax.set_xlabel('')
#%%


n_plots = len(ds_3.coords['motor'])

fig, axes = plt.subplots(n_plots, figsize=(4,20), sharex=True)

for i, ax in enumerate(axes):

    ds_3['AS'].sel(time=slice(-50,-1)).isel(motor=i).plot(ax=ax)

    if i < len(axes)-1:
        ax.set_xlabel('')

#%%


ds_3['AS'].sel(time=slice(-50,-1)).std('time').plot(hue='motor')

plt.yscale('log')
plt.ylabel("AS Pre Pulse Std Dev.")


# %%

ds


ds_mean = ds.mean('mnum')
ds_std = ds.std('mnum')

tc_dim = [dim for dim in ds_mean.dims if dim != 'time'][0]

plt.figure()

tc_vals = ds.coords[tc_dim]
fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)

for i, c in enumerate(ds.coords[tc_dim]):
    ds_sel_mean = ds_mean.sel({tc_dim: c})['AS']
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

#%%

da = ds['mag']
da.sel(mnum=0).plot()

#%%

da.plot.hist()


#%%

from xhistogram.xarray import histogram

# bins = np.linspace(-.2,.2)
bins = np.linspace(da.min(),da.max())
h = histogram(da, bins=bins, dim=['mnum'])

h.plot(row='motor')
# %%

tc_dim = 'motor'
tc_vals = h.coords[tc_dim]
fig, axes = plt.subplots(len(tc_vals), figsize=(5, 3*len(tc_vals)), sharex=True)



for i, c in enumerate(h.coords[tc_dim]):
    h_sel = h.sel({tc_dim: c})
    h_sel.plot(ax=axes[i])