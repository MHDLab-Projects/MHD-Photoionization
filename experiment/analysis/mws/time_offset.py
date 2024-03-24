#%%[markdown]

# # MWS Raw data time offset

# Noticing 05-24 is offset in time, but not sure why.

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

#%%

folder = pjoin(REPO_DIR, 'experiment/data/munged/2023-05-24/Lecroy')

fp = pjoin(folder, 'ds_seedramp_2.cdf')

ds = xr.load_dataset(fp)

ds.coords['time'].attrs['units'] = 'second'
ds = ds.pint.quantify()
ds.coords['time'] = ds.coords['time'].pint.to('microsecond')

ds = ds.mws.calc_AS_rel()

# ds['AS'].dropna('acq_time','all')
# %%

ds['AS'].mean('acq_time').plot(marker='o')

plt.xlim(-2,2)

plt.axvline(0, color='k', linestyle='--')
# %%


fp = pjoin(REPO_DIR, 'experiment/data/proc_data/ds_lecroy.cdf')

ds = xr.load_dataset(fp)

ds = ds.mws.calc_AS_rel()

#%%


ds_day = ds.groupby(ds.acq_time.dt.day).mean('acq_time')

ds_day['AS'].plot(hue='day')

plt.xlim(-2,2)
# %%
ds.time



#%%
munged_dir = pjoin(REPO_DIR, 'experiment/data/munged')
dates = ['2023-04-07', '2023-05-12','2023-05-18', '2023-05-24']

for date in dates:
    fp = pjoin(munged_dir, date, 'ds_lecroy_time.cdf')
    ds = xr.load_dataset(fp)

    ds.coords['time'].attrs['units'] = 'second'
    ds = ds.pint.quantify()
    ds.coords['time'] = ds.coords['time'].pint.to('microsecond')

    ds = ds.mws.calc_AS_rel()

    ds['AS'].mean('acq_time').plot(label='date')

# plt.xlim(-2,2)
plt.xlim(-50,-48)


#%%

ds[['i', 'q']].coords['time']