# %%
from mhdlab.analysis.standard_import import *

from mhdlab.fileio import TFxr
from mhdlab.fileio.path import gen_path_date
from mhdlab.fileio.spectral import load_aas

from mhdlab.xr_utils import interp_ds_to_var
from mhdlab.analysis.aas import calc_alpha_simple


import json 
datestr = '2023-05-24'

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

start, stop = map(pd.Timestamp, settings['calib_timewindow'])
calib_timewindow = slice(start, stop)

has_multiplexer = settings['has_multiplexer']

#%%

#TODO: make consistent with main aas_setup


#%%

# Start processing 

data_folder = os.path.join(REPO_DIR, 'experiment','data', 'munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

fp = os.path.join(data_folder,'Munged','Spectral' ,'aas.tdms')
ds_aas = load_aas(fp, convert_to_PT=False)

# Needed for scipp
ds_aas = ds_aas.assign_coords(led = ('time', ds_aas.coords['led'].astype(np.int64).values))


ds

#%%

ds2 = ds.reset_index('switch_num').reset_coords('switch_num').set_index(acq_group='time',append=True).unstack('acq_group')

tw = slice(Timestamp('2023-05-24 20:00:00'), Timestamp('2023-05-24 20:01:00'), None)

ds2 = ds2.sel(time=tw).mean('wavelength')

ds2 = ds2.sel(led=0, mp='barrel').dropna('time')

ds2

#%%


ds2['counts'].plot()
#%%

# ds2.groupby('switch_num').mean()

ds2.swa


#%%

ds2 = reduce_acq_group(ds)

# %%

import scipp as sc

ds_mean = ds2['counts_mean']
ds_mean.coords['wavelength'].attrs['units'] = 'nm'

sds = sc.compat.from_xarray(ds_mean)

sds
#%%

sds.variances = ds2['counts_std']

#%%

import matplotlib.pyplot as plt

sds['time', 0].plot()


# plt.ylim(0,100)

#%%

# sds['time', ]


sds['wavelength', 760*sc.Unit('nm'): 780*sc.Unit('nm')].plot()

#%%

ts = sds.coords['time'][2900:3000].values

sds_sel = sds['time',ts[0]*sc.Unit('ns'):ts[-1]*sc.Unit('ns')]

sds_sel['wavelength', 20].plot()

#%%

# calib_timewindow


# sds.coords['time']

start = np.datetime64(calib_timewindow.start)
stop = np.datetime64(calib_timewindow.stop)

start = sc.scalar(start, unit='ns')
stop = sc.scalar(stop, unit='ns')

sds_sel = sds['time', start:stop]['wavelength', 0 ]

sds_sel

#%%

sds_sel.plot()

#%%

sds_sel.mean('time')

#%%

s2 = sds_sel.copy()

s2.variances = None

s2

s2.plot()

#%%

s2.mean('time')