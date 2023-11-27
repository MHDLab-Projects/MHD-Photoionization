# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem

from mhdpy.analysis.xr import interp_ds_to_var
from mhdpy.process.absem import calc_alpha_simple, apply_mp


import json 
datestr = '2023-05-24'

with open(pjoin(REPO_DIR, 'experiment', 'data', 'settings.json')) as f:
    settings = json.load(f)[datestr]

start, stop = map(pd.Timestamp, settings['calib_timewindow'])
calib_timewindow = slice(start, stop)

has_multiplexer = settings['has_multiplexer']

#%%


from mhdpy.process.absem import assign_multiplexer_coord, get_value_switches, downselect_num_acq

#TODO: refactor to other functions to be able to see intermediate plots
def apply_mp2(dsst, ds_absem):
    """
    pipeline for adding a multiplexer coordinate and calculating statitstics over the switching events
    """
    ds_mp = assign_multiplexer_coord(
    ds_absem,
    mp_sent=dsst['multiplexer_send']['Position'].pint.dequantify(),
    mp_receive=dsst['multiplexer_receive']['Position'].pint.dequantify(),
    mp_port_names={1:'barrel', 2:'mw_horns'}
    )

    # Determine LED switching events
    switches = get_value_switches(ds_mp.coords['led'].values, switch_to_vals=[0,1])

    ds_mp = ds_mp.assign_coords(switch_num=('time', switches))

    # Now we remove data when the multiplexer was switching, kept to allow for accurate determination of switching events
    ds_mp = ds_mp.where(ds_mp['mp'] != 'switch').dropna('time','all')

    ds_mp = ds_mp.groupby('switch_num').apply(downselect_num_acq, num_acq=10)
    ds_mp = ds_mp.dropna('time',how='all')

    # Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
    ds = ds_mp.set_index(acq_group=['switch_num','led','mp', 'time']) # Time has to be added here or it is retained as a dimension?
    ds = ds.reset_index('time').reset_coords('time') 

    return ds


def reduce_switches(ds):

    ds_reduce_switches = xr.Dataset(coords = ds.groupby('acq_group').mean().coords)
    ds_reduce_switches['counts_mean'] = ds['counts'].groupby('acq_group').mean()
    ds_reduce_switches['counts_std'] = ds['counts'].groupby('acq_group').std()
    ds_reduce_switches['counts_first'] = ds['counts'].groupby('acq_group').first()
    ds_reduce_switches['time'] = ds['time'].groupby('acq_group').mean() # Each acquisition group assigned time based on this function

    # Get rid of switch num, as no longer needed, and rename index to be time (as determined by groupby function above)
    ds_reduce_switches = ds_reduce_switches.set_index(acq_group='time').rename(acq_group='time')
    ds_reduce_switches = ds_reduce_switches.reset_coords('switch_num')

    return ds_reduce_switches

#%%

# Start processing 

data_folder = os.path.join(REPO_DIR, 'experiment','data', 'munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

# Needed for scipp
ds_absem = ds_absem.assign_coords(led = ('time', ds_absem.coords['led'].astype(np.int64).values))

ds = apply_mp2(dsst, ds_absem)

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

ds2 = reduce_switches(ds)

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