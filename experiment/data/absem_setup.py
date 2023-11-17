# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem

from mhdpy.analysis.xr import interp_ds_to_var
from mhdpy.process.absem import calc_alpha_simple, apply_mp


import json
import argparse

parser = argparse.ArgumentParser(description='Process ABSEM data.')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

with open('settings.json') as f:
    settings = json.load(f)[datestr]

start, stop = map(pd.Timestamp, settings['calib_timewindow'])
calib_timewindow = slice(start, stop)

has_multiplexer = settings['has_multiplexer']


# Start processing 

data_folder = os.path.join('munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

if has_multiplexer:
    ds_reduce_switches = apply_mp(dsst, ds_absem)
    ds_alpha = ds_reduce_switches.set_index(temp=['time','mp','led']).unstack('temp')
    ds_alpha = ds_alpha['counts_mean']

else:
    ds_alpha = ds_absem.set_index(temp=['time','led']).unstack('temp')
    ds_alpha = ds_alpha['counts']

ds_alpha = ds_alpha.to_dataset('led').rename({0:'led_on', 1:'led_off'})
ds_alpha = interp_ds_to_var(ds_alpha, 'led_on')

#TODO: probably can do this above. 
if not has_multiplexer:
    ds_alpha = ds_alpha.assign_coords(mp='barrel').expand_dims('mp')

time_window_calib = calib_timewindow

ds_alpha = calc_alpha_simple(ds_alpha, time_window_calib)

ds_alpha.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_absem_mp.cdf'))