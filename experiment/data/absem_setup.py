# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem

from mhdpy.analysis.xr import interp_ds_to_var


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
from mhdpy.process.absem import calc_alpha_simple, reduce_switches, get_value_switches, assign_multiplexer_coord, downselect_num_acq

data_folder = os.path.join('munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

# Determine LED switching events
switches = get_value_switches(ds_absem.coords['led'].values, switch_to_vals=[0,1])
ds_absem = ds_absem.assign_coords(led_switch_num=('time', switches))

if has_multiplexer:
    ds_mp = assign_multiplexer_coord(
    ds_absem,
    mp_sent=dsst['multiplexer_send']['Position'].pint.dequantify(),
    mp_receive=dsst['multiplexer_receive']['Position'].pint.dequantify(),
    mp_port_names={1:'barrel', 2:'mw_horns'}
    )

    # Now we remove data when the multiplexer was switching, kept to allow for accurate determination of switching events
    ds_absem = ds_mp.where(ds_mp['mp'] != 'switch').dropna('time','all')

ds_absem = ds_absem.groupby('led_switch_num').apply(downselect_num_acq, num_acq=10)
ds_absem = ds_absem.dropna('time',how='all')

# Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
acq_groups = ['led_switch_num','led','time']
if 'mp' in ds_absem.coords: acq_groups.append('mp')

ds = ds_absem.set_index(acq_group=acq_groups) # Time has to be added here or it is retained as a dimension?
ds = ds.reset_index('time').reset_coords('time') 

ds_reduce_switches = reduce_switches(ds)

acq_groups.remove('led_switch_num')
ds_alpha = ds_reduce_switches.set_index(temp=acq_groups).unstack('temp')
ds_alpha = ds_alpha['counts_mean']

ds_alpha = ds_alpha.to_dataset('led').rename({0:'led_on', 1:'led_off'})
ds_alpha = interp_ds_to_var(ds_alpha, 'led_on')

#TODO: probably can do this above. 
if not 'mp' in ds_absem.coords:
    ds_alpha = ds_alpha.assign_coords(mp='barrel').expand_dims('mp')

time_window_calib = calib_timewindow

ds_alpha = calc_alpha_simple(ds_alpha, time_window_calib)

ds_alpha = ds_alpha.stack(acq=['time','mp']).reset_index('acq')

ds_alpha.to_netcdf(pjoin(data_folder, 'Munged','Spectral', 'ds_absem_mp.cdf'))