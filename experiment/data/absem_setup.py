# %%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
from mhdpy.fileio.spectral import load_absem

from mhdpy.process.absem import assign_multiplexer_coord, get_value_switches
from mhdpy.process.absem import downselect_num_acq
from mhdpy.analysis.xr import interp_ds_to_var
from mhdpy.process.absem import calc_alpha_simple



def apply_mp(dsst, ds_absem):
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


    plt.figure()

    ds_plot = ds_mp.isel(time=slice(1000,1200))

    ds_plot['led'].plot(marker='o')

    plt.twinx()
    ds_plot['switch_num'].plot(color='r', marker='o')

    plt.savefig(pjoin(DIR_FIG_OUT, 'absem_switching_sample.png'))


    plt.figure()

    ds_mp['counts'].groupby('switch_num').count().plot.hist(label='before downselect')

    ds_mp = ds_mp.groupby('switch_num').apply(downselect_num_acq, num_acq=10)
    ds_mp = ds_mp.dropna('time',how='all')

    ds_mp['counts'].groupby('switch_num').count().plot.hist(label = 'after')

    plt.legend()

    plt.savefig(pjoin(DIR_FIG_OUT, 'absem_switching_counts.png'))


    # Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
    ds = ds_mp.set_index(acq_group=['switch_num','led','mp', 'time']) # Time has to be added here or it is retained as a dimension?
    ds = ds.reset_index('time').reset_coords('time') 

    ds_reduce_switches = xr.Dataset(coords = ds.groupby('acq_group').mean().coords)
    ds_reduce_switches['counts_mean'] = ds['counts'].groupby('acq_group').mean()
    ds_reduce_switches['counts_std'] = ds['counts'].groupby('acq_group').std()
    ds_reduce_switches['counts_first'] = ds['counts'].groupby('acq_group').first()
    ds_reduce_switches['time'] = ds['time'].groupby('acq_group').mean() # Each acquisition group assigned time based on this function

    # Get rid of switch num, as no longer needed, and rename index to be time (as determined by groupby function above)
    ds_reduce_switches = ds_reduce_switches.set_index(acq_group='time').rename(acq_group='time')
    ds_reduce_switches = ds_reduce_switches.reset_coords('switch_num')

    return ds_reduce_switches


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