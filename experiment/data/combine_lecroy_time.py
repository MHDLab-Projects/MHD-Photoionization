#%%
from mhdpy.analysis.standard_import import *

import json
import argparse

parser = argparse.ArgumentParser(description='Make Lecroy time dataset')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

process_tcs = settings['mws_process_tcs']

# Start processing

data_folder = os.path.join('munged',datestr)

lecroy_munged_folder = pjoin(data_folder, 'Lecroy')

process_fns = ['ds_{}.cdf'.format(tc) for tc in process_tcs]
input_fps = [pjoin(lecroy_munged_folder, fn) for fn in process_fns] 

dss = []

for fp in input_fps:
    ds = xr.load_dataset(fp)

    #TODO: the starting time of each pulse is not the same for each run.
    # Could shift the starting time for each date, but dont sure how that would work with time_offset
    # for now just trimming a few us off the start and end of each pulse
    ds = ds.sel(time=slice(-48e-6, 48e-6))

    # This could be used to investigate data taken with the longer silicon timebase (e.g. 'nothing' data on 2023-05-24) but not necessary because resampling has already put on the same grid. Just need to slice down as above. 
    # if ds.attrs['TIMEBASE'] == '50_mV/div':
        
    dss.append(ds)

# Seems equal by eye to approach of manual first time coord setting, but xarray says not equal. Think close enough for this situation. TODO: reexamine and impelent to munging or only once in trc processing pipeline. 
ds = xr.concat(dss, 'acq_time', join='override')

ds = ds.sortby('acq_time')


# Time processing #TODO: integrate with lecroy munging or move there. 

time_offset = 0.88
ds = ds.assign_coords(time=ds.coords['time']*1e6 - time_offset)
ds.coords['time'].attrs['units'] = 'us'
ds.coords['time'].attrs['long_name'] = 'Time'

ds.to_netcdf(pjoin(data_folder, 'ds_lecroy_time.cdf'))



# %%

ds_out = ds.mws.calc_mag_phase_AS()[['mag_pp', 'mag_pp_std', 'AS']].mws._pulse_max()

mag_peak = ds.mws.calc_mag_phase_AS()['mag'].sel(time=slice(0, 0.1)).mean('time')
ds_out = ds_out.assign({'mag_peak':mag_peak})

ds_out = ds_out.rename({'acq_time':'time'})

from mhdpy.fileio.tdms import ds_to_tdms, TdmsWriter

with TdmsWriter(pjoin(data_folder, 'Processed_Data.tdms'), 'a') as tdms:
    ds_to_tdms(ds_out, group_name='mws_time', tdms_writer=tdms)


#%%
    
data_folder