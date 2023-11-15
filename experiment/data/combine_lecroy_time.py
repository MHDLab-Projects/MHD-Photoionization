
from mhdpy.analysis.standard_import import *

import json
import argparse

parser = argparse.ArgumentParser(description='Make Lecroy time dataset')

parser.add_argument('-d', '--date', type=str, default=None,
                    help='Date to process, in format YYYY-MM-DD')

datestr = parser.parse_args().date

with open('settings.json') as f:
    settings = json.load(f)[datestr]

process_tcs = settings['mws_process_tcs']

# Start processing

data_folder = os.path.join('munged',datestr)

lecroy_munged_folder = pjoin(data_folder, 'Lecroy')

process_fns = ['ds_{}.cdf'.format(tc) for tc in process_tcs]
input_fps = [pjoin(lecroy_munged_folder, fn) for fn in process_fns] 

# Seems equal by eye to approach of manual first time coord setting, but xarray says not equal. Think close enough for this situation. TODO: reexamine and impelent to munging or only once in trc processing pipeline. 
dss = [xr.load_dataset(fp) for fp in input_fps]
ds = xr.concat(dss, 'acq_time', join='override')

ds = ds.sortby('acq_time')

ds.to_netcdf(pjoin(data_folder, 'ds_lecroy_time.cdf'))


