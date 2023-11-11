
from mhdpy.analysis.standard_import import *


dates = [
'2023-04-07',
'2023-05-12',
'2023-05-18',
'2023-05-24'
]

process_tcs = [
['530', '536_Motor', '536_Motor_Power', '536_Motor_Power_Run2', '536_Power', 'SeedRamp', 'SeedRamp_Run2'],
['Motor_1', 'SeedRamp_2', 'SeedRamp_3', 'Power_1'],
['pos_1','pos_2', 'power', 'seedramp_1'],
['516', '536_power', '5x6', '5x3', 'seedramp_2']
]

for i, datestr in enumerate(dates):
    print(datestr)

    data_folder = os.path.join('munged',datestr)

    lecroy_munged_folder = pjoin(data_folder, 'Lecroy')

    process_fns = ['ds_{}.cdf'.format(tc) for tc in process_tcs[i]]
    input_fps = [pjoin(lecroy_munged_folder, fn) for fn in process_fns] 

    # Seems equal by eye to approach of manual first time coord setting, but xarray says not equal. Think close enough for this situation. TODO: reexamine and impelent to munging or only once in trc processing pipeline. 
    dss = [xr.load_dataset(fp) for fp in input_fps]
    ds = xr.concat(dss, 'acq_time', join='override')

    ds = ds.sortby('acq_time')

    ds.to_netcdf(pjoin(data_folder, 'ds_lecroy_time.cdf'))


