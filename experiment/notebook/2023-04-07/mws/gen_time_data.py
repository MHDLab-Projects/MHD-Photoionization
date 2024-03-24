"""Generate time curves to go into processed data tdms file"""

## TODO: Writing disabled for now... incorporate in post processing?
#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-04-07'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

# dsst = mhdpy.fileio.TFxr().as_dsst()

lecroy_munged_folder = pjoin(data_folder, 'Lecroy')
input_fns = [fn for fn in os.listdir(lecroy_munged_folder) if 'Silicon' not in fn]
input_fns = [fn for fn in input_fns if 'Nothing' not in fn]

from mhdpy.analysis import mws


dss = []

time_coords_from_first = True
time_coords = None

fps_in  = [os.path.join(lecroy_munged_folder, fn) for fn in input_fns]
fps_in = [fp_in for fp_in in fps_in if not os.path.isdir(fp_in)]

for i, fp_in in enumerate(fps_in):

    ds = xr.load_dataset(fp_in)

    if time_coords_from_first:
        if i == 0:
            time_coords = ds.coords['time']
        else:
            ds = ds.assign_coords(time=time_coords)

    ds = ds.sortby('acq_time')

    dss.append(ds)

    # break

ds = xr.concat(dss, 'acq_time')

ds = ds.sortby('acq_time')

ds = ds.mws.calc_AS_rel()

#%%


#%%

max_window = slice(-1,1)
pre_pulse = slice(-50,-1)

# ds['AS'].mean('acq_time').sel(time=pre_pulse).plot()
# ds['AS'].mean('acq_time').plot()

#%%

pulse_max = ds['AS'].sel(time=max_window).max('time')
pulse_max.name = 'pulse_max'

pre_pulse_avg =ds['AS'].sel(time=pre_pulse).mean('time') 
pre_pulse_avg.name = 'pre_pulse_avg'
pre_pulse_std =ds['AS'].sel(time=pre_pulse).std('time') 
pre_pulse_std.name = 'pre_pulse_std'

ds_out = xr.merge([
    pulse_max,
    pre_pulse_avg,
    pre_pulse_std
])

ds_out = ds_out.rename(acq_time =  'time')

#%%

ds_out['pulse_max'].plot()


#%%

# from mhdpy.fileio.tdms import ds_to_tdms
# from nptdms import TdmsWriter

# # fp_out = os.path.join(lecroy_munged_folder, 'test.tdms')
# fp_out = pjoin(data_folder, 'Processed_Data.tdms')

# with TdmsWriter(fp_out, 'a') as tdms_writer:
#     ds_to_tdms(ds_out, 'lecroy',tdms_writer)


