#%%

from mhdpy.analysis.standard_import import *

data_folder = mhdpy.fileio.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-05-18')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
df_cuttimes = load_df_cuttimes('cuttimes.csv').sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

#%%

process_tcs = ['pos_1','pos_2', 'power', 'seedramp_1']
process_fns = ['ds_{}.cdf'.format(tc) for tc in process_tcs]

lecroy_munged_folder = pjoin(data_folder, 'Lecroy')

input_fns = [fn for fn in os.listdir(lecroy_munged_folder) if fn in process_fns]
input_fps = [pjoin(lecroy_munged_folder, fn) for fn in input_fns] 

# Seems equal by eye to approach of manual first time coord setting, but xarray says not equal. Think close enough for this situation. TODO: reexamine and impelent to munging or only once in trc processing pipeline. 
dss = [xr.load_dataset(fp) for fp in input_fps]
ds = xr.concat(dss, 'acq_time', join='override')

ds = ds.sortby('acq_time')

ds.to_netcdf(pjoin(DIR_PROC_DATA, 'ds_lecroy_time.cdf'))
ds

#%%

from mhdpy.mws_utils.coords import gen_coords_to_assign_1, assign_coords_multi

coords_to_assign = gen_coords_to_assign_1(dsst)

for tc, row in df_cuttimes.iterrows():

    output_fp = os.path.join(DIR_PROC_DATA,'{}.cdf'.format(tc))

    sl = slice(
        np.datetime64(row['Start Time']), 
        np.datetime64(row['Stop Time'])
        )
    ds_sel = ds.sel(acq_time = sl)

    if not len(ds_sel.coords['acq_time'].values):
        print("No acquisitions for tc: {}".format(tc))
        continue

    ds_sel = assign_coords_multi(ds_sel, coords_to_assign)

    print(tc)
    print(ds_sel.coords)

    ds_sel.to_netcdf(output_fp)
# %%
