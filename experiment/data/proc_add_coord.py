#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import re
from collections import defaultdict

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
from mhdpy.fileio.tdms import ds_to_tdms, TdmsWriter
from mhdpy.coords import gen_coords_to_assign_1, assign_coords_multi

# Absorption emsssion and lecroy
ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_absem.cdf'))
ds_absem = ds_absem.set_index(acq=['time','mp']).unstack('acq')
ds_absem = ds_absem.rename(time='acq_time')
ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_lecroy.cdf'))


dsst = mhdpy.fileio.TFxr(pjoin(DIR_PROC_DATA, 'dsst.tdms')).as_dsst()
# Generate binned coordinates to assign along time dimension
coords_to_assign = gen_coords_to_assign_1(dsst)

dss = [da.to_dataset(name=name) for name, da in coords_to_assign.items()]
ds_coords = xr.merge(dss)

with TdmsWriter('proc_data/dst_coords.tdms', 'w') as tw:
    ds_to_tdms(ds_coords, 'coords', tw)

#%%

# Cuttimes 
df_cuttimes = load_df_cuttimes(pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_sequence.csv'), reduce_columns=False)

# add tc and run_num columns. tc is everything before the last underscore in 'Event' column
df_cuttimes['tc'] = df_cuttimes['Event'].str.extract(r'(.*)_\d$')
df_cuttimes['run_num'] = df_cuttimes['Event'].str.extract(r'.*_(\d)$').astype(int)


process_tcs = [
    '53x',
    '536_pos',
    '536_power',
    '5x6_pos',
    '5x3_pos',
    '536_pos_power'
    ]

# Only keep columns with tc in process_tcs
df_cuttimes = df_cuttimes[df_cuttimes['tc'].isin(process_tcs)]



#%%

def assign_coords_timewindows(ds, coords_to_assign, df_cuttimes, min_mnum=None, timeindex='acq_time', downselect_dict=None):
    """
    Splits ds into time windows and applies coords_to_assign. 

    min_mnum: minimum number of measurements that must be present in each timewindow 
    downselect_dict: dictionary of coordinates to downselect to a signle value before concatenating. indexes cannot vary between timewindows
    """

    #TODO: generalize to not require date/run_num
    required_cols = ['Start Time', 'Stop Time', 'date', 'run_num']
    for col in required_cols:
        assert col in df_cuttimes.columns, "df_cuttimes needs to have {} column for this function".format(col)

    df_cuttimes = df_cuttimes.sort_values('Start Time')#.set_index(['tc', 'run_num' ,'date'])

    dss = []
    
    for idx, row in df_cuttimes.iterrows():

        sl = slice(
            np.datetime64(row['Start Time']), 
            np.datetime64(row['Stop Time'])
            )

        ds_tw = ds.sel({timeindex: sl})

        ds_tw = assign_coords_multi(ds_tw, coords_to_assign, min_mnum=min_mnum)

        if downselect_dict is not None:
            for coord, value in downselect_dict.items():
                if coord in ds_tw.dims:
                    ds_tw = ds_tw.sel({coord: value}, method='nearest')

        ds_tw = ds_tw.assign_coords(date=[row['date']])
        ds_tw = ds_tw.assign_coords(run_num=[row['run_num']])

        ds_tw = ds_tw.stack(run=['date','run_num']) #TODO: this can not be replaced by xr_utils.stack_run, debug

        dss.append(ds_tw)

    ds_out = xr.concat(dss, 'run')

    return ds_out
    

dss_absem = defaultdict(list)
dss_lecroy = defaultdict(list)

for tc_base, df_cuttimes_tc in df_cuttimes.groupby('tc'):

    #TODO: believe these ohter coordinates are showing up because of low min_mnum
    if tc_base == '53x':
        # if 'phi' in ds.dims:
        downselect_dict = {'phi': 0.8}
    else:
        downselect_dict = None


    #TODO: can we add 'mnum_counts' as a variable in a general way and downselect later?
    ds_absem_sel = assign_coords_timewindows(ds_absem, coords_to_assign, df_cuttimes_tc, min_mnum=2, downselect_dict=downselect_dict)
    ds_lecroy_sel = assign_coords_timewindows(ds_lecroy, coords_to_assign, df_cuttimes_tc, min_mnum=10 , downselect_dict=downselect_dict)


    dss_absem[tc_base].append(ds_absem_sel)
    dss_lecroy[tc_base].append(ds_lecroy_sel)
#%%

output_dir = pjoin(DIR_PROC_DATA, 'absem')
if not os.path.exists(output_dir): os.mkdir(output_dir)

for tc in dss_absem:
    ds = xr.concat(dss_absem[tc], 'run')

    # ds['led_on'] = ds['led_on']*2

    ds = ds.unstack('run')

    ds.to_netcdf(pjoin(output_dir, '{}.cdf'.format(tc)))

# %%

output_dir = pjoin(DIR_PROC_DATA, 'lecroy')
if not os.path.exists(output_dir): os.mkdir(output_dir)

for tc in dss_lecroy:
    ds = xr.concat(dss_lecroy[tc], 'run')
    ds = ds.unstack('run')

    ds.to_netcdf(pjoin(output_dir, '{}.cdf'.format(tc)))