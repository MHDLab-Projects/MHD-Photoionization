#%%

from mhdpy.analysis.standard_import import *
import re
from collections import defaultdict

from mhdpy.io import load_df_cuttimes, extract_cuttime_list
from mhdpy.io import ds_to_tdms, TdmsWriter
from mhdpy.mws_utils.coords import gen_coords_to_assign_1, assign_coords_multi

dsst = mhdpy.io.TFxr(pjoin(DIR_PROC_DATA, 'dsst.tdms')).as_dsst()

df_cuttimes = load_df_cuttimes(pjoin(DIR_PROC_DATA, 'cuttimes.csv'), reduce_columns=False)
df_cuttimes = df_cuttimes.sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index(['Event', 'date'])

cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_absem.cdf'))
ds_absem = ds_absem.rename(time='acq_time')
ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_lecroy.cdf'))
#%%
process_tcs = [
    '53x',
    '536_pos',
    '536_power',
    '5x6_pos',
    '5x3_pos',
    '536_pos_power'
    ]


process_regex = ["({})_(\d)".format(s) for s in process_tcs]


def process_ds(ds, coords_to_assign, date, run_num, tc_base, min_mnum=10):

    ds = assign_coords_multi(ds, coords_to_assign, min_mnum=min_mnum)
    ds = ds.assign_coords(date=[date])

    #TODO: believe these ohter coordinates are showing up because of low min_mnum
    if tc_base == '53x':
        if 'phi' in ds.dims:
            ds = ds.sel(phi=0.8, method='nearest')


    run_num = int(run_num)
    ds = ds.assign_coords(run_num=[int(run_num)])#.expand_dims('run_num')
    ds = ds.stack(run=['date','run_num'])
    return ds

coords_to_assign = gen_coords_to_assign_1(dsst)

#%%

ds_coords = xr.merge(coords_to_assign.values())

with TdmsWriter('proc_data/dst_coords.tdms', 'w') as tw:
    ds_to_tdms(ds_coords, 'coords', tw)


#%%
dss_absem = defaultdict(list)
dss_lecroy = defaultdict(list)

for (tc, date), row in df_cuttimes.iterrows():

    sl = slice(
        np.datetime64(row['Start Time']), 
        np.datetime64(row['Stop Time'])
        )


    for regex in process_regex:
        m = re.match(regex, tc)

        if m:
            tc_base, run_num = m.groups()

            ds_absem_sel = ds_absem.sel(acq_time = sl)
            ds_absem_sel = process_ds(ds_absem_sel, coords_to_assign, date, run_num, tc_base, min_mnum=None) #Absem aready downselected...
            dss_absem[tc_base].append(ds_absem_sel)

            ds_lecroy_sel = ds_lecroy.sel(acq_time = sl)
            ds_lecroy_sel = process_ds(ds_lecroy_sel, coords_to_assign, date, run_num, tc_base, min_mnum=10)
            dss_lecroy[tc_base].append(ds_lecroy_sel)


#%%

output_dir = pjoin(DIR_PROC_DATA, 'absem')
if not os.path.exists(output_dir): os.mkdir(output_dir)

for tc in dss_absem:
    ds = xr.concat(dss_absem[tc], 'run')
    ds = ds.unstack('run')

    ds.to_netcdf(pjoin(output_dir, '{}.cdf'.format(tc)))

# %%

output_dir = pjoin(DIR_PROC_DATA, 'lecroy')
if not os.path.exists(output_dir): os.mkdir(output_dir)

for tc in dss_lecroy:
    ds = xr.concat(dss_lecroy[tc], 'run')
    ds = ds.unstack('run')

    ds.to_netcdf(pjoin(output_dir, '{}.cdf'.format(tc)))