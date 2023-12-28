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
df_cuttimes = load_df_cuttimes(pjoin(REPO_DIR, 'experiment', 'metadata', 'cuttimes.csv'), reduce_columns=False)

# add tc and run_num columns. tc is everything before the last underscore in 'Event' column
df_cuttimes['tc'] = df_cuttimes['Event'].str.extract(r'(.*)_\d$')
df_cuttimes['run_num'] = df_cuttimes['Event'].str.extract(r'.*_(\d)$')


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

df_cuttimes = df_cuttimes.sort_values('Start Time').set_index(['tc', 'run_num' ,'date'])


#%%
def process_ds(ds, coords_to_assign, date, run_num, tc_base, min_mnum=10):

    ds = assign_coords_multi(ds, coords_to_assign, min_mnum=min_mnum)

    #TODO: believe these ohter coordinates are showing up because of low min_mnum
    if tc_base == '53x':
        if 'phi' in ds.dims:
            ds = ds.sel(phi=0.8, method='nearest')

    ds = ds.assign_coords(date=[date])
    ds = ds.assign_coords(run_num=[int(run_num)])#.expand_dims('run_num')
    ds = ds.stack(run=['date','run_num'])
    return ds


dss_absem = defaultdict(list)
dss_lecroy = defaultdict(list)

for (tc_base, run_num, date), row in df_cuttimes.iterrows():

    sl = slice(
        np.datetime64(row['Start Time']), 
        np.datetime64(row['Stop Time'])
        )

    ds_absem_sel = ds_absem.sel(acq_time = sl)
    #TODO: can we add 'mnum_counts' as a variable in a general way and downselect later?
    ds_absem_sel = process_ds(ds_absem_sel, coords_to_assign, date, run_num, tc_base, min_mnum=2) #Absem aready downselected...
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