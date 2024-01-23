#%%
from mhdpy.analysis.standard_import import *
from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()

fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

dst_coords
# %%

fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'cuttimes_kwt.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes).sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

df_cuttimes['date'] = df_cuttimes['Start Time'].dt.date

#%%

hts = []

for idx, row in df_cuttimes.iterrows():
    timeslice = slice(row['Start Time'], row['Stop Time'])

    ht = dsst['calorimetry']['CC_heatTransfer'].sel(time=timeslice)
    ht = assign_signal(ht, dst_coords['kwt'].round(2), timeindex='time')

    hts.append(ht)

ht = xr.concat(hts, dim='time')

ht = unstack_multindexed_acq_dim(ht)

#%%

ht['CC_heatTransfer'].mean('mnum').plot(marker='o')

#%%

from mhdpy.plot.common import xr_errorbar

mean = ht['CC_heatTransfer'].mean('mnum')
std = ht['CC_heatTransfer'].std('mnum')

xr_errorbar(mean, std)