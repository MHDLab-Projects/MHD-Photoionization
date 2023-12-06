#%%


from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

import re
from collections import defaultdict

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
from mhdpy.fileio.tdms import tdms2ds
from mhdpy.mws_utils.coords import gen_coords_to_assign_1, assign_coords_multi

dsst = mhdpy.fileio.TFxr(pjoin(DIR_PROC_DATA, 'dsst.tdms')).as_dsst()

coords_to_assign = tdms2ds(pjoin(DIR_PROC_DATA, 'dst_coords.tdms'))


ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_absem.cdf'))
ds_absem = ds_absem.set_index(acq=['time','mp']).unstack('acq')
ds_absem = ds_absem.rename(time='acq_time')

# %%


coords_orig = xr.merge([
    dsst['hvof']['CC_K_massFrac_in'],
    dsst['motor']['Motor C Relative'],
    dsst['hvof']['CC_equivalenceRatio']
])

#%%


# fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'experiment_timewindows.csv')
# df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

# cuttimes = extract_cuttime_list(df_exptw)

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'cuttimes.csv')
df_cuttimes = mhdpy.fileio.load_df_cuttimes(fp_cuttimes)

cuttimes = extract_cuttime_list(df_cuttimes)

#%%

fig, axes = plt.subplots(2)

tw = cuttimes[0]

coords_orig.sel(time=tw)['CC_K_massFrac_in'].dropna('time',how='all').plot(ax=axes[0])

ta = axes[0].twinx()

coords_to_assign.sel(time=tw)['kwt'].plot(color='r', ax=ta)


ds_absem['led_on'].sel(acq_time=tw).mean('wavelength').plot(hue='mp', ax=axes[1], marker='o')
# %%
