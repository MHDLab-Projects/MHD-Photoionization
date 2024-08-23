#%%[markdown]

# # Coordinate Processing Overview

#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

import re
from collections import defaultdict

from mhdpy.fileio.ct import load_df_cuttimes
from mhdpy.fileio.tdms import tdms2ds

dsst = mhdpy.fileio.TFxr(pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')).as_dsst(convert_to_PT=False)

coords_to_assign = tdms2ds(pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms'), convert_to_PT=False)


ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'ds_absem.cdf'))
ds_absem = ds_absem.set_index(acq=['time','mp']).unstack('acq')
ds_absem = ds_absem.rename(time='acq_time')


coords_orig = xr.merge([
    dsst['hvof']['CC_K_massFrac_in'],
    dsst['motor']['Motor C Relative'],
    dsst['hvof']['CC_equivalenceRatio']
])

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'ct_sequence.csv')
df_cuttimes = mhdpy.fileio.load_df_cuttimes(fp_cuttimes)

cuttimes = df_cuttimes.ct.slice_list()

#%%[markdown]

# We begin with data for the 'test case setpoint channels' which are the meausred time-series data for the setpoints 

#%%

coords_orig

#%%[markdown]

# We bin these coordinates (with bins defined in mhdpy.coords.gen_coords_to_assign_1) to group together all datapoints associatd with setpoint plateaus
# 
# Currently the data is all replaced with the average of all data in the bin, not the nominal setpoint
#
# The data is also scaled (e.g. K mass fract to K wt %) and/or offset (e.g. Add offset to stage position )
#
# Note that this data is not used in any calculation, just generating indexes.
#
# TODO: Move this to a more logical location, bins should probably be in metadata
# 
# also shade the bins in the plot below

#%%

tw = cuttimes[14] #2023-05-24 5x3_pos

fig, axes = plt.subplots(3)

coords_orig.sel(time=tw)['CC_K_massFrac_in'].dropna('time',how='all').plot(ax=axes[0])
ta = axes[0].twinx()
coords_to_assign.sel(time=tw)['kwt'].plot(color='r', ax=ta)

axes[0].set_ylim(0.00095, 0.00105)
ta.set_ylim(0.095, 0.105)


coords_orig.sel(time=tw)['Motor C Relative'].dropna('time',how='all').plot(ax=axes[1])
ta = axes[1].twinx()
coords_to_assign.sel(time=tw)['motor'].dropna('time',how='all').plot(color='r', ax=ta)

axes[1].set_ylim(0, 250)
ta.set_ylim(0,250)

coords_orig.sel(time=tw)['CC_equivalenceRatio'].dropna('time',how='all').plot(ax=axes[2], label='original')
ta = axes[2].twinx()
coords_to_assign.sel(time=tw)['phi'].plot(color='r', ax=ta, label='binned')

axes[2].set_ylim(0.5,1.5)
ta.set_ylim(0.5,1.5)

axes[2].legend()
ta.legend(bbox_to_anchor=[0.5,0,0,1])

#%%[markdown]

# These binned coordinates are added along the time dimension multindex
#
# For example, the absorption emission data. Currently indexed by time, wavelength, and multiplexer positon (mp). Multiplexer position was added by the same process. 

#%%

ds_absem

#%%

fig, axes = plt.subplots(2)

tw = cuttimes[0]

coords_orig.sel(time=tw)['CC_K_massFrac_in'].dropna('time',how='all').plot(ax=axes[0])

ta = axes[0].twinx()

coords_to_assign.sel(time=tw)['kwt'].plot(color='r', ax=ta)


ds_absem_sel = ds_absem.sel(acq_time=tw).sel(mp='barrel')
ds_absem_sel['led_on'].mean('wavelength').plot( ax=axes[1], marker='o')

#%%

#Functions expecting a dictionary of datarrays...
coord_signal_dict = {k: coords_to_assign[k].dropna('time',how='all') for k in coords_to_assign.data_vars}

from mhdpy.coords import assign_signal

ds = assign_signal(ds_absem_sel, coord_signal_dict['kwt'], timeindex='acq_time')

ds
#%%[markdown]

# Now xarray unstack method can be used to conevrt the 'acq_time' dimension into a multiindexed dimension
#
# First we convert 'acq_time' into a multindex. 

#%%

ds = ds.set_index(acq_time='kwt', append=True)

ds

#%%

ds = ds.unstack('acq_time').rename(acq_time_level_0='acq_time')

ds

#%%[markdown]

# Now the functionality of xarray can be used

#%%

ds['led_on'].mean('acq_time').plot(hue='kwt')

#%%[markdown]

# However, when trying to make a multindimensional array with more dimensions (e.g. 2D sweep, multiple experiment dates) there is a problem with keeping the time index. 
#
# For each given acquition the time is unique, so the final dataset is 'sparse' with a lot of nans
#
# My solution is instead to convert the time index into a measurement number for each group of setpoint coordinates. So there will be less missing values because each test case will have mnum=1, for example. There will still be missing for a test case when its maximum mnum (e.g. 20) is less than the maximum across all test cases, but this still greatly compresses the dataset. 

#%%

from mhdpy.coords import assign_coords_multi

ds = assign_coords_multi(ds_absem_sel, coord_signal_dict, min_mnum=10)

ds

#%%

ds['led_on'].mean('mnum').plot(hue='kwt')

# %%[markdown]

# The final datatset for seedramps does this processing over all dates

#%%

tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))

ds_absem

#%%[markdown]

# The experiment date and run number are 'stacked' back into a mulindex, as this represents a given 'run' of a test case
#
# The 'run_plot' is a quirk needed for plotting with the line hue as the run (TOOD: link github issues)

#%%

ds_absem = ds_absem.xr_utils.stack_run()

ds_absem

#%%[markdown]

# We can calculate the number of datapoints contributing to each combination of coordinates

#%%

counts = ds_absem['led_on'].count('mnum').isel(wavelength=0) # Counts are all same for each wavelength

counts

#%%

g = counts.plot(col='mp', hue='run_plot', x='kwt', marker='o')

plt.xscale('log')
g.axes[0,0].set_ylabel('Counts')

#%%

counts.plot.hist()
plt.ylabel("Frequency of Counts")
# %%
