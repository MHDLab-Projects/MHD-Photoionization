# %% [markdown]
#  # 2022-10-05 Continous seed feeding absorption emission measurements
# 

# %% [markdown]
#  ## Load data

# %%
from mhdpy.analysis.standard_import import *

from mhdpy.analysis.xr import interp_ds_to_var
from mhdpy.io import TFxr, gen_path_date
from mhdpy.process.absem import load_absem

data_folder = gen_path_date('2023-04-07')
dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

fp = os.path.join(data_folder,'Munged','Spectral' ,'absem.tdms')
ds_absem = load_absem(fp)

#%%

# timewindow = slice(Timestamp('2023-04-07 19:44:18.315159552'), Timestamp('2023-04-07 20:41:51.547015936'), None)
# ds_absem = ds_absem.sel(time=timewindow)


#%%
ds_alpha = ds_absem.set_index(temp=['time','led']).unstack('temp')
ds_alpha = ds_alpha['counts']
ds_alpha = ds_alpha.to_dataset('led').rename({0:'led_on', 1:'led_off'})
ds_alpha = interp_ds_to_var(ds_alpha, 'led_on')

#%%

ds_alpha = ds_alpha.assign_coords(mp='barrel').expand_dims('mp')
ds_alpha 

#%%
from mhdpy.process.absem import calc_alpha_simple

time_window_calib = slice(Timestamp('2023-04-07 20:38:26.554598656'), Timestamp('2023-04-07 20:38:41.251181312'), None)
ds_alpha = calc_alpha_simple(ds_alpha, time_window_calib)

#%%
plt.figure()
ds_alpha['alpha'].mean('time').plot(hue='mp')

plt.savefig(pjoin(DIR_FIG_OUT, 'alpha_avg.png'))
#%%

ds_alpha.to_netcdf(pjoin(DIR_PROC_DATA, 'ds_alpha.cdf'))