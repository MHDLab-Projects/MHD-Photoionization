#%%[markdown]

# ## Absem Data reduction

# The determination of absorption near the peak is unreliable. Previously have
# used cut_alpha with a fixed wavelength range around each peak. This is not
# great because the width of the peak varies and the peak structure of alpha
# seems inconsistent and can even be negative. 
# 
# Experimenting with various ways of reducing data to wings
# 

#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred 
from mhdpy.analysis import absem

from mhdpy.xr_utils import XarrayUtilsAccessorCommon

# %%

tc = '53x'


ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem = ds_absem.stack(run = ['date','run_num']).dropna('run', how='all')

ds_absem = ds_absem.absem.calc_alpha()

ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem = ds_absem.drop('acq_time')

ds_absem

#%%

ds2 = ds_absem.unstack()
ds2 = ds2.xr_utils.groupby_dims([d for d in ds2.dims if d != 'wavelength']).apply(lambda x: x.absem.reduce_keep_wings()).unstack('temp')

#%%

ds3 = ds2.sel(date='2023-05-18', run_num=1, mp='barrel').groupby('kwt').apply(lambda x: x.xr_utils.assign_mnum('mnum'))
ds_cut = ds3.sel(mnum=1)

ds_cut = ds_cut.where(~ds_cut['alpha_red'].isnull()).dropna('wavelength', how='all')


# %%

ds_cut['led_off'].plot(row='kwt')
# %%

ds_alpha_fit, ds_p, ds_p_stderr = ds_cut.absem.perform_fit()
#%%

ds_alpha_fit.to_array('var').plot(row='kwt', hue='var')
# %%

ds_p['nK_m3'].plot(marker='o')

plt.xscale('log')
plt.yscale('log')
# %%
