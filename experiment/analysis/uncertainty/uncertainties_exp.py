#%%

import uncertainties as unc
from uncertainties import unumpy as unp
import numpy as np
import xarray as xr

#%%

n1 = unc.ufloat(1,0.1)
n2 = unc.ufloat(1,0.1)

n1/n2
#%%

n1 + n2



#%%

arr = np.array([unc.ufloat(1,0.1), unc.ufloat(2,0.2), unc.ufloat(3,0.3)])
arr

# %%

unp.std_devs(arr)

#%%

unp.nominal_values(arr)

#%%

np.nanmean(arr)

#%%

s = arr[0] + arr[1] + arr[2]

s/3

#%%

arr + arr

#%%

# test an array of ufloats with 0 uncertainty

arr = np.array([unc.ufloat(1,0), unc.ufloat(2,0), unc.ufloat(3,0)])

np.nanmean(arr)


#%%

# make a two dimensional test xarray dataarray with floats and two dims

arr = np.array([[1,2,3],[4,5,6],[7,8,9]])

da = xr.DataArray(arr, dims=['x','y'], coords={'x':np.arange(3), 'y':np.arange(3)})

da

# %%

# convert to uncertainties array

da_unc = unp.uarray(da.values, 1)

da_unc = xr.DataArray(da_unc, dims=da.dims, coords=da.coords)

da_unc

#%%

da2 = xr.apply_ufunc(np.nanmean, da_unc, input_core_dims=[['x']], output_core_dims=[[]], kwargs={'axis':-1})
da2

#%%

# def calc_std():

#     return xr.apply_ufunc(unp.std_devs, da_unc, input_core_dims=[['x']], output_core_dims=[[]], kwargs={'axis':-1})


da_temp = xr.apply_ufunc(unp.nominal_values, da_unc)
da_std = xr.apply_ufunc(np.std, da_temp, input_core_dims=[['x']], output_core_dims=[[]], kwargs={'axis':-1}) 

da_std

#%%

da2


#%%

# This was the old code testing the uncertianties package

# # testing uncertianties package
# #TODO: move somewhere else, but refactor fitting first?

# from uncertainties import unumpy as unp

# ds_nK_unstacked = ds_nK.unstack('run')

# arr = unp.uarray(ds_nK_unstacked['mean'].values, ds_nK_unstacked['std'].values)

# # For some reason this appears to have the same order as arr
# # dims = list(ds_nK_unstacked.indexes.keys())
# dims = ds_nK_unstacked['mean'].dims #TODO: appears conversion to datarray is needed for correct dim order..

# coords = {dim: ds_nK_unstacked.coords[dim].values for dim in dims}
# da_nK_uc = xr.DataArray(arr, dims=dims, coords=coords)

# da_nK_uc

# #%%

# da = da_nK_uc.sel(kwt=1, method='nearest').sel(mp='barrel').stack(run=['date','run_num'])

# arr = da.values

# np.nanmean(arr, axis=-1)

# #%%

# # TODO: can't make sense of the std dev

# # narr = [v.nominal_value for v in arr]
# # np.nanmean(narr)
# # np.nanstd(narr)

# narr = [v.std_dev**2 for v in arr]
# np.sqrt(np.nanmean(narr))

# #%%

# da = da_nK_uc.stack(run=['date','run_num'])

# da = xr.apply_ufunc(np.nanmean, da, input_core_dims=[['run']], output_core_dims=[[]], kwargs={'axis':-1})

# da_nominal = xr.apply_ufunc(unp.nominal_values, da)
# da_std = xr.apply_ufunc(unp.std_devs, da)

# ds_nK_uc = xr.merge([da_nominal.rename('mean'), da_std.rename('std')])
# ds_nK_uc    

