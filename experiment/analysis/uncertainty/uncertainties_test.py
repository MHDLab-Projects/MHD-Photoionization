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

