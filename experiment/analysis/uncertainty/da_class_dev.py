#%%
import numpy as np
import xarray as xr
from mhdlab.xr_utils import WeightedMeanAccessor


ds = xr.tutorial.open_dataset("air_temperature")

ds2 = ds.wma.initialize_stat_dataset('air', 'time')

ds2

#%%

ds3 = ds2.wma.calc_weighted_mean('lat')  

ds3

#%%
ds3.wma.calc_weighted_mean('lon')
# %%

# make a test dataset with random data with specified mean and std and 2 dims

mean = 1
std = 0.1

ds = xr.Dataset()

ds['data'] = xr.DataArray(np.random.normal(mean, std, (10,10)), dims=['x', 'y'])

ds['data'].plot.hist()

#%%

ds2 = ds.wma.initialize_stat_dataset('data', 'x')

ds2

#%%

ds3 = ds2.wma.calc_weighted_mean('y')
ds3



# %%

ds3['stderr']*np.sqrt(ds3['count'])
