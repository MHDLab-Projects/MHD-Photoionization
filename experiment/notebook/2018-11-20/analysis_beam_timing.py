#%%

from mhdlab.analysis.standard_import import *


ds = xr.load_dataset(pjoin('output', 'ds_beam_timing.cdf'))
da = ds['counts']

da
# %%

da.mean('estime').mean('gatedelay').plot()
# %%

da_time = da.mean('x').mean('y')
# da_time = da.max('x').max('y')

da_time

#%%

da_time.isel(estime=[0,10,20,30,40,50]).plot(hue='estime')

#%%

da_time.plot()

#%%

da_time['estime'].plot(marker='o')
# %%
tw = slice(Timestamp('2018-11-21 00:25:00'), Timestamp('2018-11-21 00:30:00'))

da_time.sel(estime=tw).plot()

#%%

laser_profile_1 = da_time.sel(estime=tw).mean('estime')

laser_profile_1.plot(marker='o')
# %%

laser_profile_1.to_series().to_csv(pjoin('output', 'laser_profile_1.csv'))

#%%

tw = slice(Timestamp('2018-11-21 00:10:00'), Timestamp('2018-11-21 00:20:00'))
laser_profile_2 = da_time.sel(estime=tw).mean('estime')

laser_profile_2.plot()
# %%


laser_profile_2.to_series().to_csv(pjoin('output', 'laser_profile_2.csv'))

# %%
