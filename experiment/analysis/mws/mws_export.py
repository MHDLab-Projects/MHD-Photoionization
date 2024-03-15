#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

fp = pjoin(REPO_DIR, 'experiment', 'data', 'munged', '2023-05-24', 'Lecroy', 'ds_seedramp_2.cdf')

# fp = r'C:\Users\ASPITARL\code\MHD-Photoionization-1\experiment\data\munged\2023-05-24\Lecroy\ds_seedramp_2.cdf'
# %%
tw = slice(Timestamp('2023-05-24 22:28:39.959077632'), Timestamp('2023-05-24 22:30:19.612336384'), None)

ds = xr.load_dataset(fp)

ds = ds.sel(acq_time=tw)
ds['acq_time'] = ds['acq_time'].astype('datetime64[ns]')

ds = ds.stack(temp=['acq_time','time'])
ds

#%%
temp_arr = ds['temp'].values

acq_times = [np.datetime64(x[0], 'ns') for x in temp_arr]
times = [x[1] for x in temp_arr]

# Convert to numpy timedelta64 objects
times_as_timedeltas = [np.timedelta64(int(x*1e9), 'ns') for x in times]
final_times = [x + y for x,y in zip(acq_times, times_as_timedeltas)]

# %%

df = ds.to_dataframe()

df['final_time'] = final_times

df = df.set_index('final_time', drop=True)

#%%

df

# %%
df.to_csv('output/mws_export.csv')
# %%