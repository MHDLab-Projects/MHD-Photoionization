#%%

from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date
data_folder = gen_path_date('2023-04-07')
dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

ds_calib_mp = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_calib.cdf'))
ds_mp_mean = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_mean.cdf'))
ds_mp_std = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_std.cdf'))

ds_mean_time = ds_mp_mean.mean('wavelength').dropna('time', how='all')

# %%


ds_mean_time['led_off'].dropna('time', how='all').plot()

#%%


ds_mean_time['led_on'].dropna('time', how='all').plot()


#%%
from pandas import Timestamp
from mhdpy.plot import dropna

plt.figure()


tw = slice(Timestamp('2023-04-07 20:05:16.134990336'), Timestamp('2023-04-07 20:35:43.603072768'), None)

ledon = ds_mp_mean['led_off']
ledon = ledon.sel(time=tw).dropna('time', how='all')
ledon = ledon.sel(wavelength = slice(765,772))


ledon = ledon.integrate('wavelength')

mean = ledon#.resample(time='60s').mean()


std = ds_mp_std['led_off'].sel(time=tw).dropna('time', how='all')
std = std.sel(wavelength = slice(765,772)).integrate('wavelength')

# std = ledon.resample(time='60s').std()

ratio = std/mean

ratio = ratio


fig, axes = plt.subplots(3 , figsize=(7,7), sharex=True)

g= ledon.plot(ax = axes[0])
dropna(g)
axes[0].set_ylabel("K Emission Intensity\n [counts nm/ms]")


ratio.plot(ax= axes[1])
axes[1].set_ylabel('K Emission \nStd. Dev./Mean Ratio')
axes[1].hlines(0.05, xmin=ratio.coords['time'].min(), xmax= ratio.coords['time'].max(), color='gray', alpha = 0.5)

da_kwt = dsst['hvof']['CC_K_massFrac_in']
da_kwt = da_kwt*100
da_kwt.sel(time=tw).plot(ax=axes[2])

axes[2].set_ylabel("Nominal K Mass \nPercentage Input [%]")
axes[2].set_ylim(-0.01,1.1)
# dsst['hvof']['CC_equivalenceRatio'].sel(time=tw).plot(ax=axes[2])

import matplotlib.dates as mdates
axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M')) 
axes[2].set_xlabel('Time')


plt.savefig(pjoin(DIR_FIG_OUT, 'absem_uncertainty_time.png'))


# %%
