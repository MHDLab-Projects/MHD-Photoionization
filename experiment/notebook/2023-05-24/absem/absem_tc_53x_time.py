# %%

from mhdpy.analysis.standard_import import *

data_folder = mhdpy.fileio.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-05-24')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()


ds = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_alpha_fit.cdf'))


ds_mp = ds

tw= slice(Timestamp('2023-05-24 22:04:25.609789952'), Timestamp('2023-05-24 22:31:23.400811264'), None)

ds_mp = ds_mp.sel(time = tw)
ds_mp = ds_mp.dropna('mp', how='all').squeeze()


ds_mp

#%%
from pint import Quantity

da_signal = dsst['hvof']['CC_K_massFrac_in']
# da_signal = da_signal + Quantity(1 + 3/16, 'inches')

da_signal.name = 'kwt'

da_signal.sel(time=tw).plot(marker='o')

#%%

from mhdpy.mws_utils.coords import new_bin_coords, bin_coords_simple
from mhdpy.analysis.xr import assign_signal


# bins = np.linspace(-12.5,237.5,11)
# coords_grid = new_bin_coords(da_signal.values, bins)


grid_values = [0.0005, 0.001, 0.0018, 0.0031, 0.005, 0.01, 0.018]
coords_grid = bin_coords_simple(da_signal.values, grid_values)

da_signal.values = coords_grid


ds_mp = assign_signal(ds_mp, da_signal, 'time')
#%%
#TODO: Taken from mws utils
from mhdpy.mws_utils.coords import assign_mnum


ds_mp
ds_mp = ds_mp.rename(time = 'acq')

ds_mp = ds_mp.set_index(acq=['kwt'])
ds_mp = ds_mp.groupby('acq').apply(lambda x: assign_mnum(x, 'acq', 1 ))
ds_mp = ds_mp.set_index(acq='mnum', append=True)
ds_mp = ds_mp.unstack('acq').rename(acq_level_0= 'kwt')

# ds_mp = ds_mp.rename(acq='pos')

#%%

ds_mp


#%%


#%%

ds#.plot(hue='var', col='phi', row='kwt')


# %%
g = ds_mp.mean('mnum')[['data','fit']].sel(mp='barrel').to_array('var').plot(col='kwt',hue='var', ylim = (1e-3,2), yscale='log', figsize=(10,6), col_wrap=3)

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')



#%%



ds_p = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_p.cdf'))


ds_p['nK_m3'].plot(hue='mp')
plt.yscale('log')


#%%
ds_p = assign_signal(ds_p, da_signal, 'time')

ds_p = ds_p.rename(time = 'acq')

ds_p = ds_p.set_index(acq=['kwt'])
ds_p = ds_p.groupby('acq').apply(lambda x: assign_mnum(x, 'acq', 1 ))
ds_p = ds_p.set_index(acq='mnum', append=True)
ds_p = ds_p.unstack('acq').rename(acq_level_0= 'kwt')

# %%

from mhdpy.plot import dropna

g = ds_p.mean('mnum')['nK_m3'].plot(marker='o', col='mp')

dropna(g)

plt.ylim(1e19,3e22)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('Stage Position (mm)')

plt.savefig(pjoin(DIR_FIG_OUT,'nK_seedramp.png'))

# %%

ds_p.to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_53x.cdf'))