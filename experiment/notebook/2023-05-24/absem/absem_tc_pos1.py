# %%
from mhdpy.analysis.standard_import import *

from mhdpy.plot.common import xr_errorbar, xr_errorbar_axes


import mhdpy

data_folder = mhdpy.fileio.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-05-24')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

ds_calib = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_calib.cdf'))
ds = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_mean.cdf'))

ds['alpha'] = 1 - ds['diff']/ds_calib['diff']


ds_mp = ds

tw = slice(Timestamp('2023-05-24 21:28:10.185809408'), Timestamp('2023-05-24 21:34:34.578646784'), None)
# tw = slice(Timestamp('2023-05-18 22:44:30.975140096'), Timestamp('2023-05-18 22:54:56.400380672'), None)
ds_mp = ds_mp.sel(time = tw)
ds_mp = ds_mp.dropna('mp', how='all').squeeze()


ds_mp

#%%
from pint import Quantity

da_signal = dsst['motor']['Motor C Relative']
da_signal = da_signal + Quantity(1 + 3/16, 'inches')

da_signal.name = 'motor'

da_signal.sel(time=tw).plot(marker='o')

#%%

from mhdpy.mws_utils.coords import new_bin_coords, bin_coords_simple
from mhdpy.analysis.xr import assign_signal

# bins = np.linspace(-12.5,237.5,11)
# coords_grid = new_bin_coords(da_signal.values, bins)

grid_values = [0,25,50,75,100,125,150,175,200,225,250]
grid_values = [g + 30 for g in grid_values] #TODO: Hacky way to add zero offset
coords_grid = bin_coords_simple(da_signal.values, grid_values)

da_signal.values = coords_grid


ds_mp = assign_signal(ds_mp, da_signal, 'time')
#%%
#TODO: Taken from mws utils
from mhdpy.mws_utils.coords import assign_mnum


ds_mp
ds_mp = ds_mp.rename(time = 'acq')

ds_mp = ds_mp.set_index(acq=['motor'])
ds_mp = ds_mp.groupby('acq').apply(lambda x: assign_mnum(x, 'acq', 1 ))
ds_mp = ds_mp.set_index(acq='mnum', append=True)
ds_mp = ds_mp.unstack('acq').rename(acq_level_0= 'pos')

# ds_mp = ds_mp.rename(acq='pos')

#%%

ds_mp


#%%


ds_mp['led_off'].mean('mnum').plot(hue='pos', col='mp')


    # ds_sel.coords['motor'].attrs = dict(long_name='Position', units = 'mm')

# %%

g = ds_mp.to_array('var').sel(wavelength=slice(700,785)).mean('mnum').plot(hue='pos', col='mp', row='var', sharey=False)
for ax in g.axes[3]: ax.set_ylim(-0.5,1.2)

# %%

ds_fit = ds_mp.mean('mnum')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_corr = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

alpha_tc = alpha_corr
wls = alpha_tc.coords['wavelength'].values


alpha_tc.plot(hue='pos', col='mp')

plt.ylim(-1,1.2)
plt.xlim(760,775)


# %%

spectral_reduction_params_fp = os.path.join(REPO_DIR, 'experiment', 'metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength', how='all')
alpha_tc_red.name = 'alpha_red'

fits, ds_p, ds_p_stderr = analysis.xr.fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, alpha_tc_red, fits]).sel(wavelength=slice(760,775))


#%%

ds = ds_alpha[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds


#%%

ds#.plot(hue='var', col='phi', row='kwt')


# %%
g = ds_alpha[['alpha','alpha_fit']].sel(mp='mw_horns').to_array('var').plot(col='pos',hue='var', ylim = (1e-3,2), yscale='log', figsize=(10,6), col_wrap=3)

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')


# %%

ds_p['nK_m3'].plot(marker='o', col='mp')
plt.ylim(1e18,3e22)
plt.yscale('log')
# plt.xscale('log')
plt.xlabel('Stage Position (mm)')

plt.savefig(pjoin(DIR_FIG_OUT,'nK_pos.png'))

# %%


ds_p.to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_pos1.cdf'))