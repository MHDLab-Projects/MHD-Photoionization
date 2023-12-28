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

from mhdpy.analysis.absem.fitting import interp_alpha
from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred, alpha_cut

# %%

tc = '53x'


ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
# Scipp cannot handle multindex, casts to a custom 'PyObject' dtype that does not go back to xarray

#TODO: having to do this on office comp?
ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

ds_absem = ds_absem.stack(run = ['date','run_num']).dropna('run', how='all')

ds_absem['diff'] = ds_absem['led_on'] - ds_absem['led_off']
ds_absem['alpha'] = 1 - ds_absem['diff']/ds_absem['calib']

ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem = ds_absem.drop('acq_time')

ds_absem
#%%

def reset_mnum(ds):
    ds = ds.dropna('mnum', how='all')

    ds = ds.assign_coords(mnum=range(len(ds.coords['mnum'])))

    return ds

ds_sel = ds_absem.sel(run=('2023-05-18', 1)).sel(mp='barrel')

ds_sel = ds_sel.groupby('kwt').apply(reset_mnum)

ds_sel = ds_sel.sel(mnum=1).dropna('kwt', how='all')

# %%

da_sel = ds_sel['led_off']

da_sel.plot(hue='kwt')

plt.xlim(764,772)

#%%

offset = da_sel.sel(wavelength=slice(750,752)).mean('wavelength')

ledoff_norm = (da_sel - offset)

ledoff_norm = ledoff_norm/ledoff_norm.max('wavelength')

ledoff_norm.plot(hue='kwt')

# ledoff_norm.sel(kwt=0.05).plot()

plt.xlim(764,772)


#%%

import numpy as np

ledoff_norm_cutoff = 0.5

def find_first_above_cutoff(data, cutoff):
    """
    Find the index of the first instance where the data is above a given cutoff.

    Parameters:
    data (numpy.ndarray): The data to search.
    cutoff (float): The cutoff value.

    Returns:
    int: The index of the first instance where the data is above the cutoff.
    """
    for i, value in enumerate(data):
        if value > cutoff:
            return i
    return None  # return None if no value is above the cutoff

def find_first_dataarray_processor(da):

    wl_idx = find_first_above_cutoff(da.values, ledoff_norm_cutoff)

    wl = da.wavelength[wl_idx]

    return wl

left_wl = ledoff_norm.groupby('kwt').apply(find_first_dataarray_processor)
right_wl = ledoff_norm.groupby('kwt').apply(lambda da: find_first_dataarray_processor(da.sortby('wavelength', ascending=False)))

right_wl


#%%


# Remove all data between the two wavelengths
ds_cut = ds_sel.where((ledoff_norm.wavelength < left_wl) | (ledoff_norm.wavelength > right_wl), drop=True)


# %%

ds_cut['led_off'].plot(row='kwt')
# %%

from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred

final_model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False,nan_policy='omit')


# %%

from mhdpy.xr_utils import fit_da_lmfit

ds_fit = ds_cut.copy()


beta = -np.log(1-ds_fit['alpha'])
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off))

wls = alpha_tc.coords['wavelength'].values
fits, ds_p, ds_p_stderr = fit_da_lmfit(alpha_tc, final_model, pars, 'wavelength', wls)

ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, fits])#.sel(wavelength=slice(760,775))

# %%

ds_alpha.to_array('var').plot(row='kwt', hue='var')
# %%

ds_p['nK_m3'].plot(marker='o')

plt.xscale('log')
plt.yscale('log')