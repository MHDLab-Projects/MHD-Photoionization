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
import pi_paper_utils as ppu

from mhdpy.analysis.absem.fitting import gen_model_alpha_blurred 
from mhdpy.analysis import absem

from mhdpy.xr_utils import XarrayUtilsAccessorCommon

dss_p = []
dss_p_stderr = []

# %%

ds_absem = ppu.fileio.load_absem('53x')

ds_absem = ds_absem.drop('acq_time')
ds_absem

seldict = dict(date='2023-05-18', run_num=1, mp='barrel')
ds_sel = ds_absem.sel(seldict).groupby('kwt').apply(lambda x: x.xr_utils.assign_mnum('mnum'))

# %%[markdown]

# ### old cut alpha method

# 

#%%

ds_fit = ds_sel.mean('mnum').absem.calc_alpha()

spectral_reduction_params_fp = os.path.join(REPO_DIR, 'experiment', 'metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_alpha_fit, ds_p, ds_p_stderr= absem.fitting.pipe_fit_alpha_1(ds_fit, spect_red_dict)
dss_p.append(ds_p.assign_coords(method='alpha_cut'))
dss_p_stderr.append(ds_p_stderr.assign_coords(method='alpha_cut'))
# %%

ds_alpha_fit['alpha_red'].plot(row='kwt')

#%%


#%%[markdown]

# ### Alpha Cut with removing alpha negative alpha

# This is the same as the old method, but with the addition of removing negative alpha values
# Removal was needed for the wings method, so this is for consistent comparison. 


#%%


ds_fit = ds_sel.mean('mnum').absem.calc_alpha()

ds_fit = ds_fit.absem.remove_beta_offset(beta_offset_wls=slice(750,755))

ds_fit = ds_fit.absem.drop_alpha_peaks_negative() # Addition to pipeline 1

ds_fit = ds_fit.absem.reduce_cut_alpha(**spect_red_dict)

model, params = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

ds_alpha_fit, ds_p, ds_p_stderr = ds_fit['alpha_red'].absem.perform_fit(model, params)
dss_p.append(ds_p.assign_coords(method='alpha_cut_no_neg'))
dss_p_stderr.append(ds_p_stderr.assign_coords(method='alpha_cut_no_neg'))

#%%
ds_fit


#%%[markdown]

# ### Reduce to wings method

# Find the wings of the spectrum by finding the first wavelength above a cutoff

#%%

# Need to calc alpha before droping for beta_offset. mw_horn data can be negative
ds_fit = ds_sel.absem.calc_alpha(beta_offset_wls=slice(750,755))

ds_fit = ds_fit.absem.drop_alpha_peaks_negative()


#%%

ds_sel['alpha'].sel(kwt=0.05).dropna('mnum').plot(row='mnum')

plt.ylim(-1.1,1.1)

#%%[markdown]

# ### Reduce to wings

#%%

ds2 = ds_fit.xr_utils.groupby_dims_wrapper(
    lambda x: x.absem.reduce_keep_wings(led_off_norm_cutoff=0.8), 
    [d for d in ds_fit.dims if d != 'wavelength']
    )

#%%

ds2.sel(kwt=0.1).isel(mnum=[0,5,10])['alpha_red'].plot(row='mnum')

#%%


# ds_cut = ds2.sel(mnum=1)
ds_cut = ds2.mean('mnum').absem.calc_alpha()
# ds_cut = ds_cut.where(~ds_cut['alpha_red'].isnull())#.dropna('wavelength', how='all')


# %%

ds_cut['alpha_red'].plot(row='kwt')
# %%

model, params = absem.gen_model_alpha_blurred(assert_xs_equal_spacing=False)

ds_alpha_fit, ds_p, ds_p_stderr = ds_cut['alpha_red'].absem.perform_fit(model, params)
dss_p.append(ds_p.assign_coords(method='wings'))
dss_p_stderr.append(ds_p_stderr.assign_coords(method='wings'))
#%%

ds_alpha_fit.to_array('var').plot(row='kwt', hue='var')

#%%[markdown]

# ### Reduce to wings with global fit

# This is the same as the previous method, but with a global fit.

#%%

da_fit = ds2.absem.calc_alpha()['alpha_red']

# da_fit = da_fit.dropna('mnum', how='all')

da_fit = da_fit.xr_utils.groupby_dims_wrapper(
    lambda x: x.dropna('mnum', how='all'),
    [d for d in da_fit.dims if d not in ['wavelength','mnum']]
    )


#%%

model, params = absem.gen_model_alpha_blurred(assert_xs_equal_spacing=False)

ds_alpha_fit, ds_p, ds_p_stderr = da_fit.absem.perform_fit(model, params, method='global')
dss_p.append(ds_p.assign_coords(method='wings_global'))
dss_p_stderr.append(ds_p_stderr.assign_coords(method='wings_global'))

# ds_alpha_fit = xr.merge([ds_alpha_fit, ds2.mean('mnum').absem.calc_alpha()['alpha']])


#%%


ds_alpha_fit.mean('mnum').to_array('var').plot(row='kwt', hue='var')

# %%

ds_p = xr.concat(dss_p, dim='method')
ds_p_stderr = xr.concat(dss_p_stderr, dim='method')

#%%

ds_p['nK_m3'].plot(marker='o', hue='method')

plt.xscale('log')
plt.yscale('log')
# %%

plot.common.xr_errorbar(ds_p['nK_m3'], ds_p_stderr['nK_m3'], huedim='method')

plt.xscale('log')
plt.yscale('log')

#%%

ds_p_stderr['nK_m3'].plot(marker='o', hue='method')
plt.ylabel('nK_m3_stderr')
# %%
