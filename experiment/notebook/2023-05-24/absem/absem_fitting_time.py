#%%

from mhdpy.analysis.standard_import import *

data_folder = mhdpy.io.gen_path('sharepoint', 'Data Share', 'MHD Lab', 'HVOF Booth', '2023-05-24')

dsst = mhdpy.io.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

ds_alpha = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_alpha.cdf'))



#%%

ds_fit = ds_alpha

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

#%%

#Downselect to remove data where there are no peaks... TODO: revisit
alpha_tc = alpha_tc.where(alpha_tc.sel(wavelength=slice(760,775)).max('wavelength') > 0.5).dropna('time','all')
# alpha_tc = alpha_tc.sel(time=alpha_tc.coords['time'].values[::100])

#%%

alpha_tc

#%%

spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength','all')
alpha_tc_red.name = 'alpha_red'

wls = alpha_tc.coords['wavelength'].values
fits, ds_p, ds_p_stderr = analysis.xr.fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
# ds_p.coords['phi'].attrs = dict(long_name='Total Mass Flow', units = 'gram/second')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, alpha_tc_red, fits]).sel(wavelength=slice(760,775))




#%%

ds = ds_alpha[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
# ds = ds.to_array('var')

ds
# %%
# ds.plot(hue='var', col='mp', row='time')
# %%
# ds_p['nK_m3'].plot(hue='mp')

ds.to_netcdf(pjoin(DIR_PROC_DATA, 'ds_alpha_fit.cdf'))
ds_p.to_netcdf(pjoin(DIR_PROC_DATA, 'ds_p.cdf'))
ds_p_stderr.to_netcdf(pjoin(DIR_PROC_DATA, 'ds_p_stderr.cdf'))