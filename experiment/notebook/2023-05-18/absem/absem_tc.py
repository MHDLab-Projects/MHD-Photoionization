# %%
from mhdpy.analysis.standard_import import *

ds_calib = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_calib.cdf'))
ds = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem_mean.cdf'))

ds['alpha'] = 1 - ds['diff']/ds_calib['diff']


ds_mp = ds

ds_mp

#%%
import pickle

with open('output/da_ct.pickle', 'rb') as f:
    da_ct = pickle.load(f)

da_ct = da_ct#.rename(kwt='ct').set_index(ct='ct')

da_ct = da_ct.rename(ct='kwt').set_index(kwt='kwt')

da_ct
# %%
from mhdpy.analysis.ct import gen_da_ct_data, assign_tc_general

ds_mp_tc = assign_tc_general(ds_mp, da_ct)

# %%

g = ds_mp_tc.to_array('var').sel(wavelength=slice(700,785)).mean('time').plot(hue='mp', col='kwt', row='var', sharey=False)
for ax in g.axes[3]: ax.set_ylim(-0.5,1.2)

# %%

ds_fit = ds_mp_tc.mean('time')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_corr = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

alpha_tc = alpha_corr
wls = alpha_tc.coords['wavelength'].values


alpha_tc.plot(hue='mp', col='kwt')

plt.ylim(-1,1.2)
plt.xlim(760,775)


# %%

spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
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

ds = ds_alpha.sel(mp='barrel')[['alpha', 'alpha_fit']].drop('mp')
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds


#%%

ds#.plot(hue='var', col='phi', row='kwt')


# %%
g = ds_alpha[['alpha','alpha_fit']].to_array('var').plot(col='kwt',row='mp',hue='var', ylim = (1e-3,2), yscale='log', figsize=(10,6))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')


# %%

ds_p['nK_m3'].plot(hue='mp', marker='o')
plt.ylim(1e18,3e22)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('K Mass Fraction')


# %%

g = ds_p['nK_m3'].plot(hue='phi', marker='o', col='mp')
plt.ylim(1e18,3e22)
plt.yscale('log')
plt.xscale('log')

for ax in g.axes.flatten():
    ax.set_xlabel('K Mass Fraction')

g.axes[0][0].set_title('Barrel Exit')
g.axes[0][1].set_title('Seed Port')

# %%

ds_p.to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_absem.cdf'))


# %%