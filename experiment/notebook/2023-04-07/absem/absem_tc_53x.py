
# %%
from mhdpy.analysis.standard_import import *
import seaborn as sns
from mhdpy.plot.common import xr_errorbar, xr_errorbar_axes

ds = xr.load_dataset(pjoin(DIR_PROC_DATA, 'ds_alpha.cdf')).squeeze()
#%%

from mhdpy.io import TFxr, gen_path_date


data_folder = gen_path_date('2023-04-07')
dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()

#%%


#%%
from mhdpy.io import load_df_cuttimes, extract_cuttime_list
df_cuttimes = load_df_cuttimes('cuttimes_53x_2.csv').sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')
cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

# da_ct = xr.DataArray(cuttimes, coords={'kwt': df_cuttimes.index.values}, dims=['kwt'])

from mhdpy.analysis.ct import gen_da_ct_data

kwt = dsst['hvof']['CC_K_massFrac_in']

da_ct = gen_da_ct_data(cuttimes, {'kwt': {'data': kwt}})

#TODO: round value not working in gen_da_ct_data....
new_kwt = [round(v,3) for v in da_ct.coords['kwt'].values]
da_ct = da_ct.set_index(ct='kwt').rename(ct='kwt')

da_ct = da_ct.assign_coords(kwt = new_kwt)

da_ct

# %%
from mhdpy.analysis.ct import gen_da_ct_data, assign_tc_general 

ds_tc = assign_tc_general(ds, da_ct)

#%%

fig, axes = plt.subplots(3, sharex=True)

for i, kwt in enumerate(ds_tc.coords['kwt']):
    da = ds_tc['led_off'].sel(kwt=kwt).dropna('time','all')

    for time in da.coords['time']:
        da.sel(time=time).plot(ax=axes[i], color='blue', alpha = 0.01)


plt.xlim(769,771)

#%%
fig, axes = plt.subplots(3, sharex=True)

for i, kwt in enumerate(ds_tc.coords['kwt']):
    da = ds_tc['led_off'].sel(kwt=kwt).dropna('time','all')

    df = da.drop('kwt').to_dataframe()

    sns.lineplot(df, x='wavelength', y='led_off', ax=axes[i], errorbar="sd")
    plt.xlim(769,771)

#%%

df = ds_tc['led_off'].to_dataframe()

sns.lineplot(df, x='wavelength', y='led_off', hue='kwt', errorbar="sd")

plt.xlim(769,771)

#%%

da_mean = ds_tc['led_off'].mean('time')
da_std = ds_tc['led_off'].std('time')

fig, axes = plt.subplots(3, sharex=True)

for i, ax in enumerate(axes):
    xr_errorbar_axes(
        da_mean.isel(kwt=i),
        da_std.isel(kwt=i),
        ax
        )

plt.xlim(769,771)

# %%

g = ds_tc.to_array('var').sel(wavelength=slice(700,785)).mean('time').plot(col='kwt', row='var', sharey=False)
for ax in g.axes[3]: ax.set_ylim(-0.5,1.2)

# %%

ds_fit = ds_tc.mean('time')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_corr = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

alpha_tc = alpha_corr
wls = alpha_tc.coords['wavelength'].values

alpha_tc.plot(hue='kwt')

plt.ylim(-1,1.2)
plt.xlim(760,775)

plt.axhline(0)

# %%

spectral_reduction_params_fp = os.path.join(REPO_DIR,'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()
print('Reducing alpha with following data reduction parameters: ')
print(spect_red_dict)
alpha_tc_red = analysis.spectral.alpha_cut(alpha_tc,**spect_red_dict).dropna('wavelength','all')
alpha_tc_red.name = 'alpha_red'

fits, ds_p, ds_p_stderr = analysis.xr.fit_da_lmfit(alpha_tc_red, final_model, pars, 'wavelength', wls)
ds_p['nK_m3'].attrs = dict(long_name='$n_{K,expt}$', units = '$\\#/m^3$')
fits.name = 'alpha_fit'

ds_alpha = xr.merge([alpha_tc, alpha_tc_red, fits]).sel(wavelength=slice(760,775))


ds_p.to_netcdf(pjoin(DIR_DATA_OUT, 'ds_p_53x.cdf'))
#%%

ds = ds_alpha[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds

# %%
g = ds_alpha[['alpha','alpha_fit']].to_array('var').plot(col='kwt',hue='var', ylim = (1e-3,2), yscale='log', figsize=(10,6))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')

# %%


xr_errorbar(ds_p['nK_m3'], ds_p_stderr['nK_m3'])

plt.ylim(5e20,3e22)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('K Mass Fraction')


plt.savefig(pjoin(DIR_FIG_OUT,'nK_53x.png'))
# %%
