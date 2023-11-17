#%%


from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.mws_utils import calc_mag_phase_AS

plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.stack(run=['date','run_num'])
ds_absem = ds_absem.assign_coords(run_plot = ('run', ds_absem.indexes['run'].values))


ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.stack(run=['date','run_num'])
ds_lecroy = ds_lecroy.assign_coords(run_plot = ('run', ds_lecroy.indexes['run'].values))


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = calc_mag_phase_AS(ds_lecroy)#[['mag', 'phase','AS']]

# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)
#%%

ds_absem.sel(mp='barrel')['alpha'].mean('mnum').mean('wavelength').dropna('run','all')

#%%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='motor', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)
# %%
from mhdpy.plot import dropna

g = ds_lecroy['AS'].mean('mnum').plot(hue='run_plot', row='motor', x='time')

dropna(g)

#%%

# plt.figure(figsize=(4,6))

# da_sel = ds_lecroy['AS'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')
da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')

da_sel = da_sel.dropna('run','all')

da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))

# plt.yscale('log')

#%%

ds_fit = ds_absem.mean('mnum')

from mhdpy.analysis.spectral import model_blurredalpha_2peak, interp_alpha
from mhdpy import analysis
final_model, pars = model_blurredalpha_2peak()

beta = -np.log(1-ds_fit['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=slice(750,755)).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

#%%

#Downselect to remove data where there are no peaks... TODO: revisit
# alpha_tc = alpha_tc.where(alpha_tc.sel(wavelength=slice(760,775)).max('wavelength') > 0.5).dropna('time','all')
# alpha_tc = alpha_tc.sel(time=alpha_tc.coords['time'].values[::100])


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
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%

g  = ds_p['nK_m3'].plot(hue='run_plot', x='motor',col='mp', marker='o')

dropna(g)

plt.yscale('log')

#%%

plt.figure(figsize=(6,3))

da_sel = ds_p['nK_m3'].sel(mp='mw_horns')

da_sel = da_sel.dropna('run', 'all')

da_sel.plot(hue='run_plot', x='motor', marker='o')

plt.yscale('log')

#%%


mws_max = ds_lecroy['AS'].mean('mnum').max('time')
mws_max.name ='AS_max'
mws_max
nK = ds_p['nK_m3'].sel(mp='mw_horns')

mws_pp = ds_lecroy['mag_pp'].mean('mnum')

mws_pp
mws_pp.name = 'mag_pp'

#%%

ds = xr.merge([mws_max, mws_pp,nK]).sortby('motor').dropna('run','all')

ds['AS_max'].attrs = dict(long_name='AS Max')
ds['mag_pp'].attrs = dict(long_name='Mag. Pre Pulse', units='V')

ds
#%%

g = ds.to_array('var').plot(hue='run_plot', row='var', x='motor', sharey=False)

dropna(g)

#%%


fig, axes = plt.subplots(3, figsize=(4,8))

ds['AS_max'].plot(hue='run_plot', x='motor', ax=axes[0], marker='o')
ds['mag_pp'].plot(hue='run_plot', x='motor', ax=axes[1],marker='o')
ds['nK_m3'].plot(hue='run_plot', x='motor', ax=axes[2],marker='o')

for ax in axes:
    ax.get_legend().remove()
    ax.set_title('')

axes[2].set_xlabel('Position [mm]')

#%%

ds.unstack('run').plot.scatter(x='nK_m3', y='mag_pp')

plt.xscale('log')
# plt.yscale('log')

# %%


fp_cfd = pjoin(os.getenv('REPO_DIR'), 'modeling\cfd\output\line_profies.csv' )

df_cfd = pd.read_csv(fp_cfd, index_col=0)

df_cfd.index = df_cfd.index - df_cfd.index.values[0]
df_cfd.index = df_cfd.index*1000


ds_cfd = xr.Dataset.from_dataframe(df_cfd)


#TODO: Convert to number density
ds_cfd_norm = ds_cfd/ds_cfd.max()


#%%

ds_cfd.to_array('var').plot(hue='var', yscale='log')

#%%

fp_cfd_beam = pjoin(os.getenv('REPO_DIR'), r'modeling\cfd\output\beam_conv.cdf' )

ds_cfd_beam = xr.load_dataset(fp_cfd_beam)


xs = ds_cfd_beam.coords['x_beam'].values
xs = (xs - xs[0])*1000

ds_cfd_beam = ds_cfd_beam.assign_coords(x_beam = xs)

#TODO: Convert to number density
ds_cfd_beam_norm = ds_cfd_beam/ds_cfd_beam.max()


#%%

da = ds['nK_m3'].dropna('run','all')

g = da.plot(hue='run_plot', x='motor', marker='o')
plot.dropna(g)

ax1 = plt.gca()
ax1.get_legend().set_title('Experiment (date, #)')
# ax1.set_yscale('log')

plt.twinx()

ds_cfd_norm['K'].plot(color='black', label ='centerline')
ds_cfd_beam_norm['K'].plot(color='grey', label = 'beam convolution', marker='o')

ax2 = plt.gca()
ax2.legend(bbox_to_anchor=[0,0,1.0,0.7])
# ax2.set_yscale('log')

plt.xlim(0,310)
# %%

g = ds['AS_max'].plot(hue='run_plot', x='motor', marker='o')

leg = plt.gca().get_legend()
leg.set_title('Experiment (date, #)')
leg.set_bbox_to_anchor([0,0,1.7,1])
plt.twinx()

ds_cfd_norm['KOH'].plot(color='black', label ='centerline')
ds_cfd_beam_norm['KOH'].plot(color='gray', label = 'beam convolution')

ax2 = plt.gca()
ax2.legend(bbox_to_anchor=[0,0,1.8,0.3])

plt.xlim(0,310)
# %%
