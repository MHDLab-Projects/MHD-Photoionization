#%%


from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '516_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]

ds_lecroy = ds_lecroy.isel(run=0) # There is only one 516

# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)

#%%

ds_absem

#%%

#%%

ds_absem.sel(mp='barrel')['alpha'].mean('mnum').mean('wavelength').dropna('run', how='all')

#%%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='motor', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)

#%%

ds_absem['alpha'].sel(mp='mw_horns').mean('mnum').mean('run').plot(hue='motor', x='wavelength')

plt.gca().get_legend().remove()

plt.ylim(-0.1,0.1)


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds_absem.mean('mnum')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_1(ds_fit, spect_red_dict)

#%%

ds = ds_alpha_fit[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%

g  = ds_p['nK_m3'].plot(hue='run_plot', x='motor',col='mp', marker='o')

dropna(g)

plt.yscale('log')


# %%

ds_p['nK_m3'].sel(mp='mw_horns').plot(hue='run_plot', x='motor', marker='o')

plt.yscale('log')

#%%[markdown]

# # 516 MWS Analysis

# only on 2023-05-24

#%%


da_max = ds_lecroy['AS'].mws._pulse_max()



da_max.mean('mnum').plot()


#%%

ds_lecroy['AS'].mean('mnum').plot(row='motor')

plt.yscale('log')


# %%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

da_fit = ds_lecroy['AS'].mean('mnum')

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(
                                da_fit, 
                                method='iterative', 
                                fit_timewindow = slice(Quantity(5, 'us'), Quantity(15, 'us'))
                                )


# %%

ds_p['decay'].plot(marker='o')

plt.ylabel('Decay Time (us)')

plt.savefig(pjoin(DIR_FIG_OUT, '516_mws_decaytime.png'))

#%%

ds_mws_fit[['AS_fit','AS_all']].to_array('var').plot(hue='var', col='motor', col_wrap=3)

plt.yscale('log')

plt.ylim(1e-4,)

plt.savefig(pjoin(DIR_FIG_OUT, '516_mws_AS_fit.png'))
# %%
