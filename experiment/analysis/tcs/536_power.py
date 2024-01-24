#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.analysis import absem

# %%

tc = '536_power'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
da_lecroy = ds_lecroy.mws.calc_mag_phase_AS()['AS']

#%%

da_lecroy.mean('mnum').plot(hue='power', row='run')

plt.yscale('log')

plt.xlim(-1,50)
plt.ylim(1e-5,1e-1)
# %%

da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')

da_max.attrs['long_name'] = 'Max AS'

da_max
#%%
from mhdpy.plot import dropna

g = da_max.plot(hue ='run_plot', x='power', marker='o')


# %%

g = da_max.plot(hue ='run_plot', x='power', marker='o')
plt.yscale('log')
plt.xscale('log')

# %%

da_fit = da_lecroy.copy()

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)

#TODO: sterr is nan where ds_p is not?
ds_p['kr'] = ds_p['kr'].where(~ds_p_stderr['kr'].isnull())

#%%

ds_p['kr'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

#%%

ds_p['dne'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')
plt.xscale('log')

plt.ylim(1e0,)

#%%

ds_mws_fit.mean('mnum').to_array('var').plot(hue='var', row='power', col='run', x='time')

plt.yscale('log')

#%%

ds_mws_fit.isel(time=0).count('mnum')['AS_all'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

plt.ylabel("mnum count")

#%%[markdown]

# # Exponential fitting

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

# da_fit = da_lecroy.mean('mnum')
da_fit = da_lecroy.copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

# ds_mws_fit = xr.merge([ds_mws_fit, da_fit.rename('AS_all')])
#%%

ds_mws_fit.mean('mnum').to_array('var').plot(hue='var', row='power', col='run', x='time', yscale='log', sharey=False)

plt.yscale('log')

plt.xlim(0,30)

# %%

ds_p['decay'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

plt.ylim(1,20)

#%%

ds_p.mean('run')['decay'].plot(marker='o')

plt.ylim(1,10)
