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

#%%
da_lecroy.mean('mnum').mean('run').plot(hue='power')

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5,1e-1)

#%%
from mhdpy.xr_utils.stats import WeightedMeanAccessor

#TODO: add weighted mean dataarray acessor and tests

ds_stat = da_lecroy.mean('mnum').to_dataset().wma.initialize_stat_dataset('AS', 'run')

ds_stat = ds_stat.sel(power=[1,0.4,0.2,0.1])

ds_stat

#%%

plt.figure()

powers = ds_stat['power'].values

fig, ax = plt.subplots()

for i, power in enumerate(powers):

    # plot with confidence interval

    ds_stat_sel = ds_stat.sel(power=power)
    ds_stat_sel['mean'].plot(label=power)

    ax.fill_between(ds_stat_sel.coords['time'], ds_stat_sel['mean'] - ds_stat_sel['std'], ds_stat_sel['mean'] + ds_stat_sel['std'], alpha=0.2)

    ax.set_title('')
    ax.set_xlabel('')

    # unit_str = '[' + ds[var].attrs['units']+ ']' if 'units' in ds[var].attrs else '' 
    # ax.set_ylabel(ds[var].attrs['long_name'] + unit_str)
    
# axes[2].set_xlabel('Position [mm]')

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5,1e-1)

plt.ylabel('AS')
plt.xlabel('Time [us]')

plt.legend(title='Power [relative]')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWS_power_time.png'))

# %%

da_max = da_lecroy.mean('mnum').sel(time=slice(-1,1)).max('time')

da_max.attrs['long_name'] = 'Max AS'

da_max
#%%
from mhdpy.plot import dropna

g = da_max.plot(hue ='run_plot', x='power', marker='o')


# %%

plt.figure()

g = da_max.plot(hue ='run_plot', x='power', marker='o')
plt.yscale('log')
plt.xscale('log')

plt.title('')
plt.xlabel('Power [relative]')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWS_power_max.png'))

#%%[markdown]

# # Thermal acvitation model

#%%


da_max2 = da_max.mean('run')
da_max2.plot(marker='o')
# plt.yscale('log')

#%%

from pint import Quantity

f_A = Quantity(0.05, 'dimensionless')
laser_power = Quantity(10, 'mJ')

powers = da_max2['power'].values*laser_power*f_A

powers

#%%[markdown]

# Laser heating calculation

#%%

cantera_results_dir = pjoin(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')

ds_TP_params = xr.load_dataset(pjoin(cantera_results_dir, 'ds_TP_params.cdf'))

ds_sel = ds_TP_params[['Cp', 'rho']].sel(T=2000, P=1e5, phi=0.8, method='nearest')



Cp = Quantity(ds_sel['Cp'].mean('Kwt').item(), ds_sel['Cp'].attrs['units'])
rho = Quantity(ds_sel['rho'].mean('Kwt').item(), ds_sel['rho'].attrs['units'].replace("m3", 'm**3'))

Cp_vol = Cp * rho
Cp_vol

heated_volume = Quantity(1, 'cm^3') # Unknown/not well defined

Cp_final = Cp_vol * heated_volume

Cp_final.to('J/K')

#%%

deltaTs = powers / Cp_final

deltaTs.to('K')

#%%[markdown]

# Max ionized potassium calculation

#%%

Ei = Quantity(4.34, 'eV/particle')

n_max = powers/Ei/heated_volume

n_max.to('particle/cm^3').magnitude


#%%

ds_TP_params['rho']


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
