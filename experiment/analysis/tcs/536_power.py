#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')


from mhdpy.analysis import mws
from mhdpy.analysis import absem

# %%

tc = '536_power'

ds_absem = ppu.fileio.load_absem(tc)

ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='relative')
da_lecroy = ds_lecroy['AS']

ds_pd = ds_lecroy[['pd1','pd2', 'dpd1']].mean('run')
#%%

from pi_paper_utils import LASER_POWER, LASER_AREA

fluence = LASER_POWER/LASER_AREA

fluence = fluence.to('mJ/cm^2')
fluence = fluence.round(3)

da_lecroy = da_lecroy.assign_coords(power=da_lecroy['power']*fluence)
da_lecroy.coords['power'].attrs['units'] = 'mJ/cm^2'
da_lecroy.coords['power'].attrs['long_name'] = 'Fluence'

ds_pd = ds_pd.assign_coords(power=ds_pd['power']*fluence)
ds_pd.coords['power'].attrs['units'] = 'mJ/cm^2'
ds_pd.coords['power'].attrs['long_name'] = 'Fluence'

#%%

da_lecroy.plot(hue='power', row='run')

plt.yscale('log')

plt.xlim(-1,50)
plt.ylim(1e-5,1e-1)

#%%
da_lecroy.mean('run').plot(hue='power')

plt.yscale('log')
plt.xlim(-1,40)
plt.ylim(1e-5,1e-1)

#%%
from mhdpy.xr_utils.stats import WeightedMeanAccessor

#TODO: add weighted mean dataarray acessor and tests

ds_stat = da_lecroy.to_dataset().wma.initialize_stat_dataset('AS', 'run')

ds_stat = ds_stat.drop(0,'power')

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
plt.xlim(-1,50)
plt.ylim(1e-5,1e-1)

plt.ylabel('AS')
plt.xlabel('Time [us]')

plt.legend(title='Fluence\n[mJ/cm^2]')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWS_power_time.png'))

# %%

# da_max = da_lecroy.sel(time=slice(-1,1)).max('time')
da_max = da_lecroy.mws._pulse_max()

da_max.attrs['long_name'] = 'Max AS'

da_max
#%%
from mhdpy.plot import dropna

g = da_max.plot(hue ='run_plot', x='power', marker='o')
plt.yscale('log')


# %%

plt.figure(figsize=(4,3))

g = da_max.plot(hue ='run_plot', x='power', marker='o')
plt.yscale('log')
plt.xscale('log')

plt.title('')
plt.xlabel('Fluence [mJ/cm^2]')

plt.ylim(5e-5,1e-1)

plt.savefig(pjoin(DIR_FIG_OUT, 'MWS_power_max.png'))
#%%

da_max2 = da_max.mean('run')
da_max2.plot(marker='o')
plt.yscale('log')
plt.xscale('log')


#%%

delta_pd1.mean('run').plot(marker='o')


#%%

from lmfit.models import PowerLawModel

model = PowerLawModel()

params = model.make_params()

params

#%%

da_fit = da_max2.sel(power=slice(1,100))

da_max2

result = model.fit(da_fit.values, x=da_fit['power'].values)

params = result.params


#%%

vals  = model.eval(params, x=da_fit['power'].values)

plt.plot(da_fit['power'], vals)

plt.plot(da_fit['power'], da_fit.values)


plt.yscale('log')
plt.xscale('log')

#%%

da_mean = da_max.mean('run')
da_std = da_max.std('run')

plt.errorbar(da_mean['power'], da_mean, yerr=da_std, marker='o', capsize=5, label='Data')

plt.plot(da_fit['power'], vals, label='Fit')

# add fit params as text

# plt.text(0.1, 0.9, "Model: Ax^b", transform=plt.gca().transAxes)
# plt.text(0.1, 0.8, f'b: {params["exponent"].value:.2f}', transform=plt.gca().transAxes)
# plt.text(0.1, 0.7, f'A: {params["amplitude"].value:.2e}', transform=plt.gca().transAxes)

plt.text(0.1, 0.9, "Model: Ax^b", transform=plt.gca().transAxes)
plt.text(0.1, 0.8, f'b: {params["exponent"].value:.2f} ± {params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.1, 0.7, f'A: {params["amplitude"].value:.2e} ± {params["amplitude"].stderr:.2e}', transform=plt.gca().transAxes)

plt.yscale('log')
plt.xscale('log')

plt.legend(loc='lower right')

plt.xlabel('Fluence [mJ/cm^2]')
plt.ylabel('$\Delta AS_{max}$')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWS_power_max_fit.png'))

#%%


ds_pd['pd1'].plot(hue='power')
plt.xlim(-1,20)
plt.yscale('log')
plt.ylim(1e-4,)
plt.title('')
plt.ylabel('PD on-peak [V]')

plt.savefig(pjoin(DIR_FIG_OUT, 'pd1_power_trace.png'))

#%%

ds_pd['dpd1'].plot(marker='o')

#%%

da_fit = ds_pd['dpd1'].sel(power=slice(1,100))

model = PowerLawModel()

params = model.make_params()

result = model.fit(da_fit.values, x=da_fit['power'].values)

result

#%%

da_fit

#%%

vals  = model.eval(result.params, x=da_fit['power'].values)

plt.plot(da_fit['power'].values, vals, label='Fit')

plt.plot(da_fit['power'].values, da_fit.values, label='Data', marker='o')

# plt.yscale('log') plt.xscale('log')

plt.xlabel('Fluence [mJ/cm^2]')
plt.ylabel('Delta PD1 [V]')

plt.legend(loc='lower right')

plt.text(0.1, 0.9, "Model: Ax^b", transform=plt.gca().transAxes)
plt.text(0.1, 0.8, f'b: {result.params["exponent"].value:.2f} ± {result.params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.1, 0.7, f'A: {result.params["amplitude"].value:.2e} ± {result.params["amplitude"].stderr:.2e}', transform=plt.gca().transAxes)


plt.savefig(pjoin(DIR_FIG_OUT, 'delta_pd1_power_fit.png'))



#%%

da_fit = da_lecroy.copy()

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit)

#TODO: sterr is nan where ds_p is not?
ds_p['krb'] = ds_p['krb'].where(~ds_p_stderr['krb'].isnull())


#%%

ds_p['krb'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

#%%

ds_p['dne'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')
plt.xscale('log')

plt.ylim(1e0,)

#%%

ds_mws_fit.to_array('var').plot(hue='var', row='power', col='run', x='time')

plt.yscale('log')

#%%

# ds_mws_fit.isel(time=0).count('mnum')['AS_all'].plot(hue='run_plot', x='power', marker='o')

# plt.yscale('log')

# plt.ylabel("mnum count")

#%%[markdown]

# # Exponential fitting

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

# da_fit = da_lecroy.mean('mnum')
da_fit = da_lecroy.copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

# ds_mws_fit = xr.merge([ds_mws_fit, da_fit.rename('AS_all')])
#%%

ds_mws_fit.to_array('var').plot(hue='var', row='power', col='run', x='time', yscale='log', sharey=False)

plt.yscale('log')

plt.xlim(0,30)

# %%

ds_p['decay'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

plt.ylim(0.1,50)

plt.xscale('log')

plt.title('')
plt.ylabel('Exponential time constant [us]')

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mws_decay_power.png'))

#%%

ds_p.mean('run')['decay'].plot(marker='o')

plt.ylim(1,10)

# %%


da_fit = da_lecroy.copy()

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3 

pipe_fit_mws_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(0, 'us'), Quantity(25, 'us'))
pipe_fit_mws_3.fit_prep_kwargs['pre_norm_cutoff'] = None

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_3(da_fit)

# %%

da_fit.mws.fit_prep()
# %%

ds_mws_fit#.sel(power=1)['AS_all'].plot(hue='run_plot', x='time')

ds_mws_fit.to_array('var').plot(hue='var', row='power', col='run', x='time', yscale='log', sharey=False)

plt.yscale('log')

plt.xlim(0,30)
# %%

ds_p_sel = ds_p.isel(power=slice(3,6))

ds_p_sel

#%%

dne = ds_p_sel['dne'].pint.quantify('1/um**3').pint.to('cm**-3')
krm = ds_p_sel['krm'].pint.quantify('1/us')


# plt.yscale('log')

#%%

fig, axes = plt.subplots(2,1, figsize=(5,4), sharex=True)

dne.plot(hue='run_plot', x='power', marker='o', ax=axes[0])
axes[0].get_legend().set_bbox_to_anchor((1, 1))
# axes[0].set_ylim(0,1e14)
krm.plot(hue='run_plot', x='power', marker='o', ax=axes[1])
axes[1].get_legend().remove()
axes[1].set_title('')

axes[0].set_ylabel('$\Delta n_e$ [cm^-3]')
axes[0].set_xlabel('')
axes[1].set_xlabel('Fluence [mJ/cm^2]')
axes[1].set_ylabel('$k_{r,m}$ [us^-1]')

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mws_dnedt_power.png'), bbox_inches='tight')


#%%

dne_mean = dne.mean('run')
dne_std = dne.std('run')

plt.errorbar(dne_mean['power'], dne_mean, yerr=dne_std, marker='o', capsize=5)


#%%

tau = 1/krm

tau_mean = tau.mean('run')
tau_std = tau.std('run')

plt.errorbar(tau_mean['power'], tau_mean, yerr=tau_std, marker='o', capsize=5)

#%%

tau
# %%
