#%%

from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from lmfit.models import PowerLawModel

# %%

tc = '536_power'

ds_aas = ppu.fileio.load_aas(tc)

ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')
da_lecroy = ds_lecroy['dAS_abs']

ds_pd = ds_lecroy.mwt.calc_time_stats()[['pd1','pd2', 'dpd1', 'dpd1_max']].mean('run')
#%%

from pi_paper_utils.constants import LASER_POWER, LASER_AREA

fluence = LASER_POWER/LASER_AREA

fluence = fluence.to('mJ/cm^2')
fluence = fluence.round(3)

da_lecroy = da_lecroy.assign_coords(power=da_lecroy['power']*fluence)
da_lecroy.coords['power'].attrs['units'] = '$mJ/cm^2$'
da_lecroy.coords['power'].attrs['long_name'] = 'Fluence'

ds_pd = ds_pd.assign_coords(power=ds_pd['power']*fluence)
ds_pd.coords['power'].attrs['units'] = '$mJ/cm^2$'
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

# %%

# da_max = da_lecroy.sel(time=slice(-1,1)).max('time')
da_max = da_lecroy.mwt._pulse_max()

da_max.attrs['long_name'] = 'Max AS'

da_max
#%%
from mhdlab.plot import dropna

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

plt.savefig(pjoin(DIR_FIG_OUT, 'MWT_power_max.png'))
#%%

da_max2 = da_max.mean('run')
da_max2.plot(marker='o')
plt.yscale('log')
plt.xscale('log')


#%%

# look at power dependence after the second order effects

da_exp = da_lecroy.sel(time=slice(3, 5)).mean('time').mean('run')

da_exp.plot()

#%%

model = PowerLawModel()

params = model.make_params()

da_fit = da_exp.sel(power=slice(0.01,100))

result = model.fit(da_fit.values, x=da_fit['power'].values)

result

#%%

# Plot both exp and max on the same plot

fig, ax = plt.subplots()

da_max2.plot(marker='o', ax=ax, label='Max AS')
da_exp.plot(marker='o', ax=ax, label='exp region')

plt.yscale('log')
plt.xscale('log')

plt.legend()


#%%

ds_pd['dpd1_max'].plot(marker='o')

#%%

#%%[markdown]

# # Exponential fitting

#%%

from mhdlab.analysis.mwt.fitting import pipe_fit_exp

# da_fit = da_lecroy.mean('mnum')
da_fit = da_lecroy.copy()

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

# ds_mwt_fit = xr.merge([ds_mwt_fit, da_fit.rename('AS_all')])
#%%

ds_mwt_fit.to_array('var').plot(hue='var', row='power', col='run', x='time', yscale='log', sharey=False)

plt.yscale('log')

plt.xlim(0,30)

# %%

ds_p['decay'].plot(hue='run_plot', x='power', marker='o')

plt.yscale('log')

plt.ylim(0.05,100)

plt.xscale('log')

plt.title('')
plt.ylabel('Exponential time constant [us]')

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mwt_decay_power.png'))

#%%

ds_p.mean('run')['decay'].plot(marker='o')

plt.ylim(1,10)

# %%


da_fit = da_lecroy.copy()

from mhdlab.analysis.mwt.fitting import pipe_fit_mwt_3 

pipe_fit_mwt_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(0, 'us'), Quantity(25, 'us'))
pipe_fit_mwt_3.fit_prep_kwargs['pre_norm_cutoff'] = None

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_3(da_fit)

# %%

da_fit.mwt.fit_prep()
# %%

ds_mwt_fit#.sel(power=1)['AS_all'].plot(hue='run_plot', x='time')

ds_mwt_fit.to_array('var').plot(hue='var', row='power', col='run', x='time', yscale='log', sharey=False)

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

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mwt_dnedt_power.png'), bbox_inches='tight')


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
