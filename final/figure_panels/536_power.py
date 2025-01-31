# Power dependence figures
# This script pulls directly from munged data, possibly should move first part to dataset generation to be consistent. 

#%%

from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

from mhdlab.xr_utils.stats import WeightedMeanAccessor
from lmfit.models import PowerLawModel
# %%

tc = '536_power'

ds_absem = ppu.fileio.load_absem(tc)

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

ds_stat = da_lecroy.to_dataset().wma.initialize_stat_dataset('dAS_abs', 'run')

ds_stat = ds_stat.drop_sel(power=0)

ds_stat

#%%

plt.figure()

powers = ds_stat['power'].values

fig, ax = plt.subplots()

for i, power in enumerate(powers):

    # plot with confidence interval

    ds_stat_sel = ds_stat.sel(power=power)
    ds_stat_sel['mean'].plot(label=f"{power:.1f}")

    ax.fill_between(ds_stat_sel.coords['time'], ds_stat_sel['mean'] - ds_stat_sel['std'], ds_stat_sel['mean'] + ds_stat_sel['std'], alpha=0.2)

    ax.set_title('')
    ax.set_xlabel('')

plt.yscale('log')
plt.xlim(-1,50)
plt.ylim(1e-5,2e-1)

plt.ylabel('$\\Delta AS$')
plt.xlabel('Time [$\\mu s$]')

plt.legend(title='Fluence\n$[mJ/cm^2]$')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWT_power_time.png'))

# %%

# da_max = da_lecroy.sel(time=slice(-1,1)).max('time')
da_max = da_lecroy.mwt._pulse_max()

da_max.attrs['long_name'] = 'Max AS'

da_max_runavg = da_max.mean('run')

#%%

model = PowerLawModel()

da_fit = da_max_runavg.sel(power=slice(0.01,100))

result = model.fit(da_fit.values, x=da_fit['power'].values)

params = result.params


#%%

x_vals = da_fit['power'].values
x_eval = np.linspace(x_vals.min(), x_vals.max(), 100)

y_eval  = model.eval(params, x=x_eval)

# plt.plot(x_eval, y_eval, label='Fit')

# plt.plot(da_fit['power'], da_fit.values)
# plt.legend()
# plt.yscale('log')
# plt.xscale('log')

#%%

plt.figure()
da_mean = da_max.mean('run')
da_std = da_max.std('run')

plt.errorbar(da_mean['power'], da_mean, yerr=da_std, fmt='o', capsize=5, label='Data')

plt.plot(x_eval, y_eval, label='Fit')

plt.text(0.1, 0.9, "Model: $Ax^b$", transform=plt.gca().transAxes)
plt.text(0.1, 0.8, f'b: {params["exponent"].value:.2f} ± {params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.1, 0.7, f'A: {params["amplitude"].value:.2e} ± {params["amplitude"].stderr:.2e}', transform=plt.gca().transAxes)

plt.yscale('log')
plt.xscale('log')

plt.legend()

plt.xlabel('Fluence [$mJ/cm^2$]')
plt.ylabel('$\\Delta AS_{max}$')

plt.savefig(pjoin(DIR_FIG_OUT, 'MWT_power_max_fit.png'))


#%%

# Photodiode 

plt.figure()
da_plot = ds_pd['dpd1'].copy()

da_plot = da_plot.drop(0,'power')

da_plot = da_plot.assign_coords(power=[float(f'{p:.1f}') for p in da_plot.coords['power'].values])
da_plot.coords['power'].attrs = ds_pd['power'].attrs

da_plot.plot(hue='power')
plt.xlim(-1,50)
plt.yscale('log')
plt.ylim(1e-4,)
plt.title('')
plt.xlabel('Time [$\\mu s$]')
plt.ylabel('$\\Delta PD$ [mV]')

plt.gca().get_legend().remove()

plt.savefig(pjoin(DIR_FIG_OUT, 'dpd1_power_trace.png'))


#%%

da_fit = ds_pd['dpd1_max'].sel(power=slice(0.01,100))

model = PowerLawModel()

params = model.make_params()

result = model.fit(da_fit.values, x=da_fit['power'].values)

result


#%%

plt.figure()
vals  = model.eval(result.params, x=da_fit['power'].values)

plt.plot(da_fit['power'].values, da_fit.values, label='Data', marker='o')

plt.plot(da_fit['power'].values, vals, label='Fit')


# plt.yscale('log') plt.xscale('log')

plt.xlabel('Fluence [$mJ/cm^2$]')
plt.ylabel('$\\Delta PD_{max}$ [mV]')

plt.legend(loc='lower right')

plt.yscale('log')
plt.xscale('log')

plt.text(0.1, 0.9, "Model: $Ax^b$", transform=plt.gca().transAxes)
plt.text(0.1, 0.8, f'b: {result.params["exponent"].value:.2f} ± {result.params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.1, 0.7, f'A: {result.params["amplitude"].value:.2e} ± {result.params["amplitude"].stderr:.2e}', transform=plt.gca().transAxes)


plt.savefig(pjoin(DIR_FIG_OUT, 'delta_pd1_power_fit.png'))


# %%
