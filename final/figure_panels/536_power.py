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

import matplotlib as mpl

# X axis limits
x_min = -1
x_max = 46
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
plt.xlim(x_min, x_max)
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

plt.text(0.05, 0.9, "Model: $Ax^b$", transform=plt.gca().transAxes)
plt.text(0.05, 0.8, f'b: {params["exponent"].value:.2f} ± {params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.05, 0.7, f'A: {params["amplitude"].value:.2e} ± {params["amplitude"].stderr:.2e}', transform=plt.gca().transAxes)

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
plt.xlim(x_min, x_max)
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

plt.text(0.05, 0.9, "Model: $Ax^b$", transform=plt.gca().transAxes)
plt.text(0.05, 0.8, f'b: {result.params["exponent"].value:.2f} ± {result.params["exponent"].stderr:.2f}', transform=plt.gca().transAxes)
plt.text(0.05, 0.7, f'A: {result.params["amplitude"].value:.2f} ± {result.params["amplitude"].stderr:.2f}', transform=plt.gca().transAxes)


plt.savefig(pjoin(DIR_FIG_OUT, 'delta_pd1_power_fit.png'))


# %%

#modified to be a single figure with subplots and mpl-created A, B, C... labels

mpl.rcParams.update({'font.size': 8})

fig, axes = plt.subplots(2, 1, figsize = (4, 5))




ms = 2 #Marker size for scatter plots
cs = 2 #Cap size for error bars

inset_loc = [0.6, 0.65, 0.3, 0.3] #arg for inset location and height

fit_export = [] #collect lines for exporting fit params to text file
#delta-AS vs time

powers = ds_stat['power'].values

for i, power in enumerate(powers):

    # plot with confidence interval

    ds_stat_sel = ds_stat.sel(power=power)
    ds_stat_sel['mean'].plot(label=f"{power:.1f}", ax = axes[0])

    axes[0].fill_between(ds_stat_sel.coords['time'], ds_stat_sel['mean'] - ds_stat_sel['std'], ds_stat_sel['mean'] + ds_stat_sel['std'], alpha=0.2)

    axes[0].set_title('')
    axes[0].set_xlabel('')

axes[0].set_yscale('log')
axes[0].set_xlim(x_min, x_max)
axes[0].set_ylim(1e-5,2e-1)

axes[0].set_ylabel('$\\Delta AS$')
axes[0].set_xlabel(r'Time [$\mathrm{\mu s}$]')

axes[0].set_xticklabels([])



# axes[0].legend(title=r'$\n Fluence \mathrm{[mJ/cm^2]}$', loc='upper left', bbox_to_anchor=(1.05, 0.2), borderpad = 0.1)


axes[0].legend(title='Fluence\n$\\mathrm{[mJ/cm^2]}$', loc='upper left', bbox_to_anchor=(1.05, 0.2), borderpad = 0.2)

#delta-AS vs fluence

ax = axes[0].inset_axes(inset_loc)
da_max = da_lecroy.mwt._pulse_max()
da_max.attrs['long_name'] = 'Max AS'
da_max_runavg = da_max.mean('run')


model = PowerLawModel()
da_fit = da_max_runavg.sel(power=slice(0.01,100))
result = model.fit(da_fit.values, x=da_fit['power'].values)
params = result.params

x_vals = da_fit['power'].values
x_eval = np.linspace(x_vals.min(), x_vals.max(), 100)
y_eval  = model.eval(params, x=x_eval)

da_mean = da_max.mean('run')
da_std = da_max.std('run')

ax.errorbar(da_mean['power'], da_mean, yerr=da_std, fmt='o', capsize = cs, markersize = ms, label='Data')

ax.plot(x_eval, y_eval, label='Fit')

ax.set_yscale('log')
ax.set_xscale('log')

# ax.legend(reverse = True)

ax.set_xlabel(r'Fluence [$\mathrm{mJ/cm^2}$]')
ax.set_ylabel('$\\Delta AS_{max}$')

#Store fit params
fit_export.append("delta-AS vs fluence fit parameters:")
fit_export.append(f'b: {params["exponent"].value:.2f} pm {params["exponent"].stderr:.2f}')
fit_export.append(f'A: {params["amplitude"].value:.2e} pm {params["amplitude"].stderr:.2e}')
fit_export.append('\n')

#delta-PD vs time
ax = axes[1]

da_plot = ds_pd['dpd1'].copy()
da_plot = da_plot.drop(0,'power')
da_plot = da_plot.assign_coords(power=[float(f'{p:.1f}') for p in da_plot.coords['power'].values])
da_plot.coords['power'].attrs = ds_pd['power'].attrs

da_plot.plot(hue='power', ax = ax, add_legend = False)
ax.set_xlim(x_min, x_max)
ax.set_yscale('log')
ax.set_ylim(1e-3,)
ax.set_title('')
ax.set_xlabel(r'Time [$\mathrm{\mu s}$]')
ax.set_ylabel('$\\Delta PD$ [mV]')


#delta-PD vs fluence



ax = axes[1].inset_axes(inset_loc)


da_fit = ds_pd['dpd1_max'].sel(power=slice(0.01,100))
model = PowerLawModel()
params = model.make_params()
result = model.fit(da_fit.values, x=da_fit['power'].values)



vals  = model.eval(result.params, x=da_fit['power'].values)

ax.plot(da_fit['power'].values, da_fit.values, label='Data', marker='o', markersize = ms)

ax.plot(da_fit['power'].values, vals, label='Fit')



ax.set_xlabel(r'Fluence [$\mathrm{mJ/cm^2}$]')
ax.set_ylabel('$\\Delta PD_{max}$ [mV]')

# ax.legend(loc='lower right')

ax.set_yscale('log')
ax.set_xscale('log')

#Store fit params
fit_export.append("delta-PD vs fluence fit parameters:")
fit_export.append(f'b: {result.params["exponent"].value:.2f} pm {result.params["exponent"].stderr:.2f}')
fit_export.append(f'A: {result.params["amplitude"].value:.2f} pm {result.params["amplitude"].stderr:.2f}')



fig.tight_layout()



fig.subplots_adjust(hspace = 0.02)

labels = ['A)','B)','C)','D)']

for ax, label in zip(axes, labels):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - .15, Y - .01, label)
    # ax.tick_params(axis='y', which = "both", colors='white')

fig.savefig(os.path.join(REPO_DIR, 'final','figures', 'Fig5_Laser_power_dependence.svg'))

with open(os.path.join(REPO_DIR, 'final','figures','Fig5_Laser_power_dependence_fitparams.txt'), 'w') as file:
    file.writelines(line + '\n' for line in fit_export)


