#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu


data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_alpha_fit = xr.open_dataset(pjoin(data_directory, 'ds_pos_alpha_fit.cdf')).xr_utils.stack_run()

ds_lecroy = ds_lecroy.sel(phi=0.8, method='nearest')

#%%

ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('particle/m^3')


#%%

ds_p

#%%

da_sel = ds_lecroy['dAS_abs'].dropna('run',how='all')

da_sel.dropna('motor', how='all')

da_sel.plot(col='motor', col_wrap=3, hue='run_plot', x='time', figsize=(6,10))

#%%


fig, axes = plt.subplots(1, 2, figsize=(6,3), sharex=True, sharey=False)


da_plot = da_sel.sel(motor = 105, method='nearest').dropna('run',how='all')
da_plot.plot(hue='run_plot', x='time', ax=axes[0])

motor_expt = da_plot.coords['motor'].item()

axes[0].set_title('Position: {:.0f} mm'.format(motor_expt))
axes[0].set_ylim(-0.03,0.35)
axes[0].set_ylabel('$\Delta AS$')


da_sel_mean = da_sel.mean('run').sel(motor=180, method='nearest')
da_sel_std = da_sel.std('run').sel(motor=180, method='nearest')

motor_expt = da_sel_mean.coords['motor'].item()

da_sel_mean.plot(x='time', ax=axes[1])
axes[1].fill_between(da_sel_mean['time'].values, (da_sel_mean-da_sel_std).values, (da_sel_mean+da_sel_std).values, color='b', alpha=0.2)
axes[1].set_title('Position: {:.0f} mm'.format(motor_expt))
axes[1].set_ylabel('')
# axes[1].set_yscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel.png'), dpi=300, bbox_inches='tight')


#%%

ds_p_sel = ds_p.sel(phi=0.8, method='nearest') 

fig, axes = plt.subplots(1,2, figsize=(5,3), sharex=True)

for i, var in enumerate(['dAS_abs_max', 'SFR_abs']):
    ax = axes[i]

    ds_stat = ds_p_sel.wma.initialize_stat_dataset(var, 'run')

    # plot with confidence interval

    # ds_stat['mean'].plot(ax=ax)
    # ax.fill_between(ds_stat.coords['motor'], ds_stat['mean'] - ds_stat['std'], ds_stat['mean'] + ds_stat['std'], alpha=0.2)

    ax.errorbar(ds_stat.coords['motor'], ds_stat['mean'], yerr=ds_stat['std'], marker='o', capsize=5)

    ax.set_title('')
    ax.set_xlabel('Stage Position [mm]')

    unit_str = '[' + ds_p_sel[var].attrs['units']+ ']' if 'units' in ds_p_sel[var].attrs else '' 
    ax.set_ylabel(ds_p_sel[var].attrs['long_name'] + unit_str)
    
# axes[2].set_xlabel('Position [mm]')

ta = axes[0].twinx()

## photodiode 
# delta_pd1 = ds_p_sel['dpd1'].dropna('run', how='all')
# delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')
# delta_pd1.plot(ax=ta, color='red', marker='o')
# ta.set_ylabel('Delta PD1 [mV]')
# ta.set_title('')

motor_sel = [ 104.8, 178]

for mp in motor_sel:
    if mp == 178:
        axes[1].axvline(mp, color='gold', linestyle='--')
    else:
        axes[1].axvline(mp, color='gray', linestyle='--')


plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_stats.png'), dpi=300, bbox_inches='tight')

#%%

fig, axes = plt.subplots(2, figsize=(4.5,5), sharex=True)


# Phi =0.8
ax1 = axes[0]

ds_p_sel = ds_p.sel(phi=0.8, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.8, method='nearest')

phi_val_expt = ds_p_sel.coords['phi'].item()

nK_barrel_mean = ds_p_sel['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds_p_sel['nK_barrel'].mean('motor').std('run')

nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: removing last datapoint. Can this be fixed with calibraiton/processing?
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,220))

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-18 Run 1', ax=ax1)
g = nK_mw_horns.isel(run=1).plot(x='motor', marker='o', label='2023-05-18 Run 2', ax=ax1)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax1)
ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax1)

ax1.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg', capsize=5, )
ax1.set_title('phi = {}'.format(phi_val_expt))


# Phi = 0.65

ax2 = axes[1]

ds_p_sel = ds_p.sel(phi=0.6, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.6, method='nearest')

phi_val_expt = ds_p_sel.coords['phi'].item()

nK_barrel_mean = ds_p_sel['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds_p_sel['nK_barrel'].mean('motor').std('run')

nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: final datapoint for 516 has no potassium peaks. However, calibration is off and fit gives a large absorption value. Removing is probably best as any absorption is not meaninful anyway. 
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,200)) 

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-24 Run 1', ax=ax2)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax2)
ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax2)

ax2.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg', capsize=5, )
ax2.set_title('phi = {}'.format(phi_val_expt))

for ax in axes:

    ax.axvline(180, color='gold', linestyle='--')
    ax.set_ylim(1e16,3e22)
    ax.set_xlim(-10,250)
    ax.set_yscale('log')
    ax.legend()

ax2.set_xlabel('Position [mm]')
ax1.set_xlabel('')

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mws_cfd.png'), dpi=300, bbox_inches='tight')

# Phi =0.65
# ax2 = axes[1]
# %%
ds_plot = ds_alpha_fit.sel(run=('2023-05-18', 1)).sel(mp='mw_horns').sel(phi=0.8, method='nearest')

# da_plot.plot(col='phi', hue='var', row='motor')
motor_sel = [80, 130, 180]

ds_plot = ds_plot.sel(motor=motor_sel, method='nearest')

fig, axes = plt.subplots(1,len(motor_sel), figsize=(5,2), sharey=True)

for i, motor in enumerate(ds_plot.coords['motor']):
    ax = axes[i]
    ds_plot.sel(motor=motor).to_array('var').plot(hue='var', ax=ax)
    ax.set_title('{:.0f} mm'.format(motor))

    ax.set_xlabel("Wavelength [nm]")

    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(760,775)

    if i != 2:
        ax.get_legend().remove()

    else:
        ax.legend(['Fit', 'Fitted Data', 'Raw Data'])

axes[0].set_ylabel('$\\alpha$')


plt.savefig(pjoin(DIR_FIG_OUT, 'pos_alpha_fit.png'), dpi=300, bbox_inches='tight')

# %%


ds_plot