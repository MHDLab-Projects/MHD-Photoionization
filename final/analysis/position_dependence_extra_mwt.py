
#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdlab.plot import xr_errorbar_axes

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
# ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
# ds_p = ds_p.drop(34.81, 'motor')

ds_lecroy_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')
ds_536 = ds_lecroy_536.mwt.calc_time_stats()

#%%

da_sel = ds_lecroy.sel(phi=0.8, method='nearest')['dAS_abs'].dropna('run',how='all')

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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_sel.png'), dpi=300, bbox_inches='tight')


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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_stats.png'), dpi=300, bbox_inches='tight')


#%%



da_plot = ds_536['T_abs'].mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,180,200,250], method='nearest')

da_plot.plot(hue='motor')

#%%
da_plot = ds_536['AS_abs'].mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,180,250], method='nearest')

da_plot.plot(hue='motor')

plt.ylabel('AS')



# %%



#%%

da_plot = ds_536['AS_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all').mean('mnum')
g = da_plot.sel(motor=[100,180], method='nearest').plot(col='motor', sharey=False)

for ax in g.axes.flatten():
    # ax.set_yscale('log')
    ax.set_ylabel('AS')

#%%
    
da_plot

#%%

AS_pp = ds_536['AS_abs'].sel(time=slice(-50,-1)).mean('time')
AS_pp.mean('mnum').mean('run').plot(label='AS Pre Pulse')

AS_laser = ds_536['AS_abs'].mwt._pulse_max()
AS_laser.mean('mnum').mean('run').plot(label='AS Laser')

plt.ylabel('AS')

plt.legend()


#%%

# Pre Pulse
AS_pp_mean = AS_pp.mean('mnum').mean('run')
AS_pp_std = AS_pp.mean('mnum').std('run')

# Laser
AS_laser_mean = AS_laser.mean('mnum').mean('run')
AS_laser_std = AS_laser.mean('mnum').std('run')

fig, ax = plt.subplots()


xr_errorbar_axes(AS_pp_mean, AS_pp_std, label='AS Pre Pulse', axes=ax)
xr_errorbar_axes(AS_laser_mean, AS_laser_std, label='AS Laser', axes=ax)

plt.ylabel('AS')

#%%

count = ds_536['T_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all').count('mnum').isel(time=0)

count.plot()

plt.ylabel("Count")





#%%

ds_p['nK_mw_horns'].mean('run').plot(hue='phi')

plt.yscale('log')

#%%

ds_p['dAS_abs_max'].mean('run').plot(hue='phi')

#%%

ds_lecroy
#%%
motor_sel = [34.81, 104.8, 178, 226.7]

ds_sel = ds_lecroy.sel(phi=0.8, method='nearest')

da_sel = ds_sel['AS_abs'].sel(motor=motor_sel, method='nearest')

fig, axes = plt.subplots(4, figsize=(3,12), sharex=True, sharey=True)

for i, motor in enumerate(da_sel.coords['motor'].values):
    ax = axes[i]
    da_sel.sel(motor=motor).plot(hue='run_plot', x='time', ax=ax)
    ax.set_title('Position: {} mm'.format(motor))
    ax.set_xlabel('')
    ax.set_ylabel('AS')

    if i == len(da_sel.coords['motor'].values) - 1:
        ax.set_xlabel('Time [us]')

    else:
        ax.get_legend().remove()

# da_sel.plot(row='motor', hue='run_plot', x='time', figsize=(8,25))

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_sel.png'), dpi=300, bbox_inches='tight')


# %%
motor_sel = [34.81, 104.8, 178, 226.7]

da_sel = ds_lecroy['AS_abs'].sel(motor=motor_sel, method='nearest')

da_sel = da_sel.sel(phi=0.8, method='nearest')

da_sel_mean = da_sel.mean('run')
da_sel_std = da_sel.std('run')

fig, axes = plt.subplots(4, figsize=(3,12), sharex=True, sharey=True)

for i, motor in enumerate(da_sel.coords['motor'].values):
    ax = axes[i]
    
    # Calculate mean and standard deviation
    mean = da_sel.sel(motor=motor).mean(dim='run')
    std = da_sel.sel(motor=motor).std(dim='run')
    
    # Plot mean and shaded region for standard deviation
    mean.plot(x='time', ax=ax)
    ax.fill_between(mean['time'].values, (mean-std).values, (mean+std).values, color='b', alpha=0.2)
    
    ax.set_title('Position: {} mm'.format(motor))
    ax.set_xlabel('')
    ax.set_ylabel('AS')

    if i == len(da_sel.coords['motor'].values) - 1:
        ax.set_xlabel('Time [us]')
    # else:
    #     ax.get_legend().remove()

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_sel_row.png'), dpi=300, bbox_inches='tight')

#%%

motor_sel = [34.81, 104.8, 178]

fig, ax = plt.subplots(figsize=(3,3))

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Define a list of colors for the traces

for i, motor in enumerate(motor_sel):
    # Calculate mean and standard deviation
    mean = da_sel.sel(motor=motor,method='nearest').mean(dim='run')
    std = da_sel.sel(motor=motor,method='nearest').std(dim='run')

    motor_val = da_sel.sel(motor=motor, method='nearest').coords['motor'].values
    
    # Plot mean and shaded region for standard deviation
    ax.plot(mean['time'].values, mean.values, color=colors[i % len(colors)], label='{:.0f} mm'.format(motor_val))
    ax.fill_between(mean['time'].values, (mean-std).values, (mean+std).values, color=colors[i % len(colors)], alpha=0.2)

ax.set_xlabel('Time [us]')
ax.set_ylabel('AS')
ax.legend(bbox_to_anchor=(1,1.1))
ax.set_ylim(0,1)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_sel.png'), dpi=300, bbox_inches='tight')
# %%
