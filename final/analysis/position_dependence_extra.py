
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdpy.plot.common import xr_errorbar_axes

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
# ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
# ds_p = ds_p.drop(34.81, 'motor')

ds_lecroy_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')

#%%


#%%

ds_536 = ds_lecroy_536.mws.calc_time_stats()


da_plot = ds_536['T_abs'].mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,180,200,250], method='nearest')

da_plot.plot(hue='motor')

#%%
da_plot = ds_536['AS_abs'].mean('mnum').mean('run')
da_plot = da_plot.sel(motor=[0,50,100,150,180,250], method='nearest')

da_plot.plot(hue='motor')

plt.ylabel('AS')


#%%

ds_plot = ds_536.sel(run=('2023-05-18', 1)).dropna('motor', how='all')

da_plot = ds_plot['AS_abs']

motor_sel = [30, 100, 180]
motor_actual = []

fig, ax = plt.subplots(figsize=(3,3))

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Define a list of colors for the traces

for i, motor in enumerate(motor_sel):
    # Calculate mean and standard deviation
    mean = da_plot.sel(motor=motor,method='nearest').mean(dim='mnum')
    std = da_plot.sel(motor=motor,method='nearest').std(dim='mnum')

    motor_val = da_plot.sel(motor=motor, method='nearest').coords['motor'].values
    motor_actual.append(motor_val)
    
    # Plot mean and shaded region for standard deviation
    ax.plot(mean['time'].values, mean.values, color=colors[i % len(colors)], label='{:.0f} mm'.format(motor_val))
    ax.fill_between(mean['time'].values, (mean-std).values, (mean+std).values, color=colors[i % len(colors)], alpha=0.2)

ax.set_xlabel('Time [us]')
ax.set_ylabel('AS')
ax.legend(bbox_to_anchor=(0.5,0.8))
ax.set_ylim(-0.05,1.05)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_time_20230518.png'), dpi=300, bbox_inches='tight')
# %%

da_plot = ds_536['AS_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all')

AS_pp = da_plot.sel(time=slice(-50,-1)).mean('mnum').mean('time')
AS_pp_fluct = da_plot.sel(time=slice(-50,-1)).std('mnum').mean('time')
AS_laser = da_plot.mws._pulse_max().mean('mnum')
AS_laser_fluct = da_plot.mws._pulse_max().std('mnum')

fig, ax = plt.subplots()

xr_errorbar_axes(AS_pp, AS_pp_fluct, axes=ax, capsize=5, label='AS Pre Laser')
xr_errorbar_axes(AS_laser, AS_laser_fluct, axes=ax, capsize=5, label='AS Laser')

for motor_val in motor_actual:
    plt.axvline(motor_val, color='k', linestyle='--', alpha=0.5)

dropna(ax)

plt.legend()

plt.ylabel('AS')

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_stats_20230518.png'))

#%%

da_plot = ds_536['SFR_abs'].dropna('run', how='all')

g=  da_plot.mean('mnum').plot(hue='run_plot', x='motor', marker='o')

for motor_val in motor_actual:
    plt.axvline(motor_val, color='k', linestyle='--', alpha=0.5)

plt.title('')
plt.xlim(0, )

dropna(g)

plt.savefig(pjoin(DIR_FIG_OUT, 'pow_mws_sfr_dates.png'))

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

AS_laser = ds_536['AS_abs'].mws._pulse_max()
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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel.png'), dpi=300, bbox_inches='tight')


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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel_row.png'), dpi=300, bbox_inches='tight')

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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel.png'), dpi=300, bbox_inches='tight')
# %%
