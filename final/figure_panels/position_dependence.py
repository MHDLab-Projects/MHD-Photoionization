#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdpy.plot import xr_errorbar_axes

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

# Position dependent AES and MWS
ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_alpha_fit = xr.open_dataset(pjoin(data_directory, 'ds_pos_alpha_fit.cdf')).xr_utils.stack_run()
# Non position dependent barrel AES
ds_p_barrel = xr.open_dataset(pjoin(data_directory, 'ds_p_barrel.cdf')).xr_utils.stack_run()

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_lecroy = ds_lecroy.sel(phi=0.8, method='nearest')

ds_lecroy_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')
ds_536 = ds_lecroy_536.mws.calc_time_stats()

#%%

ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd['nK_m3'] = ds_cfd[ppu.CFD_K_SPECIES_NAME].pint.to('particle/m^3')


#%%

fig, axes = plt.subplots(2, figsize=(5,5), sharex=True)


# Phi =0.8
ax1 = axes[0]

ds_p_sel = ds_p.sel(phi=0.8, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.8, method='nearest')

phi_val_expt = ds_p_sel.coords['phi'].item()


nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: removing last datapoint. Can this be fixed with calibraiton/processing?
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,220))

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-18 Run 1', ax=ax1)
g = nK_mw_horns.isel(run=1).plot(x='motor', marker='o', label='2023-05-18 Run 2', ax=ax1)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax1)
# ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax1)


# Barrel exit
nK_barrel_mean = ds_p_barrel['nK_barrel_cfdprofile'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_barrel['nK_barrel_cfdprofile'].std('run').sel(phi=0.8, method='nearest')
ax1.errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='o', label='Barrel Exit AES', capsize=5, )

ax1.set_title('Equivalence Ratio = {}'.format(phi_val_expt))


# Phi = 0.65

ax2 = axes[1]

ds_p_sel = ds_p.sel(phi=0.6, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.6, method='nearest')

phi_val_expt = ds_p_sel.coords['phi'].item()


nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: final datapoint for 516 has no potassium peaks. However, calibration is off and fit gives a large absorption value. Removing is probably best as any absorption is not meaninful anyway. 
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,200)) 

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-24 Run 1', ax=ax2)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax2)
# ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax2)

nK_barrel_mean = ds_p_barrel['nK_barrel_cfdprofile'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_barrel['nK_barrel_cfdprofile'].std('run').sel(phi=0.8, method='nearest')
ax2.errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='o', label='Barrel Exit AES', capsize=5, )


ax2.set_title('Equivalence Ratio = {}'.format(phi_val_expt))

for ax in axes:

    ax.axvline(180, color='gold', linestyle='--')
    ax.set_ylim(1e16,1e23)
    ax.set_xlim(-10,220)
    ax.set_yscale('log')
    ax.legend()
    ax.set_ylabel('$n_K [m^{-3}$]')

ax2.set_xlabel('Stage Position [mm]')
ax1.set_xlabel('')

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mws_cfd.png'), dpi=300, bbox_inches='tight')

# Phi =0.65
# ax2 = axes[1]

#%%

plt.figure()

ds_plot = ds_alpha_fit.sel(run=('2023-05-18', 1)).sel(phi=0.8, method='nearest')

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

#%%

plt.figure()

ds_plot = ds_536.sel(run=('2023-05-18', 1)).dropna('motor', how='all')

da_plot = ds_plot['AS_abs']

motor_sel = [30, 100, 180]
motor_actual = []

fig, ax = plt.subplots(figsize=(3,3))

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Define a list of colors for the traces

for i, motor in enumerate(motor_sel):
    # Calculate mean and standard deviation
    mean = da_plot.sel(motor=motor,method='nearest').mean(dim='mnum', keep_attrs=True)
    std = da_plot.sel(motor=motor,method='nearest').std(dim='mnum', keep_attrs=True)

    motor_val = da_plot.sel(motor=motor, method='nearest').coords['motor'].values
    motor_actual.append(motor_val)
    
    # Plot mean and shaded region for standard deviation
    ax.plot(mean['time'].pint.magnitude, mean.pint.magnitude, color=colors[i % len(colors)], label='{:.0f} mm'.format(motor_val))
    ax.fill_between(mean['time'].pint.magnitude, (mean-std).pint.magnitude, (mean+std).pint.magnitude, color=colors[i % len(colors)], alpha=0.2)

ax.set_xlabel('Time [us]')
ax.set_ylabel('AS [dimensionless]')
ax.legend(bbox_to_anchor=(0.5,0.87))
ax.set_ylim(-0.05,1.05)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_time_20230518.png'), dpi=300, bbox_inches='tight')
# %%

plt.figure()


da_plot = ds_536['AS_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all')

AS_pp = da_plot.sel(time=slice(-50,-1)).mean('mnum').mean('time').pint.dequantify()
AS_pp_fluct = da_plot.sel(time=slice(-50,-1)).std('mnum').mean('time').pint.dequantify()
AS_laser = da_plot.mws._pulse_max().mean('mnum').pint.dequantify()
AS_laser_fluct = da_plot.mws._pulse_max().std('mnum').pint.dequantify()

fig, ax = plt.subplots()

xr_errorbar_axes(AS_pp, AS_pp_fluct, axes=ax, capsize=5, label='AS Pre Laser')
xr_errorbar_axes(AS_laser, AS_laser_fluct, axes=ax, capsize=5, label='AS Laser')

for motor_val in motor_actual:
    plt.axvline(motor_val, color='k', linestyle='--', alpha=0.5)

dropna(ax)

plt.legend()
plt.xlim(-5,270)

plt.ylabel('AS [dimensionless]')
plt.xlabel('Stage Position [mm]')

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_stats_20230518.png'))

#%%
#%%

plt.figure()

da_plot = ds_536['SFR_abs'].dropna('run', how='all')

g=  da_plot.mean('mnum', keep_attrs=True).plot(hue='run_plot', x='motor', marker='o')

for motor_val in motor_actual:
    plt.axvline(motor_val, color='k', linestyle='--', alpha=0.5)

plt.title('')
plt.xlim(-5, 270)
plt.ylim(-5,125)
plt.xlabel('Stage Position [mm]')
plt.ylabel('SFR [dimensionless]')

dropna(g)

plt.savefig(pjoin(DIR_FIG_OUT, 'pow_mws_sfr_dates.png'))


# %%
