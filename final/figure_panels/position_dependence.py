#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils #Sets matplotlib style


data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

ds_lecroy = ds_lecroy.sel(phi=0.79)

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd = xr.load_dataset(fp_cfd)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('particle/m^3')


#%%

motor_sel = [34.81, 104.8, 178, 226.7]

da_sel = ds_lecroy['AS'].mean('mnum').sel(motor=motor_sel, method='nearest')

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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel.png'), dpi=300, bbox_inches='tight')


#%%

ds_p_sel = ds_p.sel(phi=0.79, method='nearest') 

fig, axes = plt.subplots(3, figsize=(3,10), sharex=True)

for i, var in enumerate(['AS_max', 'mag_pp_std', 'AS_max_std_ratio']):
    ax = axes[i]

    ds_stat = ds_p_sel.wma.initialize_stat_dataset(var, 'run')

    # plot with confidence interval

    # ds_stat['mean'].plot(ax=ax)
    # ax.fill_between(ds_stat.coords['motor'], ds_stat['mean'] - ds_stat['std'], ds_stat['mean'] + ds_stat['std'], alpha=0.2)

    ax.errorbar(ds_stat.coords['motor'], ds_stat['mean'], yerr=ds_stat['std'], marker='o', capsize=5)

    ax.set_title('')
    ax.set_xlabel('')

    unit_str = '[' + ds_p_sel[var].attrs['units']+ ']' if 'units' in ds_p_sel[var].attrs else '' 
    ax.set_ylabel(ds_p_sel[var].attrs['long_name'] + unit_str)
    
axes[2].set_xlabel('Position [mm]')

ta = axes[0].twinx()

delta_pd1 = ds_p_sel['delta_pd1'].dropna('run', how='all')
delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')
delta_pd1.plot(ax=ta, color='red', marker='o')
ta.set_ylabel('Delta PD1 [mV]')
ta.set_title('')

for mp in motor_sel:
    if mp == 178:
        axes[2].axvline(mp, color='gold', linestyle='--')
    else:
        axes[2].axvline(mp, color='gray', linestyle='--')


plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_stats.png'), dpi=300, bbox_inches='tight')

#%%

fig, axes = plt.subplots(2, figsize=(6,6), sharex=True)


# Phi =0.8
ax1 = axes[0]

ds_p_sel = ds_p.sel(phi=0.79, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.79, method='nearest')

nK_barrel_mean = ds_p_sel['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds_p_sel['nK_barrel'].mean('motor').std('run')

nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all')

#TODO: removing last datapoint. Can this be fixed with calibraiton/processing?
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,220))

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-18 Run 1', ax=ax1)
g = nK_mw_horns.isel(run=1).plot(x='motor', marker='o', label='2023-05-18 Run 2', ax=ax1)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax1)
ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm off centerline', linestyle='--', ax=ax1)

ax1.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg', capsize=5, )
ax1.set_title('phi = 0.79')


# Phi = 0.65

ax2 = axes[1]

ds_p_sel = ds_p.sel(phi=0.65, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.65, method='nearest')

nK_barrel_mean = ds_p_sel['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds_p_sel['nK_barrel'].mean('motor').std('run')

nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all')

#TODO: final datapoint for 516 has no potassium peaks. However, calibration is off and fit gives a large absorption value. Removing is probably best as any absorption is not meaninful anyway. 
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,200)) 

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-24 Run 1', ax=ax2)

ds_cfd_sel['nK_m3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax2)
ds_cfd_sel['nK_m3'].sel(offset=3).plot(color='black', label ='CFD 3mm off centerline', linestyle='--', ax=ax2)

ax2.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg', capsize=5, )
ax2.set_title('phi = 0.65')

for ax in axes:

    ax.axvline(178, color='gold', linestyle='--')
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
