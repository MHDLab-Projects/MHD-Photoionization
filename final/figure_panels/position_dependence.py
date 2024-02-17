#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

import scienceplots
import mhdpi_utils
plt.style.use(['science', 'ieee', 'mhdpi_utils.mystyle'])

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds = xr.open_dataset(pjoin(data_directory, '536_pos_ds_p_stats.cdf')).pint.quantify().xr_utils.stack_run()
ds_lecroy = xr.open_dataset(pjoin(data_directory, '536_pos_lecroy.cdf')).pint.quantify().xr_utils.stack_run()


#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'modeling', 'cfd', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd = xr.load_dataset(fp_cfd)

ds_cfd = ds_cfd.sel(kwt=1).sel(phi=0.8)

ds_cfd = ds_cfd.assign_coords(x = ds_cfd.coords['x'].values - ds_cfd.coords['x'].values[0])
ds_cfd = ds_cfd.assign_coords(x = ds_cfd.coords['x'].values*1000)

ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('1/m^3')


#%%

motor_sel = [34.81, 104.8, 178, 226.7]

da_sel = ds_lecroy['AS'].mean('mnum').sel(motor=motor_sel, method='nearest')

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

#%%

fig, axes = plt.subplots(3, figsize=(3,10), sharex=True)

for i, var in enumerate(['AS_max', 'mag_pp_std', 'AS_max_std_ratio']):
    ax = axes[i]

    ds_stat = ds.wma.initialize_stat_dataset(var, 'run')

    # plot with confidence interval

    # ds_stat['mean'].plot(ax=ax)
    # ax.fill_between(ds_stat.coords['motor'], ds_stat['mean'] - ds_stat['std'], ds_stat['mean'] + ds_stat['std'], alpha=0.2)

    ax.errorbar(ds_stat.coords['motor'], ds_stat['mean'], yerr=ds_stat['std'], marker='o', capsize=5)

    ax.set_title('')
    ax.set_xlabel('')

    unit_str = '[' + ds[var].attrs['units']+ ']' if 'units' in ds[var].attrs else '' 
    ax.set_ylabel(ds[var].attrs['long_name'] + unit_str)
    
axes[2].set_xlabel('Position [mm]')

ta = axes[0].twinx()

delta_pd1 = ds['delta_pd1'].dropna('run', how='all')
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

nK_barrel_mean = ds['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds['nK_barrel'].mean('motor').std('run')

plt.figure(figsize=(4.5,2))

da = ds['nK_mw_horns'].dropna('run', how='all')
da = da.sel(motor = slice(50,300))

g = da.isel(run=0).plot(x='motor', marker='o', label='2023-05-12 Run 1')
g = da.isel(run=1).plot(x='motor', marker='o', label='2023-05-12 Run 2')
plot.dropna(g)

ax1 = plt.gca()
ax1.set_yscale('log')
ax1.set_title('')
ax1.set_xlabel('Position [mm]')

# ax1.get_legend().set_title('Experiment (date, #)')
# ax1.get_legend().set_bbox_to_anchor([0,0,0.5,0.4])
ax1.axvline(178, color='gold', linestyle='--')


ds_cfd['nK_m3'].plot(color='black', label ='CFD centerline')

plt.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg')

plt.ylim(1e16,3e22)
plt.xlim(-10,250)



ax1.legend()
ax1.set_title('')

ax1.set_xlabel('Position [mm]')

# ax2 = plt.gca()
# ax2.legend()
# ax2.set_ylabel('Normalized $n_{K, CFD}$')
# ax2.set_yscale('log')
# ax2.set_title('')

# plt.xlim(0,310)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mws_cfd.png'), dpi=300, bbox_inches='tight')


