#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds = xr.open_dataset(pjoin(data_directory, '536_pos_ds_p_stats.cdf')).pint.quantify().xr_utils.stack_run()
ds_lecroy = xr.open_dataset(pjoin(data_directory, '536_pos_lecroy.cdf')).pint.quantify().xr_utils.stack_run()

#%%

ds.drop('run_plot')['nK_m3']#.pint.quantify()

# ds

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
