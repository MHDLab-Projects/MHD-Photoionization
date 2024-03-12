
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
# ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
# ds_p = ds_p.drop(34.81, 'motor')

#%%


ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(offset=0).sel(phi=0.8)

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('particle/m^3')

#%%


all_K_species = ['Yeq_K','Yeq_K+','Yeq_K2CO3','Yeq_KO','Yeq_KOH']
ds_cfd[[*all_K_species, 'all_K_Yeq']].to_array('var').plot(hue='var',row='kwt')

plt.yscale('log')

# plt.gca().get_legend().set_bbox_to_anchor((1,1))

plt.ylim(1e8,1e17)

plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_K_species.png'), dpi=300, bbox_inches='tight')


#%%



#%%


ds_p['nK_mw_horns'].mean('run').plot(hue='phi')

plt.yscale('log')

#%%

ds_p['AS_max'].mean('run').plot(hue='phi')

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

da_sel = ds_lecroy['AS'].sel(motor=motor_sel, method='nearest')

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

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mws_AS_sel.png'), dpi=300, bbox_inches='tight')


