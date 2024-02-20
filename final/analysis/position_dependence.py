
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_p = ds_p.drop(34.81, 'motor')

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd = xr.load_dataset(fp_cfd)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('particle/m^3')



#%%


ds_p['nK_mw_horns'].mean('run').plot(hue='phi')

plt.yscale('log')

#%%

ds_p['AS_max'].mean('run').plot(hue='phi')
#%%
motor_sel = [34.81, 104.8, 178, 226.7]

ds_lecroy = ds_lecroy.sel(phi=0.79)

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

