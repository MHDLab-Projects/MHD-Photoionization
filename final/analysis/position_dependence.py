
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

import scienceplots
import mhdpi_utils
plt.style.use(['science', 'ieee', 'mhdpi_utils.mystyle'])

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_p = ds_p.drop(34.81, 'motor')

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'modeling', 'cfd', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd = xr.load_dataset(fp_cfd)

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd = ds_cfd.assign_coords(x = ds_cfd.coords['x'].values - ds_cfd.coords['x'].values[0])
ds_cfd = ds_cfd.assign_coords(x = ds_cfd.coords['x'].values*1000)

ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd['nK_m3'] = ds_cfd['Yeq_K'].pint.to('1/m^3')



#%%


ds_p['nK_mw_horns'].mean('run').plot(hue='phi')

plt.yscale('log')

#%%

ds_p['AS_max'].mean('run').plot(hue='phi')

# ds_p
# %%
ds_p_sel = ds_p.sel(phi=0.6, method='nearest')
ds_cfd_sel = ds_cfd.sel(phi=0.6, method='nearest')

nK_barrel_mean = ds_p_sel['nK_barrel'].mean('motor').mean('run')
nK_barrel_std = ds_p_sel['nK_barrel'].mean('motor').std('run')

plt.figure(figsize=(4.5,2))

da = ds_p_sel['nK_mw_horns'].dropna('run', how='all')
da = da.sel(motor = slice(50,300))

g = da.isel(run=0).plot(x='motor', marker='o', label='2023-05-24 Run 1')
plot.dropna(g)

ax1 = plt.gca()
ax1.set_yscale('log')
ax1.set_title('')
ax1.set_xlabel('Position [mm]')

ax1.axvline(178, color='gold', linestyle='--')

ds_cfd_sel['nK_m3'].plot(color='black', label ='CFD centerline')

plt.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg')

plt.ylim(1e16,3e22)
plt.xlim(-10,250)


ax1.legend()
ax1.set_title('')

ax1.set_xlabel('Position [mm]')


plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mws_cfd.png'), dpi=300, bbox_inches='tight')


# %%
