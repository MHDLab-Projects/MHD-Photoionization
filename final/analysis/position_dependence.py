
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

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

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
