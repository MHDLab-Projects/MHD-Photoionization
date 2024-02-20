#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()



from mhdpy.pyvista_utils import CFDDatasetAccessor

fp = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)

ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)
ds_cfd = ds_cfd.sel(kwt=1)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd

#%%

goldi_pos =  Quantity(178, 'mm')
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')

#%%
# Calculate kr

ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()

from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted

ds_krb = gen_ds_krb(ds_cfd['T'], ds_cfd['rho_number'])

ds_krb['O2_C'] = calc_krbO2_weighted(ds_cfd)

ds_krb


# %%

ds_krb['']