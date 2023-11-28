#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

# fp = pjoin(REPO_DIR, 'experiment','data','munged','2023-05-24','Munged','Spectral','ds_absem_mp.cdf')
fp = pjoin(REPO_DIR, 'experiment','data','proc_data','ds_absem.cdf')

ds = xr.load_dataset(fp)

ds
# %%
# ds.isel(time=slice(0,50))['time']

ds.set_index(acq=['time','mp']).unstack('acq')