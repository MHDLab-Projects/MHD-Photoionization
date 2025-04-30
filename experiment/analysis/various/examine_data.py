#%%
from mhdlab.analysis.standard_import import *

fp = pjoin(REPO_DIR, 'experiment','data','munged','2023-05-24','Munged','Spectral','ds_aas_mp.cdf')
# fp = pjoin(REPO_DIR, 'experiment','data','proc_data','ds_aas.cdf')
# fp = pjoin(REPO_DIR, 'experiment','data','proc_data','aas', '53x.cdf')

ds = xr.load_dataset(fp)

ds

#%%

set(ds['mp'].values)

# %%

# ds['calib'].sel(date='2023-05-24', run_num=1, kwt=0.99, mp='barrel', mnum=1)

#%%

ds['calib'].isel(acq=-20)