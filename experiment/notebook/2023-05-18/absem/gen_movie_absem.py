#%%

from mhdpy.analysis.standard_import import *

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_absem.cdf')

ds = xr.load_dataset(fp)

ds = ds.absem.calc_alpha()

ds = ds.set_index(acq=['time','mp']).unstack('acq').rename(time='acq_time')

ds

#%%

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

tw = df_exptw.ct.index_slice('2023-05-18')
ds = ds.sel(acq_time=tw)

#%%
da_plot = ds['alpha'].isel(acq_time=slice(0,-1,10))
da_plot = da_plot.sel(wavelength=slice(760,780))

fp_out = 'output/absem_movie.mp4'
duration = 132

from mhdpy.plot.anim.absem import gen_movie_absem_mp

gen_movie_absem_mp(da_plot, fp_out, tw, duration)
# %%