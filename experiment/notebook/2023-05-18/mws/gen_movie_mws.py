#%%

from mhdpy.analysis.standard_import import *

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_lecroy.cdf')

ds = xr.load_dataset(fp)
ds = ds[['i','q']]

ds = ds.mws.calc_mag_phase()

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

tw = df_exptw.ct.index_slice('2023-05-18')
ds = ds.sel(acq_time=tw)

#%%

da_plot = ds['mag'].isel(acq_time=slice(0,-1,100))
da_plot = da_plot.sel(time=slice(-1,30))

duration = 132
fp_out = 'output/mws_movie.mp4'

from mhdpy.plot.anim.mws import gen_movie_mws

gen_movie_mws(da_plot, fp_out, tw, duration)
