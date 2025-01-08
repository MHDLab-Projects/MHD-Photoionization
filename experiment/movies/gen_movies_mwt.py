#%%

from mhdpy.analysis.standard_import import *
from params import duration_lookup, time_downselect_lookup

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_lecroy.cdf')

ds = xr.load_dataset(fp)
ds = ds[['i','q']]

ds = ds.mwt.calc_mag_phase()

from mhdpy.plot.anim.mwt import gen_movie_mwt

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

for date, tw in df_exptw.ct.iterslices():

    ds_plot = ds.sel(acq_time=tw)

    da_plot = ds_plot['mag']
    time_downselect = time_downselect_lookup['mwt']
    if time_downselect:
        da_plot = da_plot.isel(acq_time=slice(0,-1,time_downselect))
    da_plot = da_plot.sel(time=slice(-1,30))

    fp_out = 'output/mwt_movie_{}.mp4'.format(date)
    duration = duration_lookup[date]

    gen_movie_mwt(da_plot, fp_out, tw, duration)

# %%

ds
