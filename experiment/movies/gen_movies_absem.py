#%%

from mhdlab.analysis.standard_import import *
from params import duration_lookup, time_downselect_lookup

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_absem.cdf')

ds = xr.load_dataset(fp)

ds = ds.absem.calc_alpha()

ds = ds.set_index(acq=['time','mp']).unstack('acq').rename(time='acq_time')

ds

#%%
from mhdlab.plot.anim.absem import gen_movie_absem_mp

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

for date, tw in df_exptw.ct.iterslices():

    ds_sel = ds.sel(acq_time=tw)

    da_plot = ds_sel['alpha']
    time_downselect = time_downselect_lookup['absem']
    if time_downselect:
        da_plot = da_plot.isel(acq_time=slice(0,-1,time_downselect))
    da_plot = da_plot.sel(wavelength=slice(760,780))

    fp_out = 'output/absem_movie_{}.mp4'.format(date)
    duration = duration_lookup[date]

    gen_movie_absem_mp(da_plot, fp_out, tw, duration)

# %%
