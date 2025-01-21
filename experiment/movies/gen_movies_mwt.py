#%%

from mhdlab.analysis.standard_import import *
from params import duration_lookup, time_downselect_lookup

plt.rcParams['font.size'] = 22

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_lecroy.cdf')

ds = xr.load_dataset(fp)
ds = ds[['i','q']]

ds = ds.mwt.calc_mag_phase()

ds.coords['time'].attrs = {'units': r'$\mu s$', 'long_name': 'Time'}

ds['mag'] = ds['mag'].pint.dequantify()
ds['mag'].attrs = {'units': 'V', 'long_name': 'Microwave Transmission'}

#%%

from mhdlab.plot.anim.mwt import gen_movie_mwt

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')


da_plot = ds['mag']

for date, tw in df_exptw.ct.iterslices():

    da_sel = da_plot.sel(acq_time=tw)

    time_downselect = time_downselect_lookup['mwt']
    if time_downselect:
        da_sel = da_sel.isel(acq_time=slice(0,-1,time_downselect))
    da_sel = da_sel.sel(time=slice(-1,30))

    fp_out = 'output/mwt_movie_{}.mp4'.format(date)
    duration = duration_lookup[date]

    gen_movie_mwt(da_sel, fp_out, tw, duration)

# %%

ds
