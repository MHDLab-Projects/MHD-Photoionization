#%%

from mhdlab.analysis.standard_import import *
from params import duration_lookup, time_downselect_lookup

plt.rcParams['font.size'] = 22

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_absem.cdf')

ds = xr.load_dataset(fp)

ds = ds.absem.calc_alpha()

ds = ds.set_index(acq=['time','mp']).unstack('acq').rename(time='acq_time')

# tweaks to display info
ds = ds.rename(alpha='Absorption')
da_plot = ds['Absorption'].pint.dequantify()
del da_plot.attrs['units']
da_plot.coords['wavelength'].attrs = {'units': 'nm', 'long_name': 'Wavelength'}

#%%
from mhdlab.plot.anim.absem import gen_movie_absem_mp

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

for date, tw in df_exptw.ct.iterslices():

    da_sel = da_plot.sel(acq_time=tw)

    time_downselect = time_downselect_lookup['absem']
    if time_downselect:
        da_sel = da_sel.isel(acq_time=slice(0,-1,time_downselect))
    da_sel = da_sel.sel(wavelength=slice(760,780))

    fp_out = 'output/absem_movie_{}.mp4'.format(date)
    duration = duration_lookup[date]

    gen_movie_absem_mp(da_sel, fp_out, tw, duration, mp_names=['Barrel', 'Jet'])

# %%
