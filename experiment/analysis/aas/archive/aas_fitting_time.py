#%%

from mhdlab.analysis.standard_import import *
from mhdlab.plot import dropna
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdlab.analysis import aas

datestr = '2023-05-18'

# datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
data_folder = pjoin(REPO_DIR, 'experiment', 'data','munged',datestr)

dsst = mhdlab.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

ds_aas = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'ds_aas.cdf'))

ds_aas = ds_aas.aas.calc_alpha()

ds_aas
#%%

ds_fit = ds_aas.rename(acq='time').set_index(time='time')

# ds_fit = ds_fit.sel(time=ds_fit.coords['time'].values[::100])

ds_fit = ds_fit.aas.calc_alpha()

#Downselect to remove data where there are no peaks... TODO: revisit
ds_fit = ds_fit.where(ds_fit['alpha'].sel(wavelength=slice(760,775)).max('wavelength') > 0.5).dropna('time', how='all')

spectral_reduction_params_fp = os.path.join(REPO_DIR, 'experiment', 'metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_alpha_fit, ds_p, ds_p_stderr= aas.fitting.pipe_fit_alpha_1(ds_fit, fit_prep_kwargs={'spect_red_dict': spect_red_dict})

ds_alpha_fit['nK_m3'] = ds_p['nK_m3']

#%%

ds = ds_alpha_fit[['alpha', 'alpha_fit', 'nK_m3']].sel(wavelength=slice(760,775))
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})

ds = ds.set_index(acq=['time','mp']).unstack('acq')

ds
#%%

tw = slice(Timestamp(datestr), Timestamp(datestr) + pd.Timedelta('1 day'))

g = ds['nK_m3'].sel(time=tw).plot(hue='mp')
dropna(g)
plt.yscale('log')

#%%

# Assign coordinates to ds
#TODO: this is redundant with proc_add_coord.py, but for simlified case of individual tc. how should this be handled?

from mhdlab.fileio.tdms import tdms2ds
from mhdlab.fileio.ct import load_df_cuttimes
from mhdlab.coords import assign_coords_multi

coords_to_assign = tdms2ds(pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms'))
coords_to_assign
# %%

# Cuttimes 
df_cuttimes = load_df_cuttimes(pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_sequence.csv'))

# add tc and run_num columns. tc is everything before the last underscore in 'Event' column
df_cuttimes['tc'] = df_cuttimes['Event'].str.extract(r'(.*)_\d$')
df_cuttimes['run_num'] = df_cuttimes['Event'].str.extract(r'.*_(\d)$')

process_tcs = [
    '53x',
    ]

# Only keep columns with tc in process_tcs
df_cuttimes = df_cuttimes[df_cuttimes['tc'].isin(process_tcs)]

df_cuttimes = df_cuttimes.sort_values('Start Time').set_index(['tc', 'run_num' ,'date'])

#%%
dss = []

for (tc_base, run_num, date), row in df_cuttimes.iterrows():

    sl = slice(
        np.datetime64(row['Start Time']), 
        np.datetime64(row['Stop Time'])
        )

    ds_sel = ds.sel(time=sl)



    ds_c = assign_coords_multi(ds_sel.rename(time='acq_time'), coords_to_assign, min_mnum=1, time_coord='acq_time')

    #TODO: believe these ohter coordinates are showing up because of low min_mnum
    if tc_base == '53x':
        if 'phi' in ds_c.dims:
            ds_c = ds_c.sel(phi=0.8, method='nearest')

    ds_c = ds_c.assign_coords(date=[date])
    ds_c = ds_c.assign_coords(run_num=[int(run_num)])#.expand_dims('run_num')
    ds_c = ds_c.xr_utils.stack_run()

    dss.append(ds_c)

ds_c = xr.concat(dss, dim='run')

# Assign run_plot coordinate
ds_c = ds_c.assign_coords(run_plot = ('run', ds_c.indexes['run'].values))

# %%


g = ds_c['nK_m3'].mean('mnum').plot(col='mp', x='kwt',hue='run_plot',marker='o')

dropna(g)

plt.yscale('log')
plt.xscale('log')
# %%

# da = ds_c['data'].sel(mp='barrel')

ds_plot = ds_c[['data','fit']].sel(mp='barrel').dropna('kwt', how='all').mean('mnum')
da_plot = ds_plot.to_array('var')

da_plot

#%%

g = da_plot.plot(col='kwt',hue='var', row='run', ylim = (1e-3,2), yscale='log', figsize=(10,6))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')


# %%
