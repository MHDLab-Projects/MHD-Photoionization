#%%
from mhdpy.analysis.standard_import import *
from mhdpy.coords import gen_coords_to_assign_1, assign_coords_multi
import pi_paper_utils as ppu


# %%

fp = r'/home/leeaspitarte/code/MHD-Photoionization/experiment/data/munged/2023-05-12/Lecroy/ds_savingtest_15hz.cdf'

ds = xr.load_dataset(fp)

ds.coords['time'].attrs['units'] = 'us'

ds = ds.mws.calc_mag_phase_AS()

ds['mag'].mean('acq_time').plot()

#%%

tc = '536_pos'

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]


#%%

fp_nothing = pjoin(REPO_DIR,'final', 'dataset', 'output', 'mws_T0.csv')

df_nothing = pd.read_csv(fp_nothing, index_col=0)['0']
da_nothing = xr.DataArray(df_nothing).pint.quantify('volt')

da_nothing

#%%

ds_lecroy = ds_lecroy.unstack('run').mws.calc_mag_phase_AS(mag_0=da_nothing)#[['mag', 'phase','AS']]

ds_lecroy = ds_lecroy.xr_utils.stack_run()  

#%%
ds_lecroy


# %%

da_sel = ds_lecroy['mag_pp'].mean('mnum')

# da_sel.plot(hue='run_plot', x='motor', marker='o')

norm = da_sel.sel(motor=slice(200,300)).mean('motor')

#%%

da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')

da_sel = da_sel/norm

da_sel = da_sel.dropna('run', how='all')

da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))

#%%

mag_pp = ds_lecroy['mag_pp'].mean('mnum').dropna('run','all')

mag_pp.plot(hue='run_plot', x='motor')

#%%

da_mag_norm = mag_pp/ds_lecroy['mag_0']

da_mag_norm = da_mag_norm.plot(hue='run_plot', x='motor')
#%%[markdown]

# %%

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
# Display img on the upper left subplot

da_nothing.to_series().plot(ax=axes[0], marker='o')

axes[0].set_title('T No Torch')

plt.ylabel("$U_{Nothing} (V)$")
plt.xlabel("Date")

da_sel = ds_lecroy['mag_pp'].mean('mnum')
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()
da_sel.plot(hue='run_plot', x='motor', marker='o', ax=axes[1])

da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()

axes[1].set_title('$U/U_{Nothing}$')
axes[1].set_xlabel("Motor Position (mm)")

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_nothing_T0.png'))

# %%
