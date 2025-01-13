#%%
from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu


# %%

fp = pjoin(REPO_DIR, 'experiment', 'data','munged', '2023-05-12', 'Lecroy', 'ds_savingtest_15hz.cdf')

ds = xr.load_dataset(fp)

ds.coords['time'].attrs['units'] = 'us'

ds = ds.mwt.calc_AS_rel()

ds['mag'].mean('acq_time').plot()

#%%

tc = '536_pos'

ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')
da_nothing = ppu.fileio.load_mwt_T0()

ds_lecroy = ds_lecroy.mwt.calc_time_stats()

# %%

da_sel = ds_lecroy['mag_pp']

# da_sel.plot(hue='run_plot', x='motor', marker='o')

norm = da_sel.sel(motor=slice(200,300)).mean('motor')

#%%

da_sel = ds_lecroy['mag'].sel(motor=[50,100,150,180,225], method='nearest')

da_sel = da_sel/norm

da_sel = da_sel.dropna('run', how='all')

g = da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))

dropna(g)

#%%

mag_pp = ds_lecroy['mag_pp'].dropna('run','all')

g = mag_pp.plot(hue='run_plot', x='motor')

dropna(g)

#%%

da_mag_norm = mag_pp/ds_lecroy['mag_0']

g = da_mag_norm = da_mag_norm.plot(hue='run_plot', x='motor')

dropna(g)
#%%[markdown]

# %%

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
# Display img on the upper left subplot

da_nothing.to_series().plot(ax=axes[0], marker='o')

axes[0].set_title('T No Torch')

axes[0].set_ylabel("$U_{Nothing} (V)$")
plt.xlabel("Date")

da_sel = ds_lecroy['mag_pp']
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()
g = da_sel.plot(hue='run_plot', x='motor', marker='o', ax=axes[1])

dropna(g)

da_sel = ds_lecroy['mag'].sel(motor=[50,100,150,180,225], method='nearest')
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()

axes[1].set_title('$T_{abs} = U/U_{Nothing}$')
axes[1].set_ylabel('$T_{abs}$')
axes[1].set_xlabel("Motor Position (mm)")

plt.savefig(pjoin(DIR_FIG_OUT, 'mwt_nothing_T0.png'))

# %%
