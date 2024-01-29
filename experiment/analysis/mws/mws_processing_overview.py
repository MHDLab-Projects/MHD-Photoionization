#%%

#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

# %%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()
# %%

ds_sel = ds_lecroy.sel(run=('2023-05-24',1)).sel(kwt=1, method='nearest')

ds_sel = ds_sel.mean('mnum')

ds_sel

#%%

fig, axes = plt.subplots(3, figsize=(4,4), sharex=True)

ds_sel[['i','q']].to_array('var').plot(x='time', hue='var',ax=axes[0])
axes[0].set_ylabel('I and Q [V]')
axes[0].legend(['I','Q'])
axes[0].set_xlabel('')
ds_sel['mag'].plot(ax=axes[1])
axes[1].set_ylabel('Magnitude [V]')
axes[1].set_xlabel('')
ds_sel['AS'].plot(ax=axes[2])
axes[2].set_ylabel('1-T = AS')

for ax in axes:
    ax.set_title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_processing_overview.png'), bbox_inches='tight')


