#%%

#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdpy.analysis import mws

# %%

ds_lecroy = ppu.fileio.load_lecroy('53x', norm_mag=True)

ds_sel = ds_lecroy.sel(run=('2023-05-18',1)).sel(kwt=1, method='nearest')

ds_sel = ds_sel.mean('mnum')

ds_sel

#%%

fig, axes = plt.subplots(4, figsize=(4,6), sharex=True)

ds_sel[['i','q']].to_array('var').plot(x='time', hue='var',ax=axes[0])
axes[0].set_ylabel('I and Q [V]')
axes[0].legend(['I','Q'])
axes[0].set_xlabel('')
ds_sel['mag'].plot(ax=axes[1])
axes[1].set_ylabel('U [V]')
axes[1].set_xlabel('')

mag_0 = ds_sel['mag_0'].item()
axes[1].axhline(mag_0, color='k', linestyle='--')

ds_sel['T_abs'].plot(ax=axes[2])
axes[2].set_ylabel('T_abs')

ds_sel['AS_abs'].plot(ax=axes[3])
axes[3].set_ylabel('AS_abs')

for ax in axes:
    ax.set_title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_processing_overview.png'), bbox_inches='tight')



# %%

ds_sel

ds_sel['mag'].plot(marker='o')

plt.xlim(-1,1)