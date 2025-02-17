#%%

#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdlab.analysis import mwt

# %%

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute', avg_mnum=True)

ds_sel = ds_lecroy.sel(run=('2023-05-18',1)).sel(kwt=1, method='nearest')

ds_sel = ds_sel.pint.dequantify()

# ds_sel = ds_sel.mean('mnum')

ds_sel

#%%

fig, axes = plt.subplots(4, figsize=(4,6), sharex=True)

ds_sel[['i','q']].to_array('var').plot(x='time', hue='var',ax=axes[0])
axes[0].set_ylabel('I and Q [V]')
axes[0].legend(['I(t)','Q(t)'])
axes[0].set_xlabel('')
ds_sel['mag'].plot(ax=axes[1])
axes[1].set_ylabel('U(t) [V]')
axes[1].set_xlabel('')

mag_0 = ds_sel['mag_0'].item()
axes[1].axhline(mag_0, color='k', linestyle='--')

ds_sel['T_abs'].plot(ax=axes[2])
axes[2].set_ylabel('T(t)')

ds_sel['AS_abs'].plot(ax=axes[3], label='AS')
ds_sel['dAS_abs'].plot(ax=axes[3], label='$\Delta$AS')
axes[3].set_ylabel('AS')
axes[3].legend()

for ax in axes:
    ax.set_title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'mwt_processing_overview.png'), bbox_inches='tight')



# %%

ds_sel

ds_sel['mag'].plot(marker='o')

plt.xlim(-1,1)