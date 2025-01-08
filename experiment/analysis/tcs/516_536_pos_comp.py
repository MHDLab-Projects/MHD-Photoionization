#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

# plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_absem_536 = ppu.fileio.load_absem(tc).expand_dims('phi')
ds_lecroy_536 = ppu.fileio.load_lecroy(tc, avg_mnum=False, AS_calc='absolute').expand_dims('phi')

tc = '516_pos'

ds_absem_536 = ppu.fileio.load_absem(tc).expand_dims('phi')
ds_lecroy_516 = ppu.fileio.load_lecroy(tc, avg_mnum=False, AS_calc='absolute').expand_dims('phi')


ds_lecroy = xr.merge([ds_lecroy_536, ds_lecroy_516])

ds_lecroy = ds_lecroy.mwt.calc_time_stats()

ds_lecroy = ds_lecroy.mean('mnum', keep_attrs=True)

ds_lecroy

# %%

from mhdpy.analysis.mwt.fitting import pipe_fit_mwt_3

da_fit = ds_lecroy['dAS_abs']

pipe_fit_mwt_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(0, 'us'), Quantity(25, 'us'))

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_3(
                                da_fit, 
                                )

ds_p['decay'] = 1/ds_p['krm']

# %%


ds_mwt_fit[['AS_fit','AS_all']].mean('run').to_array('var').plot(hue='var', row='motor', col='phi')

plt.yscale('log')
plt.xlim(-1, 40)
plt.ylim(1e-3, 2)

#%%

ds_p

#%%

# %%

fig, axes= plt.subplots(2, figsize=(3,5), sharex=True)

g = ds_lecroy['SFR_abs'].mean('run').plot(hue='phi', marker='o', ax=axes[0])
dropna(g)

ds_p['decay'].mean('run').plot(marker='o', hue='phi', ax= axes[1])

axes[1].get_legend().remove()

for ax in axes:
    ax.set_title('')

axes[0].set_ylabel('SFR')
axes[1].set_ylabel('Decay Time (us)')
axes[1].set_xlabel('Stage Position (mm)')

plt.ylabel('Decay Time (us)')
plt.ylim(0, 20)

plt.savefig(pjoin(DIR_FIG_OUT, '516_536_decay_compare.png'))


# %%
