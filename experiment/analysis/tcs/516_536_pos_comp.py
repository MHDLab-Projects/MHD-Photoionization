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

from mhdpy.analysis.mwt.fitting import pipe_fit_exp

da_fit = ds_lecroy['dAS_abs']

da_fit = da_fit.sel(run=[('2023-05-18',1),('2023-05-24',1)]).dropna('motor','all')

pipe_fit_exp.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(3, 'us'), Quantity(15, 'us'))

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_exp(
                                da_fit, 
                                )



#%%

da_plot = ds_mwt_fit[['AS_fit','AS_all']].mean('run')

da_plot = da_plot.rename({'AS_fit':'Fit', 'AS_all':'Data'})

da_plot = da_plot.to_array('Variable')

g= da_plot.plot(hue='Variable', row='motor', col='phi', figsize=(3,9))

for ax in g.axs[:,0]:
    ax.set_ylabel(r'$\Delta AS (V)$')

for ax in g.axs[:, 1]:
    # get the existing right axis label
    txt = ax.texts[0]
    # remove 'motor = ' from the text
    txt.set_text(txt.get_text().split('=')[1])


# get the legend
leg = g.figlegend
leg.set_bbox_to_anchor([0.5, 0.1])

plt.yscale('log')
plt.xlim(-1, 40)
plt.ylim(1e-4, 2)

plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, '516_536_pos_fits.png'))

#%%


#%%

# %%

fig, axes= plt.subplots(2, figsize=(3,8), sharex=True)

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
plt.ylim(0, 15)

plt.savefig(pjoin(DIR_FIG_OUT, '516_536_pos_decay_compare.png'))


# %%
