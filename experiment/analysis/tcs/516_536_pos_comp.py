#%%
from mhdlab.analysis.standard_import import *
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

from mhdlab.analysis.mwt.fitting import pipe_fit_exp

da_fit = ds_lecroy['dAS_abs']

# da_fit = da_fit.sel(run=[('2023-05-18',1),('2023-05-24',1)]).dropna('motor','all')

pipe_fit_exp.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(3, 'us'), Quantity(15, 'us'))

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_exp(
                                da_fit, 
                                )




#%%


da_plot = ds_mwt_fit[['AS_fit','AS_all']]

# because there is just one phi per run selected here, averaging over run is the same as selecting the phi
da_plot = da_plot.sel(run=[('2023-05-18',1),('2023-05-24',1)]).dropna('motor','all').mean('run')

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


SFR = ds_lecroy['SFR_abs'].unstack().stack(run=['date','run_num','phi']).dropna('run','all')

decay = ds_p['decay'].unstack().stack(run=['date','run_num','phi']).dropna('run','all')

ds_plot = xr.merge([SFR, decay])
ds_plot = ds_plot.assign_coords(run_plot = ('run', [f'{r[0]}_{r[1]} Phi={r[2]}' for r in ds_plot.indexes['run'].values]))

ds_plot
# ds_p

#%%

fig, axes= plt.subplots(2, figsize=(3,8), sharex=True)

for run in ds_plot.run.values:
    ds_sel = ds_plot.sel(run=run).dropna('motor','all')

    marker = 'x' if run[2] > 0.7 else 'o'
    linestyle = '--' if run[2] > 0.7 else '-'

    label = f'{run[0]}_{run[1]} Phi={run[2]}'
    ds_sel['SFR_abs'].plot(marker=marker, label=label, ax=axes[0], linestyle=linestyle)
    ds_sel['decay'].plot(marker=marker, label=label, ax=axes[1], linestyle=linestyle)

axes[0].set_ylabel('SFR')
axes[1].set_ylabel('Decay Time (us)')

axes[1].set_xlabel('Stage Position (mm)')
axes[0].set_title('')
axes[1].set_title('')
axes[0].legend()

axes[0].set_ylim(0, 200)

# SFR maximized for phi = 0.8
decay_08 = ds_p['decay'].sel(phi=0.8, method='nearest').mean('run').dropna('motor',how='all')
decay_08 = decay_08.sel(motor=180,method='nearest')

decay_06 = ds_p['decay'].sel(phi=0.6, method='nearest').mean('run').dropna('motor', how='all')
decay_06 = decay_06.sel(motor=155,method='nearest')

ratio = (decay_08/decay_06).item()
ratio

# Make a label on axes 2 with the ratio in the upper right corner

label = r'$\frac{{\tau_{0.79}(180mm)}}{{\tau_{0.59}(155mm)}}$ = ' + f'{ratio:.2f}'

axes[1].text(0.95, 0.95, label, transform=axes[1].transAxes, ha='right', va='top', fontsize=14)




plt.savefig(pjoin(DIR_FIG_OUT, '516_536_pos_decay_compare.png'))

#%%

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



# %%


# %%
