#%%


from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mwt
from mhdpy.analysis import absem
from mhdpy.plot import dropna

# plt.rcParams.update({'font.size': 16})

# %%

tc = '516_pos'

ds_absem = ppu.fileio.load_absem(tc)

ds_lecroy = ppu.fileio.load_lecroy('516_pos', AS_calc='relative')

ds_lecroy = ds_lecroy.isel(run=0) # There is only one 516

# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)

#%%

ds_absem

#%%

#%%

ds_absem.sel(mp='barrel')['alpha'].mean('mnum').mean('wavelength').dropna('run', how='all')

#%%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='motor', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)

#%%

ds_absem['alpha'].sel(mp='mw_horns').mean('mnum').mean('run').plot(hue='motor', x='wavelength')

plt.gca().get_legend().remove()

plt.ylim(-0.1,0.1)


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2


ds_fit = ds_absem.mean('mnum')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)
ds_alpha_fit['alpha'] = ds_fit['alpha']

#%%

ds = ds_alpha_fit[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%

g  = ds_p['nK_m3'].plot(hue='run_plot', x='motor',col='mp', marker='o')

dropna(g)

plt.yscale('log')


# %%

ds_p['nK_m3'].sel(mp='mw_horns').plot(hue='run_plot', x='motor', marker='o')

plt.yscale('log')

#%%[markdown]

# # 516 MWT Analysis

# only on 2023-05-24

#%%


da_max = ds_lecroy['AS_rel'].mwt._pulse_max()



da_max.mean('mnum').plot()


#%%

ds_lecroy['AS_rel'].mean('mnum').plot(row='motor')

plt.yscale('log')


# %%

from mhdpy.analysis.mwt.fitting import pipe_fit_mwt_3

da_fit = ds_lecroy['dAS_rel'].mean('mnum')


pipe_fit_mwt_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(0, 'us'), Quantity(25, 'us'))

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_3(da_fit)

ds_p['decay'] = 1/ds_p['krm']

# %%


ds_mwt_fit[['AS_fit','AS_all']].to_array('var').plot(hue='var', col='motor', col_wrap=3)

#%%

ds_p

#%%

ds_p['decay'].plot(marker='o')

plt.ylabel('Decay Time (us)')
plt.ylim(0, 10)

plt.savefig(pjoin(DIR_FIG_OUT, '516_mwt_decaytime.png'))

#%%

ds_mwt_fit[['AS_fit','AS_all']].to_array('var').plot(hue='var', col='motor', col_wrap=3)

plt.yscale('log')

plt.ylim(1e-4,)

plt.savefig(pjoin(DIR_FIG_OUT, '516_mwt_AS_fit.png'))
# %%
