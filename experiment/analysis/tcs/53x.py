#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem
from mhdpy.coords.ct import downselect_acq_time

from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

#%%[markdown]

# # Absem

# %%

ds_absem = ppu.fileio.load_absem('53x')

#%%



ds_absem = downselect_acq_time(ds_absem, df_cuttimes_seedtcs)

# %%

ds_absem['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='kwt', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)

#%%

# # Examine calibration offset. This does not appear to be off as much as in the 5xx notebook. TODO: investigate and quantify off-peak calibraiton offset. 
# ds_absem['alpha'].mean('mnum').sel(mp='mw_horns').isel(run=-1).plot(hue='kwt')


#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit = ds_absem.mean('mnum')

ds_absem_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)


# %%
da_plot = ds_absem_fit[['alpha_red','alpha_fit']].to_array('var')

g = da_plot.sel(mp='barrel').plot(col='kwt',hue='var', row='run', ylim = (1e-3,2), yscale='log', figsize=(10,10))

for ax in g.axes.flatten():
    ax.hlines(1,760,775, linestyles='--', color='gray')

#%%

da = ds_p['nK_m3']
da.coords['kwt'].attrs = dict(long_name="Kwt", units='%')

# da = da.dropna('run', how='all')

g  = da.plot(hue='run_plot', x='kwt',col='mp', marker='o')

dropna(g)

plt.yscale('log')
plt.xscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, '53x_AES_nK.png'), dpi=300, bbox_inches='tight')

#%%

ds_nK = xr.merge([
    ds_p['nK_m3'].to_dataset(name='mean'),
    ds_p_stderr['nK_m3'].to_dataset(name='stderr'),
    ds_absem.mean('wavelength').count('mnum')['alpha'].to_dataset(name='count')
])

ds_nK['std'] = ds_nK['stderr']*np.sqrt(ds_nK['count'])

# ds_nK2 = ds_nK.wma.calc_weighted_mean('run')
# ds_nK2 = ds_nK.wma.initialize_stat_dataset('mean', 'run')
# ds_nK2

#%%[markdown]

# # Lecroy

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x', df_ct_downselect=df_cuttimes_seedtcs, AS_calc='relative', avg_mnum=True)

# ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

# ds_fit = ds_lecroy.mean('mnum')
# da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

da_fit = ds_lecroy['AS']


#%%

g = da_fit.plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')

dropna(g)

#%%[markdown]

# # Exponential Fit

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

#%%

ds_mws_fit[['AS_all','AS_fit']].to_array('var').plot(hue='var', row='kwt', col='run')

plt.yscale('log')

