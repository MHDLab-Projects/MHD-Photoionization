

#%%[markdown]

# # Investigating the shape of the mwt data for individual mnum

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.fileio.ct import load_df_cuttimes

# %%


fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute', df_ct_downselect=df_cuttimes_seedtcs)
da_sel = ds_lecroy['dAS_abs'].sel(kwt=1, method='nearest')

da_sel = da_sel.groupby('run_plot').apply(lambda x: x.xr_utils.assign_mnum('mnum'))
da_sel = da_sel.dropna('run', how='all')
#%%


mnums = da_sel.coords['mnum'].values
n_mnums = 7

# pick an evenly space selection of mnums

mnums_sel = np.linspace(mnums.min()+1, mnums.max(), n_mnums).astype(int)

mnums_sel


#%%

da_sel.sel(mnum=mnums_sel).plot(hue='mnum', col='run', col_wrap=2)

plt.xlim(-1,10)
plt.ylim(1e-3,)
plt.yscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, 'mwt_kwt_mnum_sample.png'))


#%%

da_2 = da_sel/da_sel.mwt._pulse_max()

da_2.sel(mnum=mnums_sel).plot(hue='mnum', row='run')


plt.xlim(-1,10)
plt.ylim(1e-1,)
plt.yscale('log')

#%%

da_2.mean('mnum').plot(hue='run_plot', x='time')
plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')

#%%

da_3 = da_sel.where(da_sel.mwt._pulse_max() > 5e-4)

da_3 = da_3/da_3.mwt._pulse_max()

da_3.mean('mnum').plot(hue='run_plot', x='time')
plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')
# %%

da_sel

#%%

run_max = da_sel.sel(time=slice(-2,2)).groupby('run_plot').max('time').max('mnum')

run_max

#%%

da_hp = da_sel/run_max

da_hp.isel(mnum=mnums_sel).plot(row='run', hue='mnum', x='time')

plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')
#%%

da_hp = da_hp.where(da_hp.mwt._pulse_max() > 0.3)
# %%

da_hp.mean('mnum').plot(hue='run_plot', x='time')
# plt.xlim(-1,30)
plt.ylim(1e-3,)
plt.yscale('log')
# %%
ds_lecroy = ppu.fileio.load_lecroy('536_pos', AS_calc='absolute')

ds_lecroy

#%%
ds_sel = ds_lecroy.sel(run=('2023-05-18', 1))
da_sel = ds_sel['AS_abs']
da_sel = da_sel.dropna('motor', how='all')

motor_sel = [30, 50, 80, 105, 130, 150, 180, 210]

da_sel = da_sel.sel(motor = motor_sel, method='nearest')
da_sel = da_sel.dropna('mnum', how='all')

# da_sel = da_sel.xr_utils.assign_mnum('mnum')    
da_sel = da_sel.groupby('motor').apply(lambda x: x.xr_utils.assign_mnum('mnum'))

da_sel

#%%


mnums = da_sel.coords['mnum'].values
n_mnums = 5

# pick an evenly space selection of mnums

mnums_sel = np.linspace(mnums.min()+1, mnums.max(), n_mnums).astype(int)

mnums_sel


# da_sel['AS_abs']
# %%

da_sel.sel(mnum=mnums_sel).plot(hue='mnum', col='motor', col_wrap=2)


plt.savefig(pjoin(DIR_FIG_OUT, 'mwt_536_pos_mnum_sample.png'))
# plt.
# %%
