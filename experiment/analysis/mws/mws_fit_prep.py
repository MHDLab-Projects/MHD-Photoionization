#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

# %%

tc = '53x'
ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

ds_lecroy = ds_lecroy.drop(0,'kwt')

da_sel = ds_lecroy['AS']

# %%

counts = da_sel.isel(time=0).count('mnum')
counts.plot(hue='run_plot', x='kwt', marker='o')

plt.yscale('log')
plt.xscale('log')

#%%

da_sel = da_sel.where(counts > 100)


#%%



da_sel.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')

plt.xlim(-1,30)

#%%

da_sel2 = da_sel.where(da_sel.mws._pulse_max() > 5e-4) # Targeting low power...
da_sel2 = da_sel2/da_sel.mws._pulse_max()

da_sel2

#%%
da_fit = da_sel2.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')

#%%

# Maximum that after pulse can be 
cutoff = 5e-1

aft_pulse = abs(da_sel2.sel(time=slice(40,50))).mean('time')
bef_pulse = abs(da_sel2.sel(time=slice(-50,-20))).mean('time')
during_pulse = abs(da_sel2.sel(time=slice(0,20))).mean('time')

# log spaced bins 

bins = np.logspace(np.log10(1e-4), np.log10(2), 300)

aft_pulse.plot.hist(bins=bins)
bef_pulse.plot.hist(bins=bins)
during_pulse.plot.hist(bins=bins, alpha=0.4)

plt.axvline(cutoff)

plt.yscale('log')
plt.xscale('log')

#%%

da_sel3 = da_sel2.where(aft_pulse < cutoff).where(bef_pulse < cutoff)


# %%

da_fit = da_sel3.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')
# %%

da_sel3.count('mnum').isel(time=0).plot(hue='run_plot', x='kwt', marker='o')

plt.xscale('log')
plt.yscale('log')
# %%[markdown]

# Compare this to mws.fit_prep method, should be the same. 

#%%

da_fit_prep = da_sel.mws.fit_prep()
# %%

da_fit_prep.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')