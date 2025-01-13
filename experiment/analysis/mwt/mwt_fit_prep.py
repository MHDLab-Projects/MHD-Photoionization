#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdlab.analysis import mwt

# %%

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='relative')

da_sel = ds_lecroy['dAS_rel']

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

da_sel2 = da_sel.where(da_sel.mwt._pulse_max() > 5e-4) # Targeting low power...
da_sel2 = da_sel2/da_sel.mwt._pulse_max()

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

# Compare this to mwt.fit_prep method, should be the same. 

#%%

da_fit_prep = da_sel.mwt.fit_prep()
# %%

da_fit_prep.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')