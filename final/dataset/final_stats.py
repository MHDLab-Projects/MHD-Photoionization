#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu #Sets matplotlib style
import matplotlib.ticker as mticker



# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()

ds_p_stats = ds_p_stats.pint.quantify()

# %%

#%%

m = ds_p_stats['mwt_fit_decay_exp_mean'].mean('kwt', keep_attrs=True).item()

s = ds_p_stats['mwt_fit_decay_exp_mean'].std('kwt', keep_attrs=True).item()

print(f"Mean: {m:.2f} +/- {s:.2f} us")

#%%

m = ds_p_stats['mwt_fit_decay_exp_mean'].sel(kwt=1, method='nearest').item()

s = ds_p_stats['mwt_fit_decay_exp_stderr'].sel(kwt=1, method='nearest').item()

print(f"Mean: {m:.2f} +/- {s:.2f} us")