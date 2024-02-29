#%%[markdown]
# # Lecroy Uncertainty analysis

#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor


#%%
tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('mag', 'mnum')
ds_l_unc['count'] = ds_l_unc['count'].isel(time=0)
ds_l_unc

#%%

# make facetgrid plot with confidence intervals

# Define your custom plotting function
def custom_plot(time, mean, std, count, **kwargs):
    lower_bound = mean - std
    upper_bound = mean + std
    plt.plot(time, mean, **kwargs)
    plt.fill_between(time, lower_bound, upper_bound, alpha=0.5, color='gray')  # alpha sets the transparency
    # count should be a single value, add it as text to the plot in the upper left
    if not np.isnan(count):
        plt.text(0.05, 0.9, '{}'.format(int(count)), transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')



g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

plt.xlim(-10,10)

#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('AS', 'mnum')
ds_l_unc['count'] = ds_l_unc['count'].isel(time=0)
ds_l_unc

#%%

g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

plt.xlim(-10,10)

#%%

g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'stderr', 'count')


#%%

ds_l_unc_2 = ds_l_unc.wma.calc_weighted_mean('run')

#%%

g = xr.plot.FacetGrid(ds_l_unc_2, row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

# plt.ylim(0,1)