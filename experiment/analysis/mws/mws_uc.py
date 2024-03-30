#%%[markdown]
# # Lecroy Uncertainty analysis

#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdpy.fileio.ct import load_df_cuttimes


#%%

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute', df_ct_downselect=df_cuttimes_seedtcs)

ds_lecroy = ds_lecroy.dropna('run', how='all')

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
        plt.text(0.5, 0.8, '{}'.format(int(count)), transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')



g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt')

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

plt.xlim(-10,10)


#%%

ds_l_unc = ds_lecroy.wma.initialize_stat_dataset('dAS_abs', 'mnum')

ds_l_unc = ds_l_unc.where(ds_l_unc['count'] > 0).dropna('run', how='all')

ds_l_unc['count'] = ds_l_unc['count'].isel(time=0)
ds_l_unc

#%%

g = xr.plot.FacetGrid(ds_l_unc, col='run', row='kwt', figsize=(10,10))

# Apply the custom function to each facet

g.map(custom_plot, 'time', 'mean', 'std', 'count')

for ax in g.axes[:,0]:
    ax.set_ylabel('$\Delta AS$')

for ax in g.axes[-1,:]:
    ax.set_xlabel('Time ($\mu$s)')

plt.xlim(-1,10)

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_53x_std_individual_count.png'))

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