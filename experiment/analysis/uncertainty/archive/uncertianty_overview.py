#%%
from mhdlab.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdlab.analysis import mwt
# %%[markdown]

# Summary: 

# For each test case (e.g. 1% seed, microwave horns at goldilocks) different things are being continously acquired. 1D time data (heat transfer) or 2D time data (oscilliscope, absorption spectrometer). 
# These data are chopped up into groups with the test case timewindows, and for each time window there are a given number of repeated measurements. 

# we are going to examine absorption spectrometer data for the seed ramps (53x) as an example

# %%


tc = '53x'

ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))

ds_absem = ds_absem.xr_utils.stack_run()

ds_absem


# %%

# examine one specific test case and 1% seed

ds_absem2 = ds_absem.sel(mp='barrel').sel(kwt=1, method='nearest')

ds_absem2

#%%[markdown]

# The data is transformed from being indexed by time to being indexed by the test case coordinate (e.g. K mass frac) as well as the 'measurement number' 

#%%

ds_sel = ds_absem2.sel(run=('2023-05-18',1)).dropna('mnum', how='all')

ds_sel['led_on'].plot.line(x='wavelength', hue='mnum')

# %%[markdown]

# We then calculate the absorption (alpha). 
# In this funciton the mean and standard deviation is calculated over the measurement number for led on and off, then the difference is calculated with error propogration as defined in the scipp package. 
# https://scipp.github.io/reference/error-propagation.html
# This is then divided by the calibration data (no error assumed) to get alpha.

#%%

from mhdlab.analysis.absem import calc_alpha_scipp

ds_absem_in = ds_absem2.unstack('run')

#TODO: scipp error
ds_absem_in.coords['date'] = ds_absem_in.coords['date'].astype(str)

ds_alpha = calc_alpha_scipp(ds_absem_in)

ds_alpha = ds_alpha.xr_utils.stack_run()

ds_alpha

#%%[markdown]

# Plot with the standard deviation
# %%

# Define your custom plotting function
def custom_plot(wavelength, mean, std, count, **kwargs):
    lower_bound = mean - std
    upper_bound = mean + std
    plt.plot(wavelength, mean, **kwargs)
    plt.fill_between(wavelength, lower_bound, upper_bound, alpha=0.5, color='gray')  # alpha sets the transparency
    # count should be a single value, add it as text to the plot in the upper left
    if not np.isnan(count):
        plt.text(0.05, 0.9, '{}'.format(int(count)), transform=plt.gca().transAxes, color='red', fontsize=14, fontweight='bold')



# Create a FacetGrid
g = xr.plot.FacetGrid(ds_alpha, col='run')

# Apply the custom function to each facet
g.map(custom_plot, 'wavelength', 'mean', 'std', 'count')

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

#%%[markdown]

# Plot with the standard errror

#%%


# Create a FacetGrid
g = xr.plot.FacetGrid(ds_alpha, col='run')

# Apply the custom function to each facet
g.map(custom_plot, 'wavelength', 'mean', 'stderr', 'count')

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

#%%[markdown]

# We now want to calculate the weighted mean over the run coordinate.
# This uses a custom xarray acessor that I wrote.


#%%
from mhdlab.xr_utils import WeightedMeanAccessor

ds_alpha_2 = ds_alpha.wma.calc_weighted_mean('run')

#%%

custom_plot(ds_alpha_2['wavelength'], ds_alpha_2['mean'], ds_alpha_2['std'], ds_alpha_2['count'])

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

plt.title('Standard deviation')

#%%

custom_plot(ds_alpha_2['wavelength'], ds_alpha_2['mean'], ds_alpha_2['stderr'], ds_alpha_2['count'])

plt.ylim(-0.1,1.1)
plt.xlim(765,772)

plt.title('Standard error')