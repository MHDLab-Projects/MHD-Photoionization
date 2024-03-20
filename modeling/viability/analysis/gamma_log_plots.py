
# %%
import os
import numpy as np
import xarray as xr
import xyzpy

import matplotlib.pyplot as plt

plt.rcParams.update({
    "savefig.facecolor": 'white',
    "font.size": 11, 
    'savefig.dpi': 300, 
    # 'font.sans-serif': 'arial', 
    # 'figure.figsize': (4.6, 3)
})

from dotenv import load_dotenv
load_dotenv()
REPO_DIR = os.getenv('REPO_DIR')

from pi_paper_utils import abscs, noneq

cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')

import matplotlib.colors as colors
import numpy as np

#%%

ds_P_zero = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'P_zero.cdf'))

ds_NE = xr.open_dataset(os.path.join(PI_modeling_dataset_dir, 'ds_NE.cdf')).squeeze()
gamma = ds_NE['gamma']

# Add enhancement factor
da_dsigma_tot = xr.load_dataset(os.path.join(PI_modeling_dataset_dir,'da_dsigma_tot.cdf'))['enhancement factor']
gamma = gamma*da_dsigma_tot

# %%

# Old plot, no scaling

plt.figure(figsize=(5,3))

combo_sel = dict(l_b=0, P_in=0, Kwt=0.01, phi=0.8, eta='perf', rxn='mm_sum')

cmap = plt.get_cmap('RdBu')
gamma_sel = gamma.sel(combo_sel)
g = gamma_sel.plot(vmin=-1.2,vmax=1.2,xscale='log', cmap=cmap)


#%%

norm = colors.LogNorm(vmin=0.001, vmax=1000)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm)

ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')




#%%


norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=-1, vmax=1000, base=10)

# Set the colorbar label

g = (gamma_sel - 1).plot(xscale='log', cmap=cmap, norm=norm)

#%%

from matplotlib.ticker import NullLocator

class CustomNorm(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Logarithmically interpolate values from vmin to midpoint
        log_vmin_mid = np.ma.log10(np.ma.masked_greater(value, self.midpoint) - self.vmin + 1)
        log_mid = np.ma.log10(self.midpoint - self.vmin + 1)
        result = np.where(value <= self.midpoint, 0.5 * log_vmin_mid / log_mid, 0)
        # Logarithmically interpolate values from midpoint to vmax
        log_mid_vmax = np.ma.log10(np.ma.masked_less(value, self.midpoint) - self.midpoint + 1)
        log_vmax = np.ma.log10(self.vmax - self.midpoint + 1)
        result += np.where(value > self.midpoint, 0.5 + 0.5 * log_mid_vmax / log_vmax, 0)
        return np.ma.masked_array(result)

    def inverse(self, value):
        log_mid = np.ma.log10(self.midpoint - self.vmin + 1)
        log_vmax = np.ma.log10(self.vmax - self.midpoint + 1)
        result = np.where(value <= 0.5, (np.power(10, value * 2 * log_mid) + self.vmin - 1), 0)
        result += np.where(value > 0.5, (np.power(10, (value - 0.5) * 2 * log_vmax) + self.midpoint - 1), 0)
        return np.ma.masked_array(result)

norm = CustomNorm(vmin=0.01, midpoint=1, vmax=1000)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm)


g.colorbar.ax.yaxis.set_minor_locator(NullLocator())
# Set the colorbar ticks
# g.colorbar.set_ticks([norm(0.1), norm(1), norm(10), norm(100), norm(1000)])
g.colorbar.set_ticks([0.1, 1, 10, 100, 1000])

# Set the colorbar tick labels
g.colorbar.set_ticklabels(['0.1', '1', '10', '100', '1000'])

#%%

class CustomNorm(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Logarithmically interpolate values from vmin to midpoint
        log_v_lower = np.ma.log(np.ma.masked_greater(value, self.midpoint) )
        # log_mid = np.ma.log(self.midpoint - self.vmin)
        log_vmin = np.ma.log( self.vmin)
        result = np.where(value <= self.midpoint, 0.5 * (log_v_lower - log_vmin), 0)
        # Logarithmically interpolate values from midpoint to vmax
        log_mid_vmax = np.ma.log(np.ma.masked_less(value, self.midpoint))
        log_vmax = np.ma.log(self.vmax )
        result += np.where(value > self.midpoint, 0.5 + 0.5 * (log_mid_vmax), 0)
        return result

    # def inverse(self, value):
    #     result = np.where(value <= 0.5, (np.exp(value * 2 * log_mid) + self.vmin), 0)
    #     result += np.where(value > 0.5, (np.exp((value - 0.5) * 2 * log_vmax) + self.midpoint), 0)
    #     return result

norm = CustomNorm(vmin=0.1, midpoint=1, vmax=100)

g = (gamma_sel).plot(xscale='log', cmap=cmap, norm=norm)

#%%

# Attempt at custom norm with log and linear interpolation

class CustomNorm(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Linearly interpolate negative values
        result = np.where(value < self.midpoint, 0.5 * (value - self.vmin) / (self.midpoint - self.vmin), 0)
        # Logarithmically interpolate positive values
        log_v = np.ma.log(np.ma.masked_less_equal(value, 0) - self.midpoint + 1e-10) + 1
        log_vmax = np.ma.log(self.vmax - self.midpoint + 1e-10)  +1 
        result += np.where(value >= self.midpoint, 0.5 + 0.5 * log_v / log_vmax, 0)
        return result



# Create a colormap with white at the center
cmap = colors.ListedColormap(plt.cm.RdBu(np.linspace(0, 1, 256)))

# Create a custom normalizer
norm = CustomNorm(vmin=0, midpoint=1, vmax=100)

# Plot the data with the custom colormap and normalizer
g = (gamma_sel).plot(vmin=0, vmax=100, xscale='log', cmap=cmap, norm=norm)

# Set the colorbar label
g.colorbar.set_label('$\\gamma - 1$')

# Plot the green line
ds_P_zero['P_zero'].sel(combo_sel).plot(y='T', color='green', linewidth=1, linestyle='--')

# Set the title and labels
plt.gca().set_title('')
plt.ylabel('Temperature (K)')
plt.xlabel('Pressure (Pa)')

# Set the limits
plt.xlim(0.8e4,1.2e6)
plt.ylim(1200,3500)

# Save the figure
plt.savefig('output/gamma_curve_demo.png')


#%%

# %%

beta_sel.sel(T=2000, method='nearest').plot(xscale='log')

plt.yscale('log')


