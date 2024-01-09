# %%
import numpy as np

import scipp as sc

rng = np.random.default_rng(12345)
a = sc.array(
    dims=['x', 'y'], values=rng.random((2, 4)), variances=rng.random((2, 4)), unit='m'
)
b = sc.array(
    dims=['y', 'x'], values=rng.random((4, 2)), variances=rng.random((4, 2)), unit='s'
)
a / b

# %% [markdown]
# ### Propagation of uncertainties
# 
# If variables have variances, operations correctly propagate uncertainties (the variances), in contrast to a naive implementation using NumPy:

# %%
result = a/b
result.values

#%%

a.sum().variances

# a.mean('y').mean('x')

#%%

a2 = a.copy()

a2.variances = None

a2.sum()

# a2.mean('y').mean('x')

#%%

sc.compat.to_xarray(result)

# result.variances


# %%
a.values/np.transpose(b.values)

# %%
result.variances

# %%
a.variances/np.transpose(b.variances)

# %% [markdown]
# The implementation assumes uncorrelated data and is otherwise based on, e.g., [Wikipedia: Propagation of uncertainty](https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae).
# See also [Propagation of uncertainties](../reference/error-propagation.rst) for the concrete equations used for error propagation.

# %% [markdown]
# <div class="alert alert-info">
#     <b>Note:</b>
# 
# If an operand with variances is also [broadcast](#Broadcasting) in an operation then the resulting values would be correlated.
# Scipp has no way of tracking or handling such correlations.
# Subsequent operations that combine values of the result, such as computing the mean, would therefore result in *underestimated uncertainties*.
# 
# To avoid such silent underestimation of uncertainties, operations raise [VariancesError](https://scipp.github.io/generated/classes/scipp.VariancesError.html) when an operand with variances is implicitly or explicitly broadcast in an operations.
# See our publication [Systematic underestimation of uncertainties by widespread neutron-scattering data-reduction software](https://doi.org/10.3233/JNR-220049) for more background.
# </div>


#%%
# result.save_hdf5('test.hdf5')

#%%

# import h5py
# f = h5py.File('test.hdf5', 'r')

# f.keys()
#%%

# import xarray as xr
# 
# xr.open_dataset('test.hdf5')
