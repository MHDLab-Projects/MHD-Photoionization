#%%

import numpy as np
from lmfit import Model
import matplotlib.pyplot as plt

# Define constant function
def constant(x, c):
    return np.full_like(x, c)

# Generate x values
x = np.linspace(0, 10, 100)

# Generate y values with noise
c_true = 2
noise = np.random.normal(scale=0.5, size=x.size)
y = constant(x, c_true) + noise

# Create a model
model = Model(constant)

# Fit the model to the data
result = model.fit(y, x=x, c=1)

# Print fit parameters
print(result.params)

# Plot the data and the fit
plt.scatter(x, y, label='data')
plt.plot(x, result.best_fit, label='fit')
plt.legend()
plt.show()
# %%
# Define a list of noise scales
noise_scales = np.linspace(0.1, 2, 20)

# Initialize a list to store reduced chi squared values
chi_sq = []

# Loop over noise scales
for scale in noise_scales:
    # Generate y values with noise
    y = constant(x, c_true) + np.random.normal(scale=scale, size=x.size)
    
    # Fit the model to the data
    result = model.fit(y, x=x, c=1)
    
    # Store the reduced chi squared value
    chi_sq.append(result.redchi)

# Plot noise scale vs reduced chi squared
plt.plot(noise_scales, chi_sq, marker='o')
plt.xlabel('Noise scale')
plt.ylabel('Reduced chi squared')
plt.show()


#%%

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Define a list of numbers of data points
num_points = np.arange(50, 1001, 50)

# Initialize a 2D array to store reduced chi squared values
chi_sq = np.empty((len(noise_scales), len(num_points)))

# Loop over noise scales and numbers of data points
for i, scale in enumerate(noise_scales):
    for j, n in enumerate(num_points):
        # Generate x and y values with noise
        x = np.linspace(0, 10, n)
        y = constant(x, c_true) + np.random.normal(scale=scale, size=x.size)
        
        # Fit the model to the data
        result = model.fit(y, x=x, c=1)
        
        # Store the reduced chi squared value
        chi_sq[i, j] = result.redchi

# Plot noise scale vs number of data points vs reduced chi squared
plt.imshow(chi_sq, aspect='auto', origin='lower', norm=LogNorm(),
           extent=[num_points.min(), num_points.max(), noise_scales.min(), noise_scales.max()])
plt.colorbar(label='Reduced chi squared')
plt.xlabel('Number of data points')
plt.ylabel('Noise scale')
plt.show()


#%%

# Define a list of numbers of data points
num_points = np.arange(50, 1001, 50)

# Initialize a 2D array to store reduced chi squared values
std_err = np.empty((len(noise_scales), len(num_points)))

# Loop over noise scales and numbers of data points
for i, scale in enumerate(noise_scales):
    for j, n in enumerate(num_points):
        # Generate x and y values with noise
        x = np.linspace(0, 10, n)
        y = constant(x, c_true) + np.random.normal(scale=scale, size=x.size)
        
        # Fit the model to the data
        result = model.fit(y, x=x, c=1)
        
        # Store the reduced chi squared value
        # chi_sq[i, j] = result.redchi
        std_err[i, j] = result.params['c'].stderr

# Plot noise scale vs number of data points vs reduced chi squared
plt.imshow(std_err, aspect='auto', origin='lower', norm=LogNorm(),
           extent=[num_points.min(), num_points.max(), noise_scales.min(), noise_scales.max()])
plt.colorbar(label='Constant Standard Error')
plt.xlabel('Number of data points')
plt.ylabel('Noise scale')
plt.show()

# %%

# select 5 noise scales and plot them as a function of number of data points

std_err_select = std_err[[0,4,8,12,16],:]

plt.plot(num_points, std_err_select.T)

plt.legend(noise_scales[[0,4,8,12,16]], title='Noise Scale')

plt.xlabel('Number of data points')

plt.ylabel('Constant Standard Error')

#%%[markdown]

# # Simultaneous fit of multiple datasets

# https://lmfit.github.io/lmfit-py/examples/example_fit_multi_datasets.html

#%%


import matplotlib.pyplot as plt
import numpy as np

from lmfit import Parameters, minimize, report_fit


def gauss(x, amp, cen, sigma):
    """Gaussian lineshape."""
    return amp * np.exp(-(x-cen)**2 / (2.*sigma**2))


def gauss_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    amp = params[f'amp_{i+1}']
    cen = params[f'cen_{i+1}']
    sig = params[f'sig_{i+1}']
    return gauss(x, amp, cen, sig)


def objective(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - gauss_dataset(params, i, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()


np.random.seed(2021)
x = np.linspace(-1, 2, 151)
data = []
for _ in np.arange(5):
    amp = 0.60 + 9.50*np.random.rand()
    cen = -0.20 + 1.20*np.random.rand()
    sig = 0.25 + 0.03*np.random.rand()
    dat = gauss(x, amp, cen, sig) + np.random.normal(size=x.size, scale=0.1)
    data.append(dat)
data = np.array(data)

fit_params = Parameters()
for iy, y in enumerate(data):
    fit_params.add(f'amp_{iy+1}', value=0.5, min=0.0, max=200)
    fit_params.add(f'cen_{iy+1}', value=0.4, min=-2.0, max=2.0)
    fit_params.add(f'sig_{iy+1}', value=0.3, min=0.01, max=3.0)

for iy in (2, 3, 4, 5):
    fit_params[f'sig_{iy}'].expr = 'sig_1'


out = minimize(objective, fit_params, args=(x, data))
report_fit(out.params)


plt.figure()
for i in range(5):
    y_fit = gauss_dataset(out.params, i, x)
    plt.plot(x, data[i, :], 'o', x, y_fit, '-')

#%%

out.covar.shape

np.sqrt(np.diag(out.covar))

# %%


def constant_dataset(params, i):
    """Calculate constant value from parameters for data set."""
    c = params['c']
    return c

def objective(params, data):
    """Calculate total residual for fits of constant to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - constant_dataset(params, i)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

# ... (your data generation code here) ...

x = np.linspace(-1, 2, 151)
const_val = 0.5

data = []

for _ in np.arange(5):
    dat = np.full_like(x, const_val) + np.random.normal(size=x.size, scale=0.1)
    data.append(dat)

data = np.array(data)

data.shape

fit_params = Parameters()
fit_params.add('c', value=0.5, min=0.0, max=200)

out = minimize(objective, fit_params, args=(data,))
report_fit(out.params)

plt.figure()
for i in range(5):
    y_fit = constant_dataset(out.params, i)
    plt.plot(data[i, :], 'o', np.full_like(data[i, :], y_fit), '-')
# %%


x = np.linspace(-1, 2, 100)
const_val = 0.5
noise_scale = 0.1

const_vals = [np.random.normal(loc=1, scale=1) for _ in np.arange(5)]

data = []

for const_val in const_vals:
    dat = np.full_like(x, const_val) #+ np.random.normal(size=x.size, scale=0.1)
    data.append(dat)

data = np.array(data)

fit_params = Parameters()
fit_params.add('c', value=0.5, min=0.0, max=200)

out = minimize(objective, fit_params, args=(data,))
report_fit(out)

plt.figure()
for i in range(5):
    y_fit = constant_dataset(out.params, i)
    plt.plot(data[i, :], 'o', np.full_like(data[i, :], y_fit), '-')
# %%

out.params['c'].stderr*np.sqrt(out.ndata)

#%%

np.sqrt(out.covar)

#%%
# %%
np.mean(data)


#%%


def line_dataset(params, x):
    """Calculate line from parameters for data set."""
    m = params['m']
    c = params['c']
    return m * x + c

def objective(params, x, data):
    """Calculate total residual for fits of line to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - line_dataset(params, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

# ... (your data generation code here) ...

x = np.linspace(-1, 2, 100)

m_nominal = 1.0
m_vals = [np.random.normal(m_nominal, 1) for _ in np.arange(5)]
c_val = 0.5

data = []

for m_val in m_vals:
    dat = m_val*x + c_val#+ np.random.normal(size=x.size, scale=0.1)
    data.append(dat)

data = np.array(data)

fit_params = Parameters()
fit_params.add('c', value=c_val, min=0.0, max=200)
fit_params.add('m', value=m_nominal, min=0.0, max=200)

out = minimize(objective, fit_params, args=(x, data), scale_covar=True)
report_fit(out)

plt.figure()
for i in range(5):
    y_fit = line_dataset(out.params, x)
    plt.plot(x, data[i, :], 'o', x, y_fit, '-')

#%%

print("average m value: ", np.mean(m_vals))
print("std m value: ", np.std(m_vals))

#%%
