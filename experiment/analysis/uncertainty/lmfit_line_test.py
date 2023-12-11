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
