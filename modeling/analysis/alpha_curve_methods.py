# %% [markdown]
# # Alpha = 1 Curve Methods
# 
# Exploring methods of finding smooth curves where alpha = 1
# 

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
    'font.sans-serif': 'arial', 
    # 'figure.figsize': (4.6, 3)
})

from dotenv import load_dotenv
load_dotenv()
REPO_DIR = os.getenv('REPO_DIR')

from noneq_utils import chdir_if_nb_render; chdir_if_nb_render()
from noneq_utils import abscs, noneq

cantera_data_dir = os.path.join(REPO_DIR, 'modeling','dataset','output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf')).sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf')).sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf')).sel({'phi': 0.7, 'Kwt': 0.001})


G_th = ds_TP_params['Gth']
kr = ds_TP_params['kr']

coords_P_in = np.array([0, *np.logspace(-20,10,31)])

eta = 1

r = xyzpy.Runner(noneq.calc_G_NE, constants = {'eta' :eta }, var_names = ['G_NE'])
G_NE = r.run_combos({'P_in' : coords_P_in})
G_NE = G_NE['G_NE']

ne = noneq.calc_ne(G_th, G_NE, kr)

mue_cant = ds_TP_params['mobility']*10000
u=1e5
B=5e-4

alpha = noneq.calc_alpha(ne, mue_cant, kr, u, B, eta) #Why does eta need to be specified here as well...o

beta = alpha -1 
beta = beta.sel(P_in=alpha.coords['P_in'][::4])

# %%


cmap = plt.get_cmap('RdBu')


abs(beta).plot(col = 'P_in', col_wrap = 4, vmin=-1,vmax=1, xscale='log', cmap=cmap)

# %%

abs(beta).idxmin('P').plot(hue='P_in', col_wrap=4, y = 'T', xscale='log')

# %%

min_curve = abs(beta).idxmin('P')

fig, axes = plt.subplots(2, int(len(beta.coords['P_in'])/2), figsize=(12,4), sharey=True, sharex=True)
 
for i, ax in enumerate(axes.flatten()):
    min_curve.isel(P_in=i).plot(ax=ax, y='T',xscale='log', color='g')

    beta.isel(P_in=i).plot(ax=ax, vmin=-1,vmax=1,xscale='log', cmap=cmap)

    ax.set_title('')




# %%

beta_sel = beta.sel(P_in=0)
# beta_sel = abs(beta_sel)

beta_sel.plot(vmin=-1,vmax=1,xscale='log', cmap=cmap)


# %% [markdown]
# # Linear interpolation
# 
# For each temperature, we linearly interpolate the pressure at which alpha is zero. This gives a smooth curve.

# %%

T_sel = 2000

from scipy.interpolate import InterpolatedUnivariateSpline

xs = beta_sel.sel(T=T_sel,method='nearest').coords['P']
ys = beta_sel.sel(T=T_sel,method='nearest').values

spl = InterpolatedUnivariateSpline(xs, ys)
spl.roots()

plt.plot(xs,ys, marker='o')
plt.ylim(-1,1)
plt.xscale('log')
# plt.xlim(1e5,3e5)

plt.vlines(spl.roots(), -1,1)
# plt.hlines(0, 1e5,3e5)



# %%

Ps_zero = []
Ts = beta_sel.coords['T']

for T in Ts:
    beta_sel_T = beta_sel.sel(T=T)
    xs = beta_sel_T.coords['P'].values
    ys = beta_sel_T.values
    spl = InterpolatedUnivariateSpline(xs, ys)

    roots = spl.roots() 
    if len(roots) == 0:
        Ps_zero.append(np.nan)
    elif len(roots) == 1:
        Ps_zero.append(roots[0])
    else:
        raise ValueError

# Ps_zero


# %%

beta_sel.plot(vmin=-1,vmax=1,xscale='log', cmap=cmap)
plt.plot(Ps_zero, Ts)
