"""Calculate curves for gamma =1 with linear interpolation of each pressure vs temperature curve"""
# %%
import os
import numpy as np
import xarray as xr
import xyzpy

from dotenv import load_dotenv
load_dotenv()

canterapath = os.path.join('output')
if not os.path.exists('output'): os.mkdir('output')
# %%
ds_TP_species = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(canterapath, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(canterapath, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_NE = xr.open_dataset('output/ds_NE.cdf').squeeze()
gamma = ds_NE['gamma']

# Add enhancement factor
da_dsigma_tot = xr.load_dataset('output/da_dsigma_tot.cdf')['enhancement factor']
gamma = gamma*da_dsigma_tot

# gamma.sel(Kwt=0.1, phi=1, analysis='PI_Bhall', P_in=0, T=1.6e3).plot()

beta = gamma -1 
# beta = beta.sel(P_in=gamma.coords['P_in'][::4])

# # beta = beta.sel(Kwt=0.01, phi=0.8, analysis='perf_Bconst')
# beta = beta.sel(Kwt=0.01, phi=0.8)


# # beta.coords['analysis'] = beta.coords['analysis'].astype(str)
# beta = beta.assign_coords(analysis=[0,1,2,3])

beta

# %%

from scipy.interpolate import InterpolatedUnivariateSpline


#TODO: figure out how to deal with multiple zeros for a given pressure or temperature.  

def calc_beta_zero(beta, **coords):

    beta_sel_T = beta.sel(**coords)
    xs = beta_sel_T.coords['P'].values
    ys = beta_sel_T.values
    spl = InterpolatedUnivariateSpline(xs, ys)

    roots = spl.roots() 
    if len(roots) == 0:
        return np.nan
    elif len(roots) == 1:
        return roots[0]
    else:
        # print("Got multiple roots: {}".format(roots))
        # print(beta_sel_T)
        # raise ValueError
        return np.nan

# calc_beta_zero(beta, *combo.values())


# %%

import xyzpy

r = xyzpy.Runner(calc_beta_zero,
    constants={'beta': beta},
    var_names = ['P_zero'])

combos = {
    'T' : beta.coords['T'].values,
    'phi': beta.coords['phi'].values,
    'Kwt': beta.coords['Kwt'].values,
    # 'P_in' : [0],
    # 'analysis': beta.coords['analysis'].values,
    'eta': beta.coords['eta'].values,
    'l_flow': beta.coords['l_flow'].values,
    'rxn': beta.coords['rxn'].values
}

ds_out = r.run_combos(combos)

ds_out.attrs = {}
ds_out.to_netcdf('output/P_zero.cdf')

# da_out.attrs
# %%


# da_out = ds_out['P_zero']
# da_out.plot(y='T', hue='P_in', col='phi',row='Kwt')
# plt.xscale('log')

# %%
