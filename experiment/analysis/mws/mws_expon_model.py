# %%[markdown]

# # exponential fit mws data

#%%

from mhdpy.analysis.standard_import import *
from pint import Quantity
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()

da_sel = ds_lecroy['AS'].sel(kwt=1, method='nearest')
#%%

from lmfit.models import ExponentialModel

mod = ExponentialModel()

params = mod.make_params()

params['amplitude'].value = 1
params['amplitude'].vary = True # dne peak obscures normalization 
params['decay'].value = 10

params['amplitude'].min = 0
params['decay'].min = 0

mod.take_log = False

#%%

from mhdpy.xr_utils import fit_da_lmfit

da_fit = da_sel.mean('mnum')

ds_mws_fit, ds_p, ds_p_stderr = da_fit.mws.perform_fit(
    mod, params,
    fit_timewindow=slice(Quantity(5, 'us'),Quantity(20, 'us')) 
    )

ds_mws_fit = xr.merge([ds_mws_fit, da_fit.rename('AS_all')])

#%%

ds_mws_fit.to_array('var').plot(row='run',hue='var')

plt.yscale('log')
# %%

ne0 = Quantity(2e12, 'cm**-3').to('um**-3').magnitude

# tau = 1/(2*kr*ne0)
# kr = 1/(2*tau*ne0)

ds_p['kr'] = 1/(2*ne0*ds_p['decay'])

ds_p['kr']
# %%

ds_p['decay']
