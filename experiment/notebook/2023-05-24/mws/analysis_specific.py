#%%
from mhdpy.analysis.standard_import import *

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%


tc = '536_power'

from mhdpy.mws_utils import calc_mag_phase_AS


fp_in = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy', '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)

ds = calc_mag_phase_AS(ds_in)

ds = ds.sel(date='2023-05-24').sel(run_num=1)


tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
# mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

ds = ds.mean('mnum', keep_attrs=True)

#%%


# %%

from mhdpy.mws_utils.fitting import fit_fn 
from mhdpy.analysis.xr import fit_da_lmfit
from lmfit import Model

da_fit = ds['AS'].copy()


pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

tc_dim = [dim for dim in da_fit.dims if dim != 'time'][0]

da_fit = da_fit.where(pulse_max > 5e-4) # Targeting low power...
da_fit = da_fit.dropna(tc_dim, how='all')
pulse_max = da_fit.sel(time=slice(-1,1)).max('time')

da_fit = da_fit/pulse_max

da_fit = np.log(da_fit).dropna('time', 'all')

mod = Model(lambda x, dne, ne0: fit_fn(x, dne, ne0, take_log=True))

params = mod.make_params()
params['dne'].value = 1e14
params['ne0'].value = 3e12
params['dne'].min = 0
params['ne0'].min = 0


da_fit_region = da_fit.sel(time=slice(0,30))
x_eval = da_fit.sel(time=slice(0,30)).coords['time'].values

# da_fit_coarse = da_fit_region.coarsen(time=10).mean()
fits, ds_p, ds_p_stderr = fit_da_lmfit(da_fit_region, mod, params, 'time', x_eval)

fits = np.exp(fits)
da_fit = np.exp(da_fit)

#%%

import re

plt.rcParams.update({'font.size': 14})


# plt.figure()

fig, axes = plt.subplots(3, sharex=True)

tc_dim = 'power'

for i, tc_val in enumerate(da_fit.coords[tc_dim]):

    ax = axes[i]

    ds_p_sel = ds_p.sel({tc_dim: tc_val})
    ds_p_stderr_sel = ds_p_stderr.sel({tc_dim: tc_val})

    param_str = "$ n_{e,0}" + ": {:.1e} \pm {:.1e} \#/cm^3 $  \n $ \Delta n_e: {:.1e} \pm {:.1e}  \#/cm^3$".format(
        ds_p_sel['ne0'].item(),
        ds_p_stderr_sel['ne0'].item(),
        ds_p_sel['dne'].item(), 
        ds_p_stderr_sel['dne'].item(), 
        )
    param_str = re.sub("e\+(\d\d)", "e^{\\1}", param_str)

    param_str


    # pre_pulse_std = da_fit.sel(time=slice(-50,-1)).std('time')
    # plt.hlines(pre_pulse_std.sel({tc_dim: tc_val}).item(), 0, 50, linestyle = '--')

    da_fit.sel({tc_dim: tc_val}).plot(marker='o', label='Data', ax=ax)
    fits.sel({tc_dim: tc_val}).plot(linestyle='--', label='Fit', color='red', ax=ax)

    ax.set_yscale('log')
    # ax.set_ylim(9e-3,2)
    # ax.set_xlim(-0.5,30)
    # ax.set_title('')

axes[0].legend()


plt.ylabel("Norm. Absorption + Scattering")

# plt.text(0, 0.8e-3, param_str)

# if i != len(tc_vals) -1: ax.set_xlabel('')

# if i != 0:
#     ax.set_title('')


plt.savefig('output/fit_specific.png')
# %%
