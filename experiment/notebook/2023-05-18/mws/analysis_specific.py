#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-05-18'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)
# #%%


tc = '536_power'

from mhdpy.analysis import mws

fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)
ds_in = ds_in.sel(date=datestr).sel(run_num=1)

ds = ds_in.mws.calc_AS_rel().drop('mag_pp')


tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

ds = ds.mean('mnum', keep_attrs=True)

ds

# %%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_1 

da_fit = ds['AS'].copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_1(da_fit)

ds_mws_fit

#%%

plt.rcParams.update({'font.size': 14})


# plt.figure()

tc_vals = ds_mws_fit.coords[tc_dim]

fig, axes = plt.subplots(len(tc_vals), sharex=True, figsize=(5,10))

tc_dim = 'power'

for i, tc_val in enumerate(tc_vals):

    ax = axes[i]

    ds_p_sel = ds_p.sel({tc_dim: tc_val})
    ds_p_stderr_sel = ds_p_stderr.sel({tc_dim: tc_val})

    # pre_pulse_std = da_fit.sel(time=slice(-50,-1)).std('time')
    # plt.hlines(pre_pulse_std.sel({tc_dim: tc_val}).item(), 0, 50, linestyle = '--')

    ds_mws_fit['AS'].sel({tc_dim: tc_val}).plot(marker='o', label='Data', ax=ax)
    ds_mws_fit['AS_fit'].sel({tc_dim: tc_val}).plot(linestyle='--', label='Fit', color='red', ax=ax)

    ax.set_yscale('log')
    ax.set_ylim(9e-3,2)
    ax.set_xlim(-0.5,30)
    ax.set_title('')
    ax.set_ylabel('AS')


    param_str = 'Rel. Power = {}'.format(tc_val.item())
    ax.text(17, 0.2, param_str)

axes[0].legend()


# plt.ylabel("Norm. Absorption + Scattering")


# if i != len(tc_vals) -1: ax.set_xlabel('')

# if i != 0:
#     ax.set_title('')


plt.savefig(pjoin(DIR_FIG_OUT,'fit_specific.png'))
# %%


ds_mws_fit['AS'].plot(hue='power')

ax=plt.gca()

for i, tc_val in enumerate(ds_mws_fit.coords[tc_dim]):

    ds_mws_fit['AS_fit'].sel({tc_dim: tc_val}).plot(linestyle='--', label='Fit', color='red', ax=ax)

plt.yscale('log')

ax.set_ylim(9e-3,2)
ax.set_xlim(-0.5,30)
# %%
ds_p['dne'].plot(marker='o')
# plt.yscale('log')
# plt.xscale('log')

plt.ylabel("$\Delta n_e$")
plt.xlabel("Relative Power")
plt.title('')
# %%
