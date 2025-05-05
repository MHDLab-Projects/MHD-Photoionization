"""
Analysis comparing the fitting error for a given run the error between runs. 
"""

#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

#%%

fp = pjoin(REPO_DIR, r'final\dataset\output\53x_ds_p_full.cdf')

ds = xr.load_dataset(fp)
ds = ds.xr_utils.stack_run()
ds = ds.pint.quantify()

ds['nK_m3_barrel_stderr'] = ds['nK_m3_barrel_stderr'].pint.quantify('particle/m^3').pint.to('particle/cm^3')
ds['nK_m3_barrel_stderr'].attrs['long_name'] = 'SE $n_K [#/m^3]'
ds['nK_m3_mw_horns_stderr'] = ds['nK_m3_mw_horns_stderr'].pint.quantify('particle/m^3').pint.to('particle/cm^3')
ds['mwt_fit_decay_exp_stderr'] = ds['mwt_fit_decay_exp_stderr'].pint.quantify('us')

ds['nK_m3_barrel'] = ds['nK_m3_barrel'].pint.to('particle/cm^3')
ds['nK_m3_mw_horns'] = ds['nK_m3_mw_horns'].pint.to('particle/cm^3')

ds

#%%

from mhdlab.plot import dropna

g = ds['nK_m3_barrel'].plot(hue='run_plot', x ='kwt', marker='o')
plt.xscale('log')

# plt.yscale('log')
plt.ylim(-1e15,4e15)

dropna(g)

#%%

g = ds['nK_m3_mw_horns'].plot(hue='run_plot', x ='kwt', marker='o')
plt.yscale('log')
plt.xscale('log')

dropna(g)


#%%

run_count = ds.count('run')

#%%

ds

#%%

# Standard deviation
# Don't think calculating the std. dev. of the fit error is correct.

fig, axs = plt.subplots(2, 2, figsize=(8, 8))

# AAS Std. Err. Plot
da_std = ds['nK_m3_barrel'].std('run', keep_attrs=True)
da_std.plot(ax=axs[0, 0], label='$\sigma_{run}$', marker='o')

da_fit_stderr = ds['nK_m3_barrel_stderr'].mean('run', keep_attrs=True)
da_fit_std = da_fit_stderr* np.sqrt(run_count['nK_m3_barrel'])
da_fit_std.plot(ax=axs[0, 0], label='Average $\sigma_{fit}$', marker='o')

axs[0, 0].legend()
axs[0, 0].set_title('')
axs[0, 0].set_ylabel("$\sigma$ of $n_K\ from AAS fit [\\#/cm^3]$")


# AAS Ratio Plot
ratio = da_fit_std / da_std
ratio.plot(ax=axs[0, 1], marker='o')

axs[0, 1].set_ylabel('$\sigma_{fit}/\sigma_{run}$ for $n_K$ from AAS')
axs[0, 1].set_title('')
axs[0, 1].set_ylim(0, 1)

# MWT Std. Err. Plot
da_std = ds['mwt_fit_decay_exp'].std('run', keep_attrs=True)
da_std.plot(ax=axs[1, 0], label= '$\\sigma_{run}$', marker='o')

da_fit_stderr = ds['mwt_fit_decay_exp_stderr'].mean('run')
da_fit_std = da_fit_stderr* np.sqrt(run_count['mwt_fit_decay_exp'])
da_fit_std.plot(ax=axs[1, 0], label='Average $\\sigma_{fit}$', marker='o')

axs[1, 0].legend()
axs[1, 0].set_title('')

axs[1, 0].set_ylabel("$\sigma$ of $\\tau_{exp}\  [\mu s]$")

# MWT Ratio Plot
ratio = da_fit_std / da_std
ratio.plot(ax=axs[1, 1], marker='o')

axs[1, 1].set_ylabel('$\sigma_{fit}/\sigma_{run}$ for $\\tau_{exp}$')
axs[1, 1].set_title('')
axs[1, 1].set_ylim(0,1)

# Add abcd labels
labels = ['A)', 'B)', 'C)', 'D)']
for ax, label in zip(axs.flat, labels):
    ax.text(-0.1, 1.05, label, transform=ax.transAxes, va='bottom', ha='left', fontsize=12)
    ax.set_xlabel('$K_{wt, nominal} [\\%]$')

plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, 'error_prop_analysis_stddev.png'), dpi=300, bbox_inches='tight')



# %%

# Standard Error
# Comparing the standard error of the fit to the standard deviation of the run.

fig, axs = plt.subplots(2, 2, figsize=(8, 8))

# AAS Std. Err. Plot
da_std = ds['nK_m3_barrel'].std('run', keep_attrs=True)
da_stderr = da_std / np.sqrt(run_count['nK_m3_barrel'])
da_stderr.plot(ax=axs[0, 0], label='$SE_{run}$', marker='o')

da_fit_stderr = ds['nK_m3_barrel_stderr'].mean('run', keep_attrs=True)
da_fit_stderr.plot(ax=axs[0, 0], label='Average $SE_{fit}$', marker='o')

axs[0, 0].legend()
axs[0, 0].set_title('')
axs[0, 0].set_ylabel("SE of $n_K$ $[\\#/cm^3]$")


# AAS Ratio Plot
ratio = da_fit_stderr / da_stderr
ratio.plot(ax=axs[1, 0], marker='o')

axs[1, 0].set_ylabel('$SE_{fit}/SE_{run}$ for $n_K$')
axs[1, 0].set_title('')
axs[1, 0].set_ylim(0, 1)

# MWT Std. Err. Plot
da_std = ds['mwt_fit_decay_exp'].std('run', keep_attrs=True)
da_stderr = da_std / np.sqrt(run_count['mwt_fit_decay_exp'])
da_stderr.plot(ax=axs[0, 1], label= '$SE_{run}$', marker='o')

da_fit_stderr = ds['mwt_fit_decay_exp_stderr'].mean('run')
da_fit_stderr.plot(ax=axs[0, 1], label='Average $SE_{fit}$', marker='o')

axs[0, 1].legend()
axs[0, 1].set_title('')

axs[0, 1].set_ylabel("SE of $\\tau_{exp}\  [\mu s]$")

# MWT Ratio Plot
ratio = da_fit_std / da_std
ratio.plot(ax=axs[1, 1], marker='o')

axs[1, 1].set_ylabel('$SE_{fit}/SE_{run}$ for $\\tau_{exp}$')
axs[1, 1].set_title('')
axs[1, 1].set_ylim(0,1)

# Add abcd labels
labels = ['A)', 'B)', 'C)', 'D)']
for ax, label in zip(axs.flat, labels):
    ax.text(-0.1, 1.05, label, transform=ax.transAxes, va='bottom', ha='left', fontsize=12)
    ax.set_xlabel('$K_{wt, nominal} [\\%]$')

plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, 'error_prop_analysis_stderr.png'), dpi=300, bbox_inches='tight')
