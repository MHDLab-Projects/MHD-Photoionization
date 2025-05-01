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

run_count = ds.count('run')

#%%

da=ds['nK_m3_barrel'].std('run', keep_attrs=True)
da_stderr = da/run_count['nK_m3_barrel']
da_stderr.plot(label='Date Std. Err.', marker='o')

da_fit_stderr = ds['nK_m3_barrel_stderr'].mean('run')
da_fit_stderr.plot(label = 'Average Fit Std. Err.', marker='o')

plt.legend()
plt.title('AAS')

#%%

ratio = da_fit_stderr/da_stderr
ratio.plot(marker='o')

plt.ylabel('Ratio Fit SE/Run SE, AAS')

#%%

da = ds['mwt_fit_decay_exp'].std('run', keep_attrs=True)
da_stderr = da/run_count['mwt_fit_decay_exp']
da_stderr.plot(label='Date Std. Err.', marker='o')

da_fit_stderr = ds['mwt_fit_decay_exp_stderr'].mean('run')
da_fit_stderr.plot(label = 'Average Fit Std. Err.', marker='o')

plt.legend()
plt.title('MWT')

#%%


ratio = da_fit_stderr/da_stderr
ratio.plot(marker='o')

plt.ylabel('Ratio Fit SE/Run SE, MWT')


