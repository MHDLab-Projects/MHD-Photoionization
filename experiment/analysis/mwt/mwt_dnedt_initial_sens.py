#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis import mwt
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem
from mhdpy.coords.ct import downselect_acq_time

#%%


from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='absolute')

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

ds_fit = ds_lecroy.mean('mnum')
da_fit = ds_fit['dAS_abs']

# %%

da_fit

#%%


from mhdpy.analysis.mwt.fitting import pipe_fit_mwt_2

dnes = [1, 100, 10000, 1000000]

dss_p = []

for dne in dnes:
    print(dne)

    pipe_fit_mwt_2.params['dne'].value = dne

    ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_2(da_fit)

    dss_p.append(ds_p.assign_coords(dne_init=[dne]))


# %%

ds_p_all = xr.concat(dss_p, 'dne_init')
ds_p_all['krm'] = ds_p_all['krb']*ds_p_all['ne0']

#%%

da_plot= ds_p_all[['krm','dne']].to_array('fit_var').dropna('run',how='all')

g = da_plot.plot(row='run', col='fit_var', hue='dne_init', marker='o', sharey=False)

dropna(g)
# %%
