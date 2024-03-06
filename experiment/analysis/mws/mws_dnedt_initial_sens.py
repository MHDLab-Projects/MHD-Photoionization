#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem
from mhdpy.coords.ct import downselect_acq_time

#%%


from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x')

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

ds_fit = ds_lecroy.mean('mnum')
da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

# %%

da_fit

#%%


from mhdpy.analysis.mws.fitting import gen_model_dnedt_v2

pre_norm_cutoff=5e-4
fit_timewindow=slice(Quantity(0, 'us'),Quantity(25, 'us'))
method='iterative'
take_log=False



params
#%%k

dnes = [1, 100, 10000, 1000000]

dss_p = []

for dne in dnes:
    print(dne)

    mod, params = gen_model_dnedt_v2(take_log=take_log)
    params['dne'].value = dne

    da_fit2 = da_fit.mws.fit_prep(pre_norm_cutoff=pre_norm_cutoff, remove_negative=False, min_mnum=None)

    ds_mws_fit, ds_p, ds_p_stderr = da_fit2.mws.perform_fit(mod, params, fit_timewindow=fit_timewindow, method=method)

    dss_p.append(ds_p.assign_coords(dne_init=[dne]))


# %%

dne

# %%

ds_p_all = xr.concat(dss_p, 'dne_init')

#%%

da_plot= ds_p_all[['krm','dne']].to_array('fit_var').dropna('run',how='all')

da_plot.plot(row='run', col='fit_var', hue='dne_init', marker='o', sharey=False)

dropna(g)
# %%
