#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis.mws.fitting import calc_dnedt
from scipy.integrate import solve_ivp

tw_goldi = slice(Timestamp('2023-05-18 22:28:40.081841408'), Timestamp('2023-05-18 22:29:19.447003648'), None)

from mhdpy.fileio.path import gen_path
sp_dir = gen_path('sharepoint')

# %%

results_dir = pjoin(sp_dir, 'Data Share', 'Photoionization', 'Lecroy Resample Check', 'munged_nors', '2023-05-18', 'Lecroy')
fp_raw = pjoin(results_dir, 'ds_pos_1.cdf')

ds_raw = xr.load_dataset(fp_raw)

ds_raw = ds_raw.sel(acq_time=tw_goldi)

ds_raw.pint.dequantify().to_netcdf('output/ds_raw.cdf')

#%%


results_dir = pjoin(sp_dir, 'Data Share', 'Photoionization', 'Lecroy Resample Check', 'munged_10', '2023-05-18', 'Lecroy')
fp_raw = pjoin(results_dir, 'ds_pos_1.cdf')

ds_raw = xr.load_dataset(fp_raw)

ds_raw = ds_raw.sel(acq_time=tw_goldi)

ds_raw.pint.dequantify().to_netcdf('output/ds_10.cdf')

#%%

results_dir = pjoin(sp_dir, 'Data Share', 'Photoionization', 'Lecroy Resample Check', 'munged_100', '2023-05-18', 'Lecroy')
fp_raw = pjoin(results_dir, 'ds_pos_1.cdf')

ds_raw = xr.load_dataset(fp_raw)

ds_raw = ds_raw.sel(acq_time=tw_goldi)

ds_raw.pint.dequantify().to_netcdf('output/ds_100.cdf')

#%%

results_dir = pjoin(sp_dir, 'Data Share', 'Photoionization', 'Lecroy Resample Check', 'piecewise_10_1000')
fp_raw = pjoin(results_dir, 'ds_pos_1.cdf')

ds_raw = xr.load_dataset(fp_raw)

ds_raw = ds_raw.sel(acq_time=tw_goldi)

ds_raw.pint.dequantify().to_netcdf('output/ds_pw_10_1000.cdf')

#%%

fp_rs = pjoin(REPO_DIR, 'experiment/data/munged/2023-05-18/Lecroy/ds_pos_1.cdf')

ds_rs = xr.load_dataset(fp_rs)

ds_rs_sel = ds_rs.sel(acq_time=tw_goldi)


ds_rs_sel.pint.dequantify().to_netcdf('output/ds_1000.cdf')

#%%

ds_rs_sel.mean('acq_time')['AS'].plot()