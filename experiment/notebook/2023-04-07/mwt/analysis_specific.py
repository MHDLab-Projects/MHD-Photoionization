#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import re

datestr = '2023-04-07'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')


dsst = mhdlab.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)
# #%%


tc = '536_power'

from mhdlab.analysis import mwt


fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

ds_in = xr.load_dataset(fp_in)
ds_in = ds_in.sel(date=datestr).sel(run_num=1)

ds = ds_in.mwt.calc_AS_rel().drop('mag_pp')


tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

mnum_counts
#%%

ds_avg = ds.mean('mnum', keep_attrs=True)

# %%

from mhdlab.analysis.mwt.fitting import pipe_fit_mwt_1

da_fit = ds_avg['AS'].copy()

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_1(da_fit)

#%%

plt.rcParams.update({'font.size': 14})


plt.figure()


tc_dim = 'power'
tc_val = 1


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

ds_mwt_fit['AS'].sel({tc_dim: tc_val}).plot(marker='o', label='Data')
ds_mwt_fit['AS_fit'].sel({tc_dim: tc_val}).plot(linestyle='--', label='Fit', color='red')

plt.yscale('log')
plt.ylim(5e-4,2)
plt.xlim(-0.5,30)

plt.legend()


plt.title('')
plt.ylabel("Norm. Absorption + Scattering")

plt.text(0, 0.8e-3, param_str)

# if i != len(tc_vals) -1: ax.set_xlabel('')

# if i != 0:
#     ax.set_title('')


plt.savefig(pjoin(DIR_FIG_OUT,'fit_specific.png'))
# %%

tw_power = slice(Timestamp('2023-04-07 15:41:23.104546816'), Timestamp('2023-04-07 15:50:16.422889984'), None)

da_power_cal = dsst['lasen_meter2']['Power']

print("Laser Power Calibration")

print("Mean: {}, Std: {}", 
      da_power_cal.sel(time=tw_power).mean('time').item(),
      da_power_cal.sel(time=tw_power).std('time').item()
)


#%%

from pint import Quantity

mwt_offset_min = Quantity(1 + 2/16, 'inches')
mwt_offset_max = Quantity(1 + 4/16, 'inches')
pos_536_goldilocks = Quantity(150, 'mm')

print("SFR-maximized position min: {}".format((mwt_offset_min + pos_536_goldilocks).to('mm')))
print("SFR-maximized position max: {}".format((mwt_offset_max + pos_536_goldilocks).to('mm')))
