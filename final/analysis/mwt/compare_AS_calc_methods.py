#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_mwt_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_lecroy.cdf')).xr_utils.stack_run()
ds_mwt_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')

#%%

ds_sel = ds_mwt_53x.sel(kwt=1, method='nearest').sel(run=('2023-05-18', 1))

ds_sel = ds_sel.pint.quantify()

ds_sel

#%%

from pi_paper_utils.constants import MAG_MWT_BLOCK
from pi_paper_utils.fileio import load_mwt_T0



mag_0 = load_mwt_T0().sel(date=('2023-05-18'))
time_window_delta=slice(-2,-1)
mag_block=MAG_MWT_BLOCK


# Square
ds_sel['T_abs'] = (ds_sel['mag']**2-mag_block**2)/(mag_0**2 - mag_block**2)

ds_sel['T_abs'].attrs = dict(long_name='Absolute Transmission')

ds_sel['AS_abs'] = 1 - ds_sel['T_abs']

AS_pre_pulse_avg = ds_sel['AS_abs'].sel(time=time_window_delta).mean('time')
ds_sel['dAS_abs'] =  ds_sel['AS_abs'] - AS_pre_pulse_avg
ds_sel['dAS_abs'].attrs = dict(long_name='$\Delta AS_{abs}$', units='dimensionless') 

ds_square = ds_sel.copy()


# Linear
ds_sel['T_abs'] = (ds_sel['mag']-mag_block)/(mag_0 - mag_block)

ds_sel['T_abs'].attrs = dict(long_name='Absolute Transmission')

ds_sel['AS_abs'] = 1 - ds_sel['T_abs']

AS_pre_pulse_avg = ds_sel['AS_abs'].sel(time=time_window_delta).mean('time')
ds_sel['dAS_abs'] =  ds_sel['AS_abs'] - AS_pre_pulse_avg
ds_sel['dAS_abs'].attrs = dict(long_name='$\Delta AS_{abs}$', units='dimensionless') 

ds_lin = ds_sel.copy()

#%%

plt.figure()

ds_square['T_abs'].plot(label='Square')
ds_lin['T_abs'].plot(label='Linear')

plt.legend()

plt.xlim(-1,30)

plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, '536_T_lin_square_compare.png'), dpi=300, bbox_inches='tight')

#%%

plt.figure()

ds_square['AS_abs'].plot(label='Square')
ds_lin['AS_abs'].plot(label='Linear')

plt.legend()
# %%

plt.figure()

ds_square['dAS_abs'].plot()
ds_lin['dAS_abs'].plot()

# %%

plt.figure()

ds_square['dAS_abs'].plot(label='Square')
ds_lin['dAS_abs'].plot(label='Linear')  

plt.legend()

plt.yscale('log')

plt.xlim(-1,30)
plt.ylim(1e-4,1e0)
plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, '536_dAS_lin_square_compare.png'), dpi=300, bbox_inches='tight')



# %%
