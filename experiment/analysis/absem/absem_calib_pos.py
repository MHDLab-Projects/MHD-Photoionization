
#%%[markdown]

# # Examing Calibration of ABSEM signal with motor position

# Calibration was taken with position dependence on 2023-05-24, and an extra datapoint on 2023-05-18. Looking at the dependence 
# The positino dependence doesn't appear to be reproducible, and is on the order of 1%. Just using calibration at the goldilocks position of 150mm for now. 

#%%
from mhdlab.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

from mhdlab.fileio.tdms import TFxr
from mhdlab.fileio.spectral import load_absem


# datestr = '2023-05-24'
# data_folder = os.path.join(REPO_DIR, 'experiment', 'data' , 'munged',datestr)
# dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

fp_dsst = os.path.join(DIR_EXPT_PROC_DATA, 'dsst.tdms')

dsst = TFxr(fp_dsst).as_dsst(convert_to_PT=False)

fp = os.path.join(DIR_EXPT_PROC_DATA, 'ds_absem.cdf')
ds_absem = xr.load_dataset(fp)
ds_absem = ds_absem.set_index(acq=['time','mp']).unstack('acq')
# %%

tws = {
    '2023-05-18_bef': slice(Timestamp('2023-05-18 20:38:07.167729408'), Timestamp('2023-05-18 20:46:11.524545024'), None),
    '2023-05-24_bef': slice(Timestamp('2023-05-24 18:06:29.319461376'), Timestamp('2023-05-24 18:15:49.467337728'), None),
    '2023-05-24_aft': slice(Timestamp('2023-05-24 22:35:43.373361152'), Timestamp('2023-05-24 22:44:56.927864320'), None),

}


from mhdlab.coords import assign_signal, unstack_multindexed_acq_dim
motor = dsst['motor']['Motor C Relative'].rename('motor')

dss = []

for tw_calib_pos, tw in tws.items():
    ds_sel = ds_absem.sel(time=tw).sel(mp='mw_horns')
    ds_sel = assign_signal(ds_sel, motor, 'time')
    ds_sel = unstack_multindexed_acq_dim(ds_sel, min_mnum=2).dropna('motor', how='all')

    ds_sel = ds_sel.assign_coords(tc=tw_calib_pos)

    dss.append(ds_sel)

ds_sel = xr.concat(dss, 'tc')

#%%

# ds_sel.isel(wavelength=0).groupby('motor').count('mnum')

da_led = ds_sel['led_on'].mean('wavelength').mean('mnum')

da_led_norm = da_led/da_led.max('motor')

g= da_led_norm.plot(hue='tc', marker='o')

dropna(g)

plt.ylim(0.97,)
# plt.ylim(550,)