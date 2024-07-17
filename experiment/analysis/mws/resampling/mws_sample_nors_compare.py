"""
Comparing mws data between processed and raw unsampled. 

Ran a test case sequence (2023-05-18, 536_pos_1) through the lecroy munging without resampling 
"""

#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis.mws.fitting import calc_dnedt
from scipy.integrate import solve_ivp

#%%

ds_lecroy = ppu.fileio.load_lecroy('536_pos')

ds_final = ds_lecroy.sel(run=('2023-05-18',1)).sel(motor=180, method='nearest')
mag_0 = ppu.fileio.load_mws_T0()


# %%


# ds_final['AS'].mean('mnum').plot()

# plt.yscale('log')

# %%

dss = []

for tc in ['raw', 'pw_10_1000']:
    fp_in = pjoin('output', 'ds_{}.cdf'.format(tc))

    ds = xr.load_dataset(fp_in)
    ds = ds.assign_coords(tc=tc)

    ds.coords['time'] = ds.coords['time'].pint.quantify('s').pint.to('us')
    ds = ds.mws.calc_AS_rel()

    dss.append(ds)

#%%

for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS_rel'].mean('acq_time').plot(label=tc)

plt.legend()

plt.yscale('log')

plt.ylim(1e-3, )
plt.xlim(-1,20)

#%%


fig, axes = plt.subplots(1,2, figsize=(10,5))

for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS_rel'].mean('acq_time').plot(label=tc, ax=axes[0])


axes[0].set_xlim(0,20)
axes[0].legend(['No Resample', '10 peak/1000 off peak'])


for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS_rel'].mean('acq_time').plot(label=tc, ax=axes[1])

for ax in axes:
    ax.set_title('')

axes[1].set_xlim(0.8,2.1)
axes[1].set_ylim(0.03,0.2)

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_sample_nors_compare.png'))
# %%
