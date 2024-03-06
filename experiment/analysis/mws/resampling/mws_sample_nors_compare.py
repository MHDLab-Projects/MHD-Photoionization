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

# %%


ds_final['AS'].mean('mnum').plot()

plt.yscale('log')

# %%

dss = []

for tc in ['raw', 'pw_10_1000']:
    fp_in = pjoin('output', 'ds_{}.cdf'.format(tc))

    ds = xr.load_dataset(fp_in)
    ds = ds.assign_coords(tc=tc)

    ds.coords['time'] = ds.coords['time'].pint.quantify('s').pint.to('us')
    ds = ds.mws.calc_mag_phase_AS()

    dss.append(ds)

#%%

for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS'].mean('acq_time').plot(label=tc)

plt.legend()

plt.yscale('log')

plt.ylim(1e-3, )
plt.xlim(-1,20)

#%%


for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS'].mean('acq_time').plot(label=tc)


plt.xlim(0,2)
plt.legend()

#%%



for ds in dss:
    tc = ds.coords['tc'].item()
    ds['AS'].mean('acq_time').plot(label=tc)

plt.xlim(0.8,2.1)
plt.ylim(0.03,0.1)
plt.legend()
# %%
