#%%
from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

ds_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')

ds_536 = ds_536.pint.quantify()

ds_536 =ds_536.mws.calc_time_stats()

# %%

da_AS_pp = ds_536['AS_abs'].sel(time=slice(-50,-1))

da_AS_pp['time'] = da_AS_pp['time'].pint.to('s')

da_AS_pp = da_AS_pp.pint.dequantify()

da_AS_pp
# %%

import xrft

# da_fft = xrft.fft(da_AS_pp, dim='time')
da_fft = xrft.power_spectrum(da_AS_pp, dim='time', scaling='density')

motor_sel = [30,100,180]

da_fft = da_fft.mean('run').mean('mnum').sel(motor=motor_sel, method='nearest')#.plot(hue='motor')


da_fft
# %%

da_fft.real.plot(hue='motor', x='freq_time')

plt.xscale('log')
# plt.xlim(0.1,1)

plt.yscale('log')

plt.ylabel("Power Spectral Density")

plt.title('')

plt.xlabel('Frequency (Hz)')

plt.savefig(pjoin(DIR_FIG_OUT, 'AS_fft.png'))



# %%
