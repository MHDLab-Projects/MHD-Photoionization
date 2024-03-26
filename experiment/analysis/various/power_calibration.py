"""
Power calibration for the laser power meters.
Two measurements: full power on both pyroelectric and photovoltaic power meters to calibrate PV meter, and a sweep of the filter wheel while measuring PE.
the PV power meter was used throughout the experiment, whereas the PE power meter was only used for calibration before and after experiment
"""

#%%

from mhdpy.analysis.standard_import import *
from mhdpy.fileio.tdms import TFxr
from mhdpy.coords.ct import downselect_acq_time
dir_dataset = pjoin(REPO_DIR, 'experiment','data', 'proc_data')

dsst = TFxr(pjoin(dir_dataset,'dsst.tdms')).as_dsst()

dsst

#%%

pow_pv = dsst['lasen_meter1']['Power']
pow_pe = dsst['lasen_meter2']['Power']
fw = dsst['filterwheel']['Filter Position']

#%%
from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list

fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_power_calibration.csv')

df_cuttimes = load_df_cuttimes(fp_cuttimes)

df_cuttimes
# %%

df_fullpower_times = df_cuttimes[df_cuttimes['Event'].str.contains('fullpower')]

df_fullpower_times

pow_pe_fullpower = downselect_acq_time(pow_pe, df_fullpower_times, timeindex='time')


#%%
pow_pe_fullpower = pow_pe_fullpower.assign_coords(date = ('time', pow_pe_fullpower['time'].dt.date.values))

pow2 = pow_pe_fullpower.set_index(time='date', append=True)

pow_mean = pow2.groupby('date').mean('time')
pow_std = pow2.groupby('date').std('time')

#TODO: the time windows on 2023-05-18 in UTC cause this to regiser as 05-19. Hack to get right date to display
datestrs = pow_mean['date'].astype(str).str.replace('2023-05-19', '2023-05-18').values

#%%

plt.errorbar(datestrs, pow_mean, pow_std, capsize=5, marker='o')

plt.ylabel("Pulse Energy (J)")

plt.xticks(rotation=90)

avg_power = pow_pe_fullpower.mean('time').item()

plt.title("Average Energy: {:2f} J".format(avg_power))

plt.savefig(pjoin(DIR_FIG_OUT, 'laser_energy_date.png'))

#%%

avg_power

#%%
pow_mean['date']


#%%

from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

df_sweep_times = df_cuttimes[df_cuttimes['Event'].str.contains('sweep')]

# Needing to downselect time to avoid error in assign_signal. Beleive it is due to filterwheel data not starting at the same time as the power meter data. TODO: check this an add assertion to assign_signal
tw_sweeps = slice(df_sweep_times['Start Time'].iloc[0], df_sweep_times['Stop Time'].iloc[-1])

pow_pe_sweep = assign_signal(pow_pe.sel(time=tw_sweeps), fw, 'time')

# Downselect to only the sweep times. 
pow_pe_sweep = downselect_acq_time(pow_pe_sweep, df_sweep_times, timeindex='time')

pow_pe_sweep = unstack_multindexed_acq_dim(pow_pe_sweep, 'time')

pow_pe_sweep = pow_pe_sweep['Power']

pow_pe_sweep

#%%


pow_pe_sweep.mean('mnum').plot(marker='o')

plt.yscale('log')



# %%
