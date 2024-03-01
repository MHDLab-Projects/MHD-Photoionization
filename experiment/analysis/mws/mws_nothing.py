#%%[markdown]

# # Investigating mws data with no sample. Measurement of air transmission T0

#%%

from mhdpy.analysis.standard_import import *
from mhdpy.coords import gen_coords_to_assign_1, assign_coords_multi
import pi_paper_utils as ppu

#have time axes show up in PST
plt.rcParams['timezone'] = 'US/Pacific'

# %%

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()

da_motor = dsst['motor']['Motor C Relative'].rename(time='acq_time')
da_motor    

# fp_dst = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
# dst = mhdpy.fileio.TFxr(fp_dst).as_dsst()['coords']

# da_motor_coords = dst['motor'].rename(time='acq_time')

dst_coords = gen_coords_to_assign_1(dsst)

da_motor_coords = dst_coords['motor'].rename(time='acq_time')
da_motor_coords = da_motor_coords.rename('motor')

#%%

dir = pjoin(REPO_DIR, 'experiment', 'data','munged')

dates = {
    '2023-04-07': ['ds_Nothing_init_mp.cdf'],
    # '2023-05-12': ['ds_Shutdown.cdf'], # I believe this was taken with the jet still on. #TODO: check 
    '2023-05-12': ['ds_Init_Nothing.cdf'],
    '2023-05-18': ['ds_nothing_after.cdf', 'ds_Nothing.cdf', 'ds_Nothing_2.cdf', 'ds_Nothing_3.cdf'],
    '2023-05-24': ['ds_nothing_2.cdf', 'ds_Nothing_Init.cdf']
}

#%%

# Examine time trace of a test case
# time coords are different between files as data was taken near silicon ( which requires a larger timebase). TODO: align time axes below and join

fn = dates['2023-05-18'][1]

fp = pjoin(dir, '2023-05-18', 'Lecroy', fn)

ds = xr.load_dataset(fp)

ds.coords['time'].attrs['units'] = 'us'

ds = ds.mws.calc_mag_phase_AS()

ds['mag'].mean('acq_time').plot()

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_nothing_time_trace_2023-05-18.png'))


#%%

dss = []

import re

for date, fn_list in dates.items():
    for fn in fn_list:
        fp = pjoin(dir, date, 'Lecroy', fn)
        ds = xr.load_dataset(fp)
        time_offset = 0.93
        ds = ds.assign_coords(time=ds.coords['time']*1e6 - time_offset)
        ds.coords['time'].attrs['units'] = 'us'
        ds.coords['time'].attrs['long_name'] = 'Time'

        tc = re.search(r'ds_(.*).cdf', fn).group(1)
        # ds = ds.mean('acq_time')
        ds = ds.assign_coords(date=date).expand_dims('date')
        ds = ds.assign_coords(tc=tc).expand_dims('tc')
        ds = ds.stack(temp=['date','tc'])

        ds = ds.mws.calc_mag_phase_AS()
        ds = ds.mean('time')

        # motor = da_motor.interp(acq_time=ds['acq_time'], method='nearest')
        motor = da_motor_coords.interp(acq_time=ds['acq_time'], method='nearest')
        ds['motor'] = motor

        dss.append(ds)


# ds = ds.mws.calc_mag_phase_AS()

ds_concat = xr.concat(dss, dim='temp')

#%%

ds_concat


#%%

# ds_concat['acq_time'].plot(row='temp')

# ds_concat['acq_time']

fig, axes = plt.subplots(len(ds_concat['temp']), 2, figsize=(10,15))

for i, (temp, ds) in enumerate(ds_concat.groupby('temp')):
    ax = axes[i, 0]
    ds['mag'].plot(ax=ax)
    ax.set_title(temp)

    ax = axes[i, 1]
    ds['motor'].plot(ax=ax)


#%%

# Drop test cases where motor was not at goldilocks position

from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

ds_sel = ds_concat#.unstack('temp').sel(date='2023-05-24').dropna('tc', how='all')

ds_2 = assign_signal(ds_sel, da_motor_coords.rename(acq_time='time'), 'acq_time')

ds_2 = unstack_multindexed_acq_dim(ds_2, 'acq_time')

ds_2['mag'].mean('mnum').plot(row='temp', marker='o')


#%%


ds_sel = ds_2.unstack('temp').sel(date='2023-05-24').dropna('tc', how='all')

ds_sel

ds_sel.mean('mnum')['mag'].plot(hue='tc', marker='o')

plt.ylabel("U (V)")
plt.xlabel("Stage Postion (mm)")

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_nothing_motor_2023-05-24.png'))


#%%

#%%

fig, axes = plt.subplot_mosaic([[date for date in dates.keys()]], figsize=(15,10))

for ds in dss:

    date = ds['date'].item()
    tc = ds['tc'].item()

    ax = axes[date]
    ax.plot(ds['acq_time'], label=tc)


for ax in axes.values():
    ax.legend()


#%%

for ds in dss:

    date = ds['date'].item()
    tc = ds['tc'].item()

    ds_offset = ds#.mean('time')
    ds_offset['acq_time'] = ds_offset['acq_time'] - ds_offset['acq_time'].min()

    label = f'{date} {tc}'

    plt.plot(ds_offset['acq_time'], ds_offset['mag'], label=label)

plt.legend()

#%%


ds_sel = ds_2.sel(motor=178, method='nearest').mean('mnum').dropna('temp', how='all')

ds_sel = ds_sel.unstack('temp')

ds_sel['mag']

#%%

m = ds_sel['mag'].groupby('date').mean('tc')
std = ds_sel['mag'].groupby('date').std('tc')

plt.errorbar(m['date'], m, yerr=std, fmt='o')

#%%

# Need to make a correction factor for 2023-05-12 because the Nothing data was only taken at mp = 35. At this location it appears there is a slight distortion, probably from torch. Measure this factor for other dates and apply to get 

ds35 = ds_2.mean('mnum').sel(motor=35, method='nearest')

ds178 = ds_2.mean('mnum').sel(motor=178, method='nearest')

rat = ds178/ds35

corr_factor = rat['mag'].mean('temp').item()

corr_factor


#%%

# Calculate value for 178 motor position for 05-12

ds_35_0512 = ds_2.sel(temp=('2023-05-12', 'Init_Nothing')).sel(motor=35, method='nearest').mean('mnum')

val_178_05_12 = ds_35_0512['mag']*corr_factor

val_178_05_12 = val_178_05_12.item()

val_178_05_12


#%%

s = m.to_pandas()
s.loc['2023-05-12'] = val_178_05_12

s = s.sort_index()

s.to_csv(pjoin(REPO_DIR, 'experiment', 'analysis', 'mws', 'output', 'mws_T0.csv'))

#%%

s.plot(marker='o')


# %%

fp = r'/home/leeaspitarte/code/MHD-Photoionization/experiment/data/munged/2023-05-12/Lecroy/ds_savingtest_15hz.cdf'

ds = xr.load_dataset(fp)

ds.coords['time'].attrs['units'] = 'us'

ds = ds.mws.calc_mag_phase_AS()

ds['mag'].mean('acq_time').plot()

#%%

tc = '536_pos'

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()#[['mag', 'phase','AS']]


#%%

fp_nothing = pjoin(REPO_DIR, 'experiment','analysis','mws','output','mws_T0.csv')

df_nothing = pd.read_csv(fp_nothing, index_col=0)['0']
da_nothing = xr.DataArray(df_nothing).pint.quantify('volt')

da_nothing

#%%

ds_lecroy = ds_lecroy.unstack('run').mws.calc_mag_phase_AS(mag_0=da_nothing)#[['mag', 'phase','AS']]

#%%
ds_lecroy


# %%

da_sel = ds_lecroy['mag_pp'].mean('mnum')

# da_sel.plot(hue='run_plot', x='motor', marker='o')

norm = da_sel.sel(motor=slice(200,300)).mean('motor')

#%%

da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')

da_sel = da_sel/norm

da_sel = da_sel.dropna('run', how='all')

da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))

#%%[markdown]

# Ignore 2023-05-12. No explicity T0 calibration 



#%%

da_sel = ds_lecroy['mag_pp'].mean('mnum')

da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()

da_sel.plot(hue='run_plot', x='motor', marker='o')

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_nothing_motor_2023-05-12.png'))

#%%

da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')

#%%

# da_nothing
da_sel.unstack('run')#.mws.calc_mag_phase_AS(mag_0=da_nothing)

#%%


da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()
da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))


# %%

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
# Display img on the upper left subplot

da_nothing.to_series().plot(ax=axes[0], marker='o')

axes[0].set_title('T No Torch')

plt.ylabel("$U_{Nothing} (V)$")
plt.xlabel("Date")

da_sel = ds_lecroy['mag_pp'].mean('mnum')
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()
da_sel.plot(hue='run_plot', x='motor', marker='o', ax=axes[1])

da_sel = ds_lecroy['mag'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')
da_sel = (da_sel.unstack('run')/da_nothing).xr_utils.stack_run()

axes[1].set_title('$U/U_{Nothing}$')
axes[1].set_xlabel("Motor Position (mm)")

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_nothing_T0.png'))
