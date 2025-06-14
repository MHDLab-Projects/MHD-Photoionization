#%%

from mhdlab.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

from mhdlab.analysis import aas

#%%

import pytz
import datetime

fp_calib = pjoin(DIR_EXPT_PROC_DATA, 'ds_calib.cdf')

ds_calib = xr.load_dataset(fp_calib)



timestamps = ds_calib['time'].values
timestamps = [str(ts) for ts in timestamps]
# Create a timezone object for PST
pst = pytz.timezone('America/Los_Angeles')

# Parse the timestamps, convert to PST, and format as HH:MM:SS
times_pst = [datetime.datetime.strptime(ts[:26], '%Y-%m-%dT%H:%M:%S.%f').replace(tzinfo=pytz.UTC).astimezone(pst).strftime('%H:%M:%S') for ts in timestamps]

ts = ds_calib['time'].to_series()
ds_calib = ds_calib.assign_coords(date = ('time', ts.apply(lambda x: x.date())))
ds_calib = ds_calib.assign_coords(time_pst= ('time', times_pst))


ds_calib_plot = ds_calib.set_index(time = ('time_pst', 'date')).unstack('time')

dates = set(ds_calib.coords['date'].values)
dates = sorted(dates)

ds_calib_plot

#%%

fig, axes = plt.subplots(len(dates),2, sharex=True, figsize=(10,10))

# for date, ds in ds_calib.groupby('date'):
for i, date in enumerate(dates):
    ds = ds_calib_plot.unstack().sel(date=date).dropna('time_pst', how='all')

    ax = axes[i,0]
    da = ds['diff'].sel(mp='barrel')
    da.plot(hue='time_pst', ax=ax)
    ax.get_legend().remove()
    ax.set_xlabel('')


    ax = axes[i,1]

    if 'mw_horns' in ds.coords['mp'].values:
        da = ds['diff'].sel(mp='mw_horns')
        da.plot(hue='time_pst', ax=ax)
    # Format the legend labels to only show time
    # ax.get_legend().remove()
    ax.set_xlabel('')
    # ax.set_title(date)

plt.savefig(pjoin(DIR_FIG_OUT, 'aas_calib_dates.png'))

#%%



#%%

ds_calib = ds_calib.reset_coords('time_pst' ,drop=True).set_index(dt = ('time', 'date'))

#%%

da = ds_calib['diff'].mean('wavelength')

fig, axes = plt.subplots(1,2, sharex=True, figsize=(10,5))

for date, ds in ds_calib.groupby('date'):
    da_sel = da.sel(date=date)
    da_sel = da_sel.assign_coords(time = da_sel.coords['time'] - da_sel.coords['time'].min())
    # da_sel.plot(hue='mp', col='mp')
    da_sel.sel(mp='barrel').plot(ax=axes[0], marker='o', label=date)
    if 'mw_horns' in da_sel.coords['mp'].values:
        da_sel.sel(mp='mw_horns').plot(ax=axes[1], marker='o', label='date')
    # da_sel.plot(hue='mp', label=date)

axes[0].legend()
#%%

fp = pjoin(REPO_DIR, 'experiment','data','proc_data','ds_aas.cdf')

ds = xr.load_dataset(fp)

ds['diff'] = ds['led_on'] - ds['led_off']

ds = ds.set_index(acq=['time','mp']).unstack('acq')


# assign a new coordinate to ds that is the date of the time coordinate. ds is an xarray Dataset

ts = ds['time'].to_series()
ds['date'] = ts.apply(lambda x: x.date())
ds = ds.set_index(dt = ('time', 'date'), append=True)

ds
# %%

ds

# %%

ds_sel =ds.sel(mp='barrel').dropna('dt', how='all')

ds_sel = ds_sel.unstack('dt')

# make an axes for each date

fig, axes = plt.subplots(len(ds_sel.coords['date'].values),1, sharex=False, figsize=(10,10))

for i, date in enumerate(ds_sel.coords['date'].values):

    da = ds_sel.sel(date=date)['diff'].dropna('time', how='all')
    da.plot(x='time', ax=axes[i], robust=True)
    axes[i].set_title(date)
    # axes[i].get_legend().remove()


# g = ds_sel['led_on'].plot(col='date', sharex=False)



# %%

# bin the wavelengths into 10nm bins

#TODO: get all dates working
date_sel = pd.Timestamp('2023-04-07').date()

bins = [
    725, 735, 745, 755, 765, 775
]

ds2 = ds_sel.groupby_bins('wavelength', bins=bins).mean('wavelength')

ds2 = ds2.sel(date=date_sel)


ds2[['diff', 'calib']].to_array('var').plot(row='wavelength_bins', hue='var')

# plt.yscale('log')
# %%

rat = ds2['diff']/ds2['calib']

rat.plot(row='wavelength_bins')

plt.ylim(0.8,1.2)


#%%

rat.sel(wavelength_bins=pd.Interval(725,735)).plot()

# type(rat.coords['wavelength_bins'].values[0])

# rat.coords['wavelength_bins'].values[0].right



#%%

ds_aas = ppu.fileio.load_aas('53x', avg_mnum=True)

ds_aas = ds_aas.sel(mp='barrel')

dss_aas_methods = [ds_aas['alpha'].assign_coords(method='raw')]

# %%


da = ds_aas['alpha']

g = da.plot(hue='run_plot', x='wavelength', row='kwt')

plt.ylim(-0.1,1.1)

# plt.xlim(760,780)

for ax in g.axes.flatten():
    ax.axhline(0, color='k', linestyle='--')
    ax.axhline(1, color='k', linestyle='--')


# %%
from mhdlab.analysis.aas.fitting import gen_model_alpha_blurred
from mhdlab import analysis
final_model, pars = gen_model_alpha_blurred()


sl1 = slice(750,755)
sl2 = slice(780,785)

beta = -np.log(1-ds_aas['alpha'])/pars['L'].value
beta_off = beta.sel(wavelength=sl1).mean('wavelength')
alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

dss_aas_methods.append(alpha_tc.assign_coords(method='const'))


#%%

da = alpha_tc

g = da.plot(hue='run_plot', x='wavelength', row='kwt')

plt.ylim(-0.1,1.1)
# plt.xlim(760,780)

for ax in g.axes.flatten():
    ax.axhline(0, color='k', linestyle='--')
    ax.axhline(1, color='k', linestyle='--')


#%%

ds_sel = ds_aas

beta = -np.log(1-ds_sel['alpha'])/pars['L'].value

# Make a linear interpolation of beta between the two slices

beta1 = beta.sel(wavelength=sl1).mean('wavelength')
beta2 = beta.sel(wavelength=sl2).mean('wavelength')

beta1_wl_center = beta.sel(wavelength=sl1).coords['wavelength'].mean()
beta2_wl_center = beta.sel(wavelength=sl2).coords['wavelength'].mean()

beta1 = beta1.assign_coords(wavelength=beta1_wl_center)
beta2 = beta2.assign_coords(wavelength=beta2_wl_center)

beta_temp = xr.concat([beta1, beta2], 'wavelength')

beta_temp = beta_temp.sortby('wavelength')

beta_off = beta_temp.interp(wavelength=beta.coords['wavelength'], kwargs={'fill_value':'extrapolate'})

alpha_tc = 1 - np.exp(-(beta - beta_off)*pars['L'].value)

dss_aas_methods.append(alpha_tc.assign_coords(method='linear'))

#%%

beta_off.plot(hue='run_plot', x='wavelength', row='kwt')

# %%
da = alpha_tc#.sel(mp='barrel')

g = da.plot(hue='run_plot', x='wavelength', row='kwt')


plt.ylim(-0.1,1.1)
# plt.xlim(760,780)

for ax in g.axes.flatten():
    ax.axhline(0, color='k', linestyle='--')
    ax.axhline(1, color='k', linestyle='--')
# %%

# put in same form as previous code
ds_fit = xr.concat(dss_aas_methods, 'method').to_dataset(name='alpha')
ds_fit['alpha_red'] = ds_fit['alpha']

model, params = aas.gen_model_alpha_blurred(assert_xs_equal_spacing=False)

ds_alpha_fit, ds_p, ds_p_stderr = ds_fit['alpha_red'].aas.perform_fit(model, params)
# %%


ds_p['nK'].plot(row='run', hue='method', x='kwt', marker='o')
# %%
