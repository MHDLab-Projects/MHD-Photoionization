#%%[markdown]

# # Lecroy timewindow 

# Examining the time window for the Lecroy data.
#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws

# %%

tc = '53x'
ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

ds_lecroy = ds_lecroy.drop(0,'kwt')

da_sel = ds_lecroy['AS']

#%%

da_fit_prep = da_sel.mws.fit_prep()
# %%

da_fit_prep.mean('mnum').plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')
plt.xlim(-1,50)

#%%

mean = abs(da_fit_prep).mean('mnum')
std = da_fit_prep.std('mnum')
rat = mean/std

#%%

rat.plot(hue='run_plot', row='kwt', x='time')

plt.yscale('log')

plt.xlim(-1,50)

#%%


rat_mean = rat.mean('run')

rat_mean.plot(hue='kwt', x='time')

plt.yscale('log')
plt.xlim(-1,50)



#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2 

dss_p = []

for tw_max in [5,10,20,30,40,50]:
    print('tw_max = {}'.format(tw_max))
    ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(
                                    da_sel, 
                                    fit_timewindow=slice(Quantity(0, 'us'),Quantity(tw_max, 'us')),
                                    take_log=False
                                    )

    dss_p.append(ds_p.assign_coords(tw_max=tw_max))


#%%

ds_p = xr.concat(dss_p, 'tw_max')

#%%

ds_p['kr'].plot(hue='run_plot', row='kwt', x='tw_max', marker='o')

plt.yscale('log')
