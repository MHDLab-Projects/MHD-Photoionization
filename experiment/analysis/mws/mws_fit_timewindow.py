#%%[markdown]

# # Lecroy timewindow 

# Examining the time window for the Lecroy data.
#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu    

from mhdpy.analysis import mws

# %%

ds_lecroy = ppu.fileio.load_lecroy('53x', AS_calc='relative')
da_sel = ds_lecroy['AS']

#%%

da_fit_prep = da_sel.mws.fit_prep(min_mnum=100)
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

da_fit = da_sel.mean('mnum')

for tw_max in [5,10,20,30,40,50]:
    print('tw_max = {}'.format(tw_max))

    pipe_fit_mws_2.perform_fit_kwargs.update({'fit_timewindow':slice(Quantity(0,'us'),Quantity(tw_max,'us'))})
    ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(
                                    da_fit, 
                                    )

    dss_p.append(ds_p.assign_coords(tw_max=tw_max))


#%%

ds_p = xr.concat(dss_p, 'tw_max')

#%%

ds_p['krb'].plot(hue='run_plot', row='kwt', x='tw_max', marker='o')

plt.yscale('log')

# %%
ds_p