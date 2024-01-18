
#%%[markdown]

# exploring dnedt calculation manual parameters
#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis.mws.fitting import calc_dnedt
from scipy.integrate import solve_ivp

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()

ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS()


da_sel = ds_lecroy['AS'].sel(kwt=1, method='nearest')

# sol = solve_ivp(lambda t, y: calc_dnedt(t, y, ne0, kr), t_span=(min(t_eval),max(t_eval)), y0 = [dne], t_eval=t_eval)


# %%

da_fit = da_sel.mean('mnum')
da_fit = da_fit/da_fit.mws._pulse_max()

t_eval = da_fit.sel(time=slice(0,30)).time.values
dne=10
ne0 = 2
kr = .05

sol = solve_ivp(lambda t, y: calc_dnedt(t, y, ne0, kr), t_span=(min(t_eval),max(t_eval)), y0 = [dne], t_eval=t_eval)

sol

#%%

y = sol.y[0]
y = y/max(y)


plt.yscale('log')

da_fit.plot(hue='run_plot', x='time')
# da_fit.plot(hue='run_plot')
plt.plot(t_eval, y, color='black', linestyle='--')

plt.xlim(-1,30)

#%%

da_fit = da_sel.copy()

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit)

da_fit_norm = da_fit/da_fit.mws._pulse_max()

# ds_mws_fit = xr.merge([ds_mws_fit, da_fit_norm.rename('AS_all')])

#%%

ds_mws_fit.mean('mnum').to_array('var').plot(x='time', row='run', hue='var')

plt.yscale('log')


#%%

ds = ds_mws_fit.sel(run=('2023-05-24', 1)).mean('mnum')
ds_p_sel = ds_p.sel(run=('2023-05-24', 1))
ds_p_stderr_sel = ds_p_stderr.sel(run=('2023-05-24', 1))

ds['AS_all'].plot(label='Data (all)')
ds['AS_sel'].plot(label='Data (fitted)')
ds['AS_fit'].plot(label='Fit', color='black', linestyle='--')

plt.legend()

plt.yscale('log')
plt.ylabel('AS')

plt.xlim(-1,50)

kr = ds_p_sel['kr'].item()
kr_stderr = ds_p_stderr_sel['kr'].item()
dne = ds_p_sel['dne'].item()
dne_stderr = ds_p_stderr_sel['dne'].item()

# add text to plot of kr and dne with uncertianties

kr_str = '$k_r$ = {:.3f} +/- {:.3f}  um^3/us'.format(kr, kr_stderr)
dne_str = '$\Delta n_e$ = {:.2f} +/- {:.2f} #/um^3'.format(dne, dne_stderr)

plt.text(0.2, 0.95, kr_str, transform=plt.gca().transAxes)
plt.text(0.2, 0.90, dne_str, transform=plt.gca().transAxes)

#%%[markdown]

# # check with exponential fit

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

da_fit = da_sel.copy()

da_fit = da_fit.mws.fit_prep()


ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(
    da_fit,
    fit_timewindow=slice(Quantity(5,'us'),Quantity(30,'us')),
    )

#%%

ds = ds_mws_fit.sel(run=('2023-05-24', 1)).mean('mnum')
ds_p_sel = ds_p.sel(run=('2023-05-24', 1))
ds_p_stderr_sel = ds_p_stderr.sel(run=('2023-05-24', 1))

ds['AS_all'].plot(label='Data (all)')
ds['AS_sel'].plot(label='Data (fitted)')    
ds['AS_fit'].plot(label='Fit', color='black', linestyle='--')

plt.legend()

plt.yscale('log')
plt.ylabel('AS')

plt.xlim(-1,50)

decay = ds_p_sel['decay'].item()
decay_stderr = ds_p_stderr_sel['decay'].item()

# add text to plot of kr and dne with uncertianties

decay_str = '$\\tau$ = {:.3f} +/- {:.3f}  us'.format(decay, decay_stderr)

plt.text(0.2, 0.95, decay_str, transform=plt.gca().transAxes)

# %%
