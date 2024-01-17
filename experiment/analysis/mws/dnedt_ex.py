
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
da_fit = da_sel.mean('mnum')
da_fit = da_fit/da_fit.mws._pulse_max()

# sol = solve_ivp(lambda t, y: calc_dnedt(t, y, ne0, kr), t_span=(min(t_eval),max(t_eval)), y0 = [dne], t_eval=t_eval)


# %%

t_eval = da_fit.sel(time=slice(0,30)).time.values
dne=10
ne0 = 2
kr = .05

sol = solve_ivp(lambda t, y: calc_dnedt(t, y, ne0, kr), t_span=(min(t_eval),max(t_eval)), y0 = [dne], t_eval=t_eval)

sol

#%%

y = sol.y[0]
y = y/max(y)

plt.plot(t_eval, y)

plt.yscale('log')

da_fit.plot(hue='run_plot', x='time')
# da_fit.plot(hue='run_plot')

plt.xlim(-1,30)

#%%


