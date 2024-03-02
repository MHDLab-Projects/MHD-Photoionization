
#%%[markdown]

# exploring dnedt calculation manual parameters
#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdpy.analysis.mws.fitting import calc_dnedt
from scipy.integrate import solve_ivp

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x')
da_sel = ds_lecroy['AS'].sel(kwt=1, method='nearest')

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

plt.figure(figsize=(4,3))

ds = ds_mws_fit.sel(run=('2023-05-24', 1)).mean('mnum')
ds_p_sel = ds_p.sel(run=('2023-05-24', 1))
ds_p_stderr_sel = ds_p_stderr.sel(run=('2023-05-24', 1))

ds['AS_all'].plot(label='Data (all)')
ds['AS_sel'].plot(label='Data (fitted)')
ds['AS_fit'].plot(label='Fit', color='black', linestyle='--')


plt.yscale('log')
plt.ylabel('AS')

plt.xlim(-1,50)
plt.ylim(0.5e-3,2)

kr = ds_p_sel['kr'].item()
kr_stderr = ds_p_stderr_sel['kr'].item()
dne = ds_p_sel['dne'].item()
dne_stderr = ds_p_stderr_sel['dne'].item()

# add text to plot of kr and dne with uncertianties

kr_str = '$k_r$ = {:.3f} +/- {:.3f}  um^3/us'.format(kr, kr_stderr)
dne_str = '$\Delta n_e$ = {:.2f} +/- {:.2f} #/um^3'.format(dne, dne_stderr)

plt.text(0.2, 0.93, kr_str, transform=plt.gca().transAxes)
plt.text(0.2, 0.85, dne_str, transform=plt.gca().transAxes)

plt.legend(bbox_to_anchor=(1, 0.6), loc='right')
6
plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mws_dnedt.png'), bbox_inches='tight')

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

plt.figure()

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

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mws_exp.png'))

# %%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_3(
    da_fit.mean('mnum'),
    fit_timewindow=slice(Quantity(0,'us'),Quantity(30,'us')),
    )

#%%

plt.figure()

ds = ds_mws_fit.sel(run=('2023-05-24', 1))

dne = ds_p['dne'].pint.quantify('particle/um**3').pint.to('particle/cm**3')
dne = dne.sel(run=('2023-05-24', 1)).item()

ds = ds*dne

ds['AS_all'].plot(label='Data (all)')
ds['AS_sel'].plot(label='Data (fitted)')
ds['AS_fit'].plot(label='Fit', color='black', linestyle='--')

plt.legend(bbox_to_anchor=(1, 0.6), loc='right')

plt.yscale('log')
plt.ylabel('$\Delta n_e (cm^{-3})$')

plt.title('')

plt.xlim(-1,50)
plt.ylim(1e11)


tau = 1/ds_p['krm'].pint.quantify('us**-1')
tau = tau.sel(run=('2023-05-24', 1)).item()



# add text to plot of tau and dne with uncertianties

tau_str = '$\\tau$ = {:.3f}  us'.format(tau.magnitude)
dne_str = '$\Delta n_e$ = {:.2e} #/cm^3'.format(dne.magnitude)

plt.text(0.2, 0.92, tau_str, transform=plt.gca().transAxes, fontsize=12)
plt.text(0.2, 0.85, dne_str, transform=plt.gca().transAxes, fontsize=12)

plt.savefig(pjoin(DIR_FIG_OUT, 'fit_mws_dnedt_v2.png'))

#%%

dne

#%%

n_e_profile = ds['AS_all']*dne
n_e_profile.name = 'delta n_e'

n_e_profile.plot()

plt.yscale('log')


