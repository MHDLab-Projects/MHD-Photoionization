
#%%[markdown]

# exploring dnedt calculation manual parameters
#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdpy.analysis.mws.fitting import calc_dnedt_v2
from scipy.integrate import solve_ivp

ds_lecroy = ppu.fileio.load_lecroy('53x')
da_sel = ds_lecroy.mws.calc_AS_rel()['AS'].sel(kwt=1, method='nearest')
da_fit = da_sel.mean('mnum')
da_fit = da_fit/da_fit.mws._pulse_max()

# %%

t_eval = np.linspace(0, 100, 1000)

dne = Quantity(1e14, 'cm**-3').to('um**-3').magnitude
krm = (1/Quantity(7, 'us')).magnitude
krb = Quantity(1e-8, 'cm**3/s').to('um**3/us').magnitude

sol = solve_ivp(lambda t, y: calc_dnedt_v2(t, y, krm, krb), t_span=(min(t_eval),max(t_eval)), y0 = [dne], t_eval=t_eval)
y = sol.y[0]
y = y/max(y)

# plt.plot(t_eval, y, color='black', linestyle='--')
# plt.yscale('log')

# plt.figure(figsize=(5,5))

fig, axes = plt.subplots(2,1, figsize=(5,5), sharex=False)

plt.sca(axes[0])
da_fit.plot(hue='run_plot', x='time')
plt.plot(t_eval, y, color='black', linewidth=2)

plt.gca().get_legend().remove()

plt.yscale('log')
plt.xlim(-1,5)
plt.ylim(1e-1, 1.1)

plt.title('krm: {:.1e} us^-1, krb: {:.1e} um^3/us, dne: {:.1e}'.format(krm, krb, dne))

plt.sca(axes[1])
da_fit.plot(hue='run_plot', x='time')
plt.plot(t_eval, y, color='black', linewidth=2)

plt.yscale('log')

plt.xlim(-1,30)
plt.ylim(1e-3, 1.1)
plt.gca().get_legend().remove()
plt.title('')


#%%

da_fit = da_sel.copy().mean('mnum')

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3

pipe_fit_mws_3.params['A'].vary = True

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_3(da_fit)

da_fit_norm = da_fit/da_fit.mws._pulse_max()

# ds_mws_fit = xr.merge([ds_mws_fit, da_fit_norm.rename('AS_all')])

#%%

ds_mws_fit.to_array('var').plot(x='time', row='run', hue='var')

plt.yscale('log')

#%%

ds_p['A']


#%%

plt.figure(figsize=(4,3))

ds = ds_mws_fit.sel(run=('2023-05-24', 1))
ds_p_sel = ds_p.sel(run=('2023-05-24', 1))
ds_p_stderr_sel = ds_p_stderr.sel(run=('2023-05-24', 1))

ds['AS_all'].plot(label='Data (all)')
ds['AS_sel'].plot(label='Data (fitted)')
ds['AS_fit'].plot(label='Fit', color='black', linestyle='--')


plt.yscale('log')
plt.ylabel('AS')

plt.xlim(-1,50)
plt.ylim(0.5e-3,2)

kr = ds_p_sel['krm'].item()
kr_stderr = ds_p_stderr_sel['krm'].item()
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

da_fit = da_sel.copy().mean('mnum')

pipe_fit_exp.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(5,'us'),Quantity(30,'us'))

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

#%%

plt.figure()

ds = ds_mws_fit.sel(run=('2023-05-24', 1))
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

pipe_fit_mws_3.perform_fit_kwargs['fit_timewindow'] = slice(Quantity(0,'us'),Quantity(30,'us'))

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_3(da_fit)

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

plt.xlim(-10,50)
plt.ylim(1e9)


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



# %%
