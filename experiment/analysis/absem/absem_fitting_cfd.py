# %% [markdown]
#  # Absem fitting with CFD profiles

# uses euler method to solve the differential equation for I(x)
# using normalized number density profile nK(x) and maximum number density nK_max

# %%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

# %%

tc = '536_pos'

ds_absem = ppu.fileio.load_absem(tc)

ds_absem

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd_profiles = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_beam_Yeq.cdf')

ds_cfd = xr.load_dataset(fp_cfd_profiles)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd.coords['dist'] = ds_cfd.coords['dist'].pint.to('cm') # Fitting expects cm
ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)

da_cfd = ds_cfd['Yeq_K']

da_cfd



# %%

# make an artifical tophat profile. zero outside of the beam, constant inside the beam 

L = Quantity(1, 'cm')
#TODO: why centered at 5.5?
L_center = Quantity(5.5, 'cm')

def tophat_profile(L, L_center, x):
    nK = np.zeros_like(x)
    nK[(x > L_center - L/2) & (x < L_center + L/2)] = 1
    return nK

# x = np.linspace(0, 2, 100)
x = da_cfd.coords['dist'].values

nK_profile = tophat_profile(L.magnitude, L_center.magnitude, x) 

plt.plot(x, nK_profile)


# %%[markdown]

# # I(x) solution functions
#
# Solve differential equation dI/dx = -mu_atten(x) I(x)
#
# We want to solve the differential equation for I(x) using the normalized number density profile nK_CFD(x).
#
# Given a max nK_max (Fit parameter) we calculate the absorption coefficient mu_atten(x) = kapppa * nK_max * nK_profile(x).
#
# x is the distance along the beam axis.
#
# Below are three methods to solve the differential equation for I(x) using the mu_atten(x) profile.
#
# We first test everything with the tophat profile, and just selecting a single wavelength and nK_max

# %%

from scipy.integrate import solve_ivp
from mhdpy.analysis.absem.fitting import kappa_2peak
from pint import Quantity

from mhdpy.analysis.absem.fitting import calc_I_profile_euler, alpha_deq_solve


wl = Quantity(770, 'nm')

kappa = kappa_2peak(wl)
nK = Quantity(1e21, '1/m^3')

mu_atten_max = (kappa*nK).to('1/cm').magnitude
mu_atten_max_profile = mu_atten_max * nK_profile
mu_atten_max_profile = pd.Series(mu_atten_max_profile, index=x)

I = calc_I_profile_euler(mu_atten_max_profile)

tau = I[-1]/I[0]

print(tau)

line1 = plt.plot(x, I, label='I(x)')
plt.gca().set_ylabel('I(x)')
plt.gca().set_xlabel('x [cm]')
ta=  plt.twinx()
ta.set_ylabel('mu_atten_max(x)')
line2 = plt.plot(x, mu_atten_max_profile, 'r', label='$\mu_{atten}$ (x)')

lines = line1 + line2
labels = [l.get_label() for l in lines]
plt.legend(lines, labels, loc=(0.2,1.05))

plt.xlim(3,8)

plt.savefig(pjoin(DIR_FIG_OUT, 'euler_method_tophat_demo.png'))

#%%

mu_atten_max_profile

# %%[markdown]

# Now pull in the CFD data to make a normalized number density profile nK_CFD(x) that is used to calculate the absorption coefficient mu_atten(x) = kappa * nK(x).


# %%


g = da_cfd.plot(hue='motor', col='kwt')

plt.yscale('log')

for ax in g.axes.flatten():
    ax.plot(x, nK_profile*nK.to('1/cm^3'), 'k--')

plt.ylim(1e13,1e16)

plt.xlim(3,7)

# plt.gca().get_legend().remove()


# %%
da_cfd_sel = da_cfd.sel(kwt=1)

da_cfd_mu_atten = (da_cfd_sel/da_cfd_sel.max())*mu_atten_max

da_cfd_mu_atten 

Is = []

for motor in da_cfd_mu_atten.coords['motor'].values:
    da_sel = da_cfd_mu_atten.sel(motor=motor)
    mu_atten_profile = pd.Series(da_sel.values, index=da_sel.coords['dist'].values)

    I = calc_I_profile_euler(mu_atten_profile)

    Is.append(I)


da_I_pos = xr.DataArray(Is, coords=da_cfd_mu_atten.coords, dims=['motor', 'dist'])


# %%

da_I_pos.plot(hue='motor')

plt.ylabel('I(x)')
plt.xlabel('x_beam [cm]')

plt.gca().get_legend().set_title('Stage position [mm]')

# %%

# Compare to the experimental data

da_sel = ds_absem['alpha'].mean('mnum').mean('run').sel(mp='mw_horns')
da_sel.coords['motor'] = da_sel.coords['motor']

da_sel_wl = da_sel.sel(wavelength=770, method='nearest')
da_sel_max = da_sel.max('wavelength')

da_sel_wl.plot(label='exp 770 nm', marker='o')
da_sel_max.plot(label='exp max', marker='o')


da_tau = da_I_pos[:,-1]/da_I_pos[:,0]

da_alpha = -da_tau + 1

da_alpha.plot(marker='o')

da_alpha.plot(label='cfd 770 nm', marker='o')

plt.legend()

plt.yscale('log')
plt.ylabel('alpha')
plt.xlabel('Stage position [mm]')

# .sel(wavelength=770, method='nearest')

# %%[markdown]

# # Full spectrum fiitting

# Now we perform the above for each individual wavelength in the spectrum

#%%

# wls = ds_absem.coords['wavelength'].values


# %%

from mhdpy.analysis.absem.fit_prep import interp_alpha

ds_test = ds_absem.mean('mnum').mean('run').sel(mp='barrel')

ds_test = ds_test.mean('motor')

ds_test = ds_test.absem.reduce_keep_wings()

# ds_test = ds_test.dropna('wavelength')

ds_test['alpha_red'].plot(marker='o')


# %%

wls = ds_test.coords['wavelength'].values

nK = Quantity(1e-5, '1/nm^3')

x = np.linspace(0, 2, 100)
nK_profile = tophat_profile(L.magnitude, 1, x)
nK_profile = pd.Series(nK_profile, index=x)

alpha = alpha_deq_solve(wls, nK_profile=nK_profile, nK_max=nK)
ds_test['alpha'].plot(marker='o')
plt.plot(wls, alpha)

plt.ylim(0,1)


# %%


# %%

ds_fit = ds_absem.mean('mnum').mean('run').sel(mp='mw_horns').dropna('motor')

ds_fit

#%%

# quick fix for calibraiton

ds_fix = ds_fit.copy()

rat = ds_fit['diff']/ds_fit['calib']

rat = rat.sel(wavelength=slice(750,755)).mean('wavelength')

ds_fix['calib'] = ds_fix['calib'] * rat

ds_fix = ds_fix.absem.calc_alpha()

ds_fix['alpha'].plot(row='motor')

plt.ylim(-0.1,1.1)

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2
from mhdpy.analysis.absem.fitting import pipe_fit_alpha_num_1

ds_absem_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fix)

ds_p_beer = ds_p.copy()


# %% [markdown]
# ## Tophat profile

# %%
from lmfit import Model

mod = Model(alpha_deq_solve, nK_profile=nK_profile)

params = mod.make_params(nK_max=1e-5, )

params.add('nK_m3', expr='nK_max*1e27')
params['nK_m3'].vary = False



from mhdpy.analysis.absem.fit_prep import pipe_fit_prep_alpha_2

da_fit = pipe_fit_prep_alpha_2(ds_fix)

ds_absem_fit, ds_p, ds_p_fit = da_fit.absem.perform_fit(mod, params)

ds_p_tophat = ds_p.copy()

# ds_alpha_fit = ds_alpha_fit.dropna('wavelength')

# ds_alpha_fit, ds_p, ds_p_stderr = ds_fit

# %%

ds_p['nK_m3'].plot()

plt.yscale('log')


# %%

ds_absem_fit.to_array('var').plot(hue='var', row='motor')

plt.ylim(0,1)

plt.xlim(763,775)

# %% [markdown]
# ## CFD Profiles

#%%


da_cfd_nK = ds_cfd['Yeq_K'].sel(kwt=1)

da_cfd_nK_norm = da_cfd_nK/da_cfd_nK.max('dist')

#TODO: motor coordinates not exactly the same, why? 
da_cfd_nK_norm = da_cfd_nK_norm.sel(motor=ds_fit.coords['motor'].values, method='nearest')
da_cfd_nK_norm = da_cfd_nK_norm.assign_coords(motor=ds_fit.coords['motor'].values) 
da_cfd_nK_norm


#%%


ds_absem_fit, ds_p, ds_p_fit = pipe_fit_alpha_num_1(ds_fix, perform_fit_kwargs={'nK_profile': da_cfd_nK_norm})

ds_p_cfd = ds_p.copy()

#%%
ds_absem_fit.to_array('var').plot(hue='var', row='motor')

plt.ylim(0,1)

plt.xlim(763,775)

# %%

# Centerline profile
fp_cfd = pjoin(os.getenv('REPO_DIR'), 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd_cl = xr.load_dataset(fp_cfd)
ds_cfd_cl = ds_cfd_cl.cfd.quantify_default()
ds_cfd_cl = ds_cfd_cl.cfd.convert_all_rho_number()

ds_cfd_cl = ds_cfd_cl.sel(kwt=1).sel(phi=0.8)

ds_cfd_cl['nK_m3'] = ds_cfd_cl['Yeq_K'].pint.to('particle/m^3')

# %%.


ds_p_cfd['nK_m3'].plot(marker='o', label='cfd profile')
ds_p_tophat['nK_m3'].plot(marker='o', label='tophat profile')
ds_p_beer['nK_m3'].plot(marker='o', label='beer lambert')

ds_cfd_cl['nK_m3'].sel(offset=0).plot(label='cfd centerline')
ds_cfd_cl['nK_m3'].sel(offset=1).plot(label='cfd 0.1')
ds_cfd_cl['nK_m3'].sel(offset=2).plot(label='cfd 0.3')

plt.yscale('log')

plt.legend()

# plt.ylim(1e21,1e22)
plt.ylim(1e17,1e22)


# %%
