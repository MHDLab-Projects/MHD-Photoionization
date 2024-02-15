# %% [markdown]
#  # Absem fitting with CFD profiles

# uses euler method to solve the differential equation for I(x)
# using normalized number density profile nK(x) and maximum number density nK_max

# %%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem

# %%

tc = '536_pos'

ds_absem = xr.load_dataset(pjoin(DIR_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem = ds_absem.xr_utils.stack_run()

ds_absem = ds_absem.absem.calc_alpha()
ds_absem = ds_absem.sel(wavelength=slice(750,790))

ds_absem

# %%

# make an artifical tophat profile. zero outside of the beam, constant inside the beam 

L = Quantity(1, 'cm')
L_center = Quantity(1, 'cm')

def tophat_profile(L, L_center, x):
    nK = np.zeros_like(x)
    nK[(x > L_center - L/2) & (x < L_center + L/2)] = 1
    return nK

x = np.linspace(0, 2, 100)

nK_profile = tophat_profile(L.magnitude, 1, x) 

plt.plot(x, nK_profile)


# %%[markdown]

# # I(x) solution functions
#
# Solve differential equation dI/dx = -kappa(x) I(x)
#
# We want to solve the differential equation for I(x) using the normalized number density profile nK_CFD(x).
#
# Given a max nK_max (Fit parameter) we calculate the absorption coefficient kappa(x) = Q(nK_max) * nK_profile(x).
#
# x is the distance along the beam axis.
#
# Below are three methods to solve the differential equation for I(x) using the kappa(x) profile.
#
# We first test everything with the tophat profile, and just selecting a single wavelength and nK_max

# %%

from scipy.integrate import solve_ivp
from mhdpy.analysis.absem.fitting import Q_2peak
from pint import Quantity



def dI_dx(x, I, kappa_profile):

    # differential equation for I(x)
    # dI/dx = -kappa(x) * I(x)

    kappa = np.interp(x, kappa_profile.index, kappa_profile)
    return -kappa * I

def calc_I_profile_deq(kappa_profile):
    """
    Solve the differential equation for I(x) using the kappa(x) profile
    """

    x = kappa_profile.index

    I_0 = 1
    sol = solve_ivp(dI_dx, [x[0], x[-1]], [I_0], args=(kappa_profile,), t_eval=x, method='BDF', max_step=0.1)

    I = sol.y[0]

    return I

from scipy.integrate import cumtrapz

def calc_I_profile_num(kappa_profile):
    """
    numerically integrate the differential equation for I(x) using the kappa(x) profile
    """
    x = kappa_profile.index

    # Compute -kappa(x) * I(x)
    dI_dx = -kappa_profile.values

    # Compute the cumulative integral
    I = cumtrapz(dI_dx, x, initial=0)

    # Add the initial condition
    I += 1

    return I

def calc_I_profile_euler(kappa_profile):
    """
    use the euler method to solve the differential equation for I(x) using the kappa(x) profile
    """


    x = kappa_profile.index
    I = np.zeros_like(kappa_profile.values)
    I[0] = 1  # initial condition

    # Euler method
    for i in range(1, len(x)):
        dI_dx = -kappa_profile.values[i-1] * I[i-1]
        deltax = x[i] - x[i-1]
        I[i] = I[i-1] + dI_dx * deltax

        if I[i] < 0:
            I[i] = 0

    return I


# %%


wl = Quantity(770, 'nm')

Q = Q_2peak(wl)
nK = Quantity(1e21, '1/m^3')

kappa_max = (Q*nK).to('1/cm').magnitude
kappa_profile = kappa_max * nK_profile
kappa_profile = pd.Series(kappa_profile, index=x)

I = calc_I_profile_euler(kappa_profile)

tau = I[-1]/I[0]

print(tau)

plt.plot(x, I)
plt.twinx()
plt.plot(x, kappa_profile, 'r')


# %%[markdown]

# Now pull in the CFD data to make a normalized number density profile nK_CFD(x) that is used to calculate the absorption coefficient kappa(x) = Q * nK(x).

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd_profiles = pjoin(REPO_DIR, 'modeling', 'cfd', 'output', 'line_profiles_beam_Yeq.cdf')

ds_cfd = xr.load_dataset(fp_cfd_profiles)
# ds_cfd = ds_cfd.coarsen(dist=5000, boundary='trim').mean()
ds_cfd.coords['dist'] = ds_cfd.coords['dist'] * 100
ds_cfd.coords['pos'] = ds_cfd.coords['pos'] * 10

ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

da_cfd = ds_cfd['Yeq_K']

da_cfd


# %%


da_cfd.plot(hue='pos', col='kwt')

plt.yscale('log')

plt.ylim(1e17,1e23)

# plt.gca().get_legend().remove()


# %%
da_cfd_sel = da_cfd.sel(kwt=1)

da_cfd_kappa = (da_cfd_sel/da_cfd_sel.max())*kappa_max

da_cfd_kappa 

Is = []

for pos in da_cfd_kappa.coords['pos'].values:
    da_sel = da_cfd_kappa.sel(pos=pos)
    kappa_profile = pd.Series(da_sel.values, index=da_sel.coords['dist'].values)

    I = calc_I_profile_deq(kappa_profile)

    Is.append(I)


da_I_pos = xr.DataArray(Is, coords=da_cfd_kappa.coords, dims=['pos', 'dist'])


# %%

da_I_pos.plot(hue='pos')

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

from scipy.ndimage import gaussian_filter1d

from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d

def alpha_deq_solve(x,
    nK_profile=None, #must be a kwarg...
    nK_max =1e-7,
    ):
    """
    Parameter
    ---------
    x : float array
        wavelength [nm]
    nK_profile : float array
        normalized number density profile [dimensionless]
    nK_max : float
        maximum number density [1/nm^3]

    Returns
    -------
    alpha : float array
        absorption []
    """

    if nK_profile is None:
        raise ValueError('nK_profile must be provided')

    Qs = Q_2peak(x)

    nK_max = Quantity(nK_max, '1/nm^3')

    # kappa0 = Quantity(0, '1/nm')

    alphas = np.zeros_like(x)

    for i, Q in enumerate(Qs):

        kappa_max = (Q*nK_max).to('1/cm').magnitude

        kappa_profile = kappa_max * nK_profile

        I = calc_I_profile_euler(kappa_profile)

        tau = Quantity(I[-1]/I[0], 'dimensionless')

        alpha = -tau + Quantity(1, 'dimensionless')
        alpha = alpha.magnitude

        alphas[i] = alpha


    # New method of spectrometer blurring, interpolate here then sample back. Allows for fitting fewer wavelengths
    # TODO: incorporate this into other fitting methods
    spect_res = 0.026
    blur_pixels_amount = 1 # Was a parameter... 
    wls_interp = np.arange(min(x), max(x), spect_res)

    # Interpolate to wls_interp
    f = interp1d(x, alphas, kind='linear', fill_value='extrapolate')
    alphas_interp = f(wls_interp)

    # Apply Gaussian blur
    blur_pixels_amount = 1
    alphas_blur = gaussian_filter1d(alphas_interp, blur_pixels_amount)

    # Interpolate back to x
    f = interp1d(wls_interp, alphas_blur, kind='linear', fill_value='extrapolate')
    alphas = f(x)


    return alphas

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

alpha = alpha_deq_solve(wls, nK_profile, nK_max=nK)
ds_test['alpha'].plot(marker='o')
plt.plot(wls, alpha)

plt.ylim(0,1)


# %%

from lmfit import Model

mod = Model(alpha_deq_solve, nK_profile=nK_profile)

params = mod.make_params(nK_max=1e-5, )

params.add('nK_m3', expr='nK_max*1e27')
params['nK_m3'].vary = False


# %%

from mhdpy.xr_utils import fit_da_lmfit

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

ds_fit = ds_fit.absem.remove_beta_offset(beta_offset_wls=slice(750,755))

ds_fit = ds_fit.absem.drop_alpha_peaks_negative()

# ds_fit = ds_fit.groupby('motor').apply(lambda x: x.absem.reduce_keep_wings())

alpha_fit = ds_fit['alpha']
#%%

alpha_fit.plot(row='motor')

plt.yscale('log')

# %% [markdown]
# ## Tophat profile

# %%


ds_absem_fit, ds_p, ds_p_fit = alpha_fit.absem.perform_fit(mod, params, method='iterative')

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

# %%

#TODO: pass a different profile in for each iteration

dss_p = []
dss_absem_fit = []

for pos in ds_fit.coords['motor'].values:


    nK_profile = da_cfd_kappa.sel(pos=pos, method='nearest').to_series()
    nK_max = nK_profile.max()
    nK_profile = nK_profile/nK_max

    mod = Model(alpha_deq_solve, nK_profile=nK_profile)

    params = mod.make_params(nK_max=1e-5, )
    params['nK_max'].min = 0

    params.add('nK_m3', expr='nK_max*1e27')
    params['nK_m3'].vary = False

    alpha_fit_sel = alpha_fit.sel(motor=pos)

    ds_absem_fit_sel, ds_p_sel, ds_p_fit_sel = alpha_fit_sel.absem.perform_fit(mod, params, method='iterative')

    dss_p.append(ds_p_sel)
    dss_absem_fit.append(ds_absem_fit_sel)

    # break

ds_p = xr.concat(dss_p, dim='motor')
ds_absem_fit = xr.concat(dss_absem_fit, dim='motor')

#%%
ds_absem_fit.to_array('var').plot(hue='var', row='motor')

plt.ylim(0,1)

plt.xlim(763,775)

# %%

# Centerline profile
fp_cfd = pjoin(os.getenv('REPO_DIR'), 'modeling', 'cfd', 'output', 'line_profiles_torchaxis_Yeq.cdf' )

ds_cfd_cl = xr.load_dataset(fp_cfd)

ds_cfd_cl = ds_cfd_cl.sel(kwt=1)

ds_cfd_cl = ds_cfd_cl.assign_coords(x = ds_cfd_cl.coords['x'].values - ds_cfd_cl.coords['x'].values[0])
ds_cfd_cl = ds_cfd_cl.assign_coords(x = ds_cfd_cl.coords['x'].values*1000)

ds_cfd_cl = ds_cfd_cl.cfd.quantify_default()
ds_cfd_cl = ds_cfd_cl.cfd.convert_all_rho_number()


ds_cfd_cl['nK_m3'] = ds_cfd_cl['Yeq_K'].pint.to('1/m^3')




# %%


ds_p['nK_m3'].plot(marker='o', label='cfd profile')
ds_p_tophat['nK_m3'].plot(marker='o', label='tophat profile')

ds_cfd_cl['nK_m3'].plot(label='cfd centerline')


plt.yscale('log')

plt.legend()

plt.ylim(1e18,1e23)

# %%



