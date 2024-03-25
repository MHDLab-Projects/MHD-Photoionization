"""
Calcualtions of theoretical laser heating and ionization
"""


#%%

from mhdpy.analysis.standard_import import *
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')


tc = '536_power'
ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy['power']

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp_cfd_profiles = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_beam_Yeq.cdf')

ds_cfd = xr.load_dataset(fp_cfd_profiles)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd.coords['dist'] = ds_cfd.coords['dist'].pint.to('cm') # Fitting expects cm
ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)

da_cfd = ds_cfd['Yeq_KOH']

goldi_pos = Quantity(180, 'mm')

da_cfd = da_cfd.sel(motor=goldi_pos.magnitude, method='nearest').sel(kwt=1)

da_cfd.plot()

#%%

abscs = Quantity(1e-17, 'cm^2/particle')
kappa = (da_cfd*abscs).pint.to('1/cm')
kappa.plot()

#%%

f_A_nint = np.trapz(kappa, x=da_cfd.coords['dist'].pint.magnitude)

f_A_nint


#%%

from pi_paper_utils import LASER_POWER, LASER_AREA

absorption_length = Quantity(3, 'cm')
n_KOH = da_cfd.max('dist').item()

f_A = 1 - np.exp(-abscs*absorption_length*n_KOH)
f_A

#%%


heated_volume = LASER_AREA * absorption_length

heated_volume

#%%

powers = ds_lecroy['power'].values*LASER_POWER*f_A_nint

powers

#%%[markdown]

# Laser heating calculation

#%%

cantera_results_dir = pjoin(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')

ds_TP_params = xr.load_dataset(pjoin(cantera_results_dir, 'ds_TP_params.cdf'))

ds_sel = ds_TP_params[['Cp', 'rho']].sel(T=2000, P=1e5, phi=0.8, method='nearest')



Cp = Quantity(ds_sel['Cp'].mean('Kwt').item(), ds_sel['Cp'].attrs['units'])
rho = Quantity(ds_sel['rho'].mean('Kwt').item(), ds_sel['rho'].attrs['units'].replace("m3", 'm**3'))

Cp_vol = Cp * rho
Cp_vol


Cp_final = Cp_vol * heated_volume

Cp_final.to('J/K')

#%%

deltaTs = powers / Cp_final

deltaTs.to('K')

#%%[markdown]

# Max ionized potassium calculation

#%%

Ei = Quantity(4.34, 'eV/particle')

n_max = powers/Ei/heated_volume

n_max.to('particle/cm^3').magnitude


# %%
