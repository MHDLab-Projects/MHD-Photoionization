#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')


tc = '536_power'
ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

ds_lecroy['power']

#%%

from mhdpy.pyvista_utils import CFDDatasetAccessor

fp = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)

ds_cfd = ds_cfd.sel(offset=0)
# ds_cfd = ds_cfd.sel(kwt=1)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd

goldi_pos =  Quantity(178, 'mm')


#%%

laser_power = Quantity(10, 'mJ')
laser_area = Quantity(0.4, 'cm^2')


#TODO: this needs to be calculated from CFD KOH
# f_A = Quantity(0.05, 'dimensionless')
abscs = Quantity(1e-17, 'cm^2/particle')
absorption_length = Quantity(2, 'cm')
n_KOH = ds_cfd.sel(x=goldi_pos.magnitude, method='nearest').sel(kwt=1, phi=0.8)['Yeq_KOH'].item()

f_A = 1 - np.exp(-abscs*absorption_length*n_KOH)
f_A

#%%


heated_volume = laser_area * absorption_length

heated_volume

#%%

powers = ds_lecroy['power'].values*laser_power*f_A

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
