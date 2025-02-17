
#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdlab.analysis.absem.fitting import pipe_fit_alpha_num_1

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

# Position dependent AES and MWT
ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_alpha_fit = xr.open_dataset(pjoin(data_directory, 'ds_pos_alpha_fit.cdf')).xr_utils.stack_run()
# Non position dependent barrel AES
ds_p_barrel = xr.open_dataset(pjoin(data_directory, 'ds_p_barrel.cdf')).xr_utils.stack_run()

ds_cfd_cl = ppu.fileio.load_cfd_centerline()
ds_cfd_cl = ds_cfd_cl.sel(kwt=1).sel(phi=0.8)

ds_cfd_beam = ppu.fileio.load_cfd_beam()
ds_cfd_beam = ds_cfd_beam.sel(offset=0).sel(kwt=1).sel(phi=0.8)


#%%

plt.figure()

fig, ax = plt.subplots(4, figsize=(4,8))  

plt.sca(ax[0])

ds_plot = ds_cfd_beam.sel(mp='barrel').isel(motor=0).drop('motor')

ds_plot[ppu.CFD_K_SPECIES_NAME].plot(label='K')
ds_plot[ppu.CFD_KOH_SPECIES_NAME].plot(label='KOH')
ds_plot[ppu.CFD_allK_SPECIES_NAME].plot(label='All K')

# plt.yscale('log')
plt.ylabel('Concentration\n(particle/ml)')
plt.legend()

plt.sca(ax[1])

ds_plot['T'].plot(label='T')

plt.ylabel('Temperature (K)')

plt.sca(ax[2])

ds_plot['p'].pint.to('atm').plot(label='p')

plt.ylabel('Pressure (atm)')

plt.sca(ax[3])


rho_num = ds_plot.cfd.calc_rho_number()

rho_num.plot(label='rho')

plt.ylabel('Density (particle/ml)')

# plt.yscale('log')

for a in ax.flatten():
    a.set_xlim(4,6)
    a.set_title('')

ax[-1].set_xlabel('Distance (cm)')
plt.savefig(pjoin('output','figures','cfd_barrel_beam_K_KOH_compare.png'))   

#%%



#%%

plt.figure()

da_plot = ds_cfd_cl['p'].pint.to('atm').sel(offset=0)

da_plot.plot()

plt.xlim(0,200)

plt.savefig(pjoin('output','figures','cfd_centerline_p.png'))

#%%

min_p = ds_cfd_cl['p'].min().item()
max_p = ds_cfd_cl['p'].max().item()

min_p, max_p


#%%

stack_dims = ['date','run_num','motor','phi']

tc = '536_pos'
ds_absem_536_pos = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))
ds_absem_536_pos = ds_absem_536_pos.expand_dims('phi').stack(temp=stack_dims)

ds_absem = xr.concat([ds_absem_536_pos], dim='temp').unstack('temp')
ds_absem = ds_absem.xr_utils.stack_run()
ds_absem = ds_absem.absem.calc_alpha().sel(wavelength=slice(750,790))

#%%


#TODO: motor coordinates not exactly the same, why?
ds_cfd_beam = ds_cfd_beam.sel(motor=ds_absem.coords['motor'].values, method='nearest')
ds_cfd_beam = ds_cfd_beam.assign_coords(motor=ds_absem.coords['motor'].values)
da_cfd_beam = ds_cfd_beam[ppu.CFD_K_SPECIES_NAME] / ds_cfd_beam[ppu.CFD_K_SPECIES_NAME].max('dist')


#%%
# fit with 2 atm pressure

# pipe_fit_alpha_num_1.params['p'].value = 2
ds_fit = ds_absem.mean('mnum').sel(mp='barrel').dropna('run', how='all')
ds_fit = ds_fit.mean('motor', keep_attrs=True)

nK_profile = da_cfd_beam.sel(mp='barrel').isel(motor=0)
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds_p_base = ds_p.copy()

#%%


pipe_fit_alpha_num_1.params['p'].value = max_p.to('atm').magnitude
ds_fit = ds_absem.mean('mnum').sel(mp='barrel').dropna('run', how='all')
ds_fit = ds_fit.mean('motor', keep_attrs=True)

nK_profile = da_cfd_beam.sel(mp='barrel').isel(motor=0)
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds_p_max_p = ds_p.copy()

#%%

pipe_fit_alpha_num_1.params['p'].value = min_p.to('atm').magnitude
ds_fit = ds_absem.mean('mnum').sel(mp='barrel').dropna('run', how='all')
ds_fit = ds_fit.mean('motor', keep_attrs=True)

nK_profile = da_cfd_beam.sel(mp='barrel').isel(motor=0)
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds_p_min_p = ds_p.copy()

#%%

# Tophat profile with 1 atm pressure

pipe_fit_alpha_num_1.params['p'].value = 1
# tophat 

L = Quantity(1, 'cm')
L_center = Quantity(5, 'cm')

def tophat_profile(L, L_center, x):
    nK = np.zeros_like(x)
    nK[(x > L_center - L/2) & (x < L_center + L/2)] = 1
    return nK

x = da_cfd_beam.coords['dist'].values
nK_profile = tophat_profile(L.magnitude, L_center.magnitude, x) 

# plt.plot(x, nK_profile)

# make a dataset with same strucutre as da_cfd_beam, but with nK_profile as data along the 'dist' dimension

nK_profile_tophat = xr.DataArray(nK_profile, coords={'dist':x}, dims=['dist'])

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':nK_profile_tophat})

ds_alpha_fit['alpha_full'] = ds_fit['alpha']
ds = ds_alpha_fit[['alpha_red', 'alpha_fit', 'alpha_full']].rename({'alpha_red': 'fitted', 'alpha_fit':'fit', 'alpha_full':'full'})
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

ds_p_tophat = ds_p.copy()


# %%

plt.figure()


ds_cfd_cl['nK_m3'] = ds_cfd_cl[ppu.CFD_K_SPECIES_NAME].pint.to('particle/m^3')
ds_cfd_cl['nK_m3'] = ds_cfd_cl['nK_m3'].sel(offset=0).drop('offset')

ds_cfd_cl['nK_m3'].plot(color='black', label ='CFD centerline', linestyle='-.')


nK_barrel_mean = ds_p_base['nK_m3'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_base['nK_m3'].std('run').sel(phi=0.8, method='nearest')
plt.gca().errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='o', label='1 Atm', capsize=5, )

max_p_str = '{:.2f} atm'.format(max_p.to('atm').magnitude)
min_p_str = '{:.2f} atm'.format(min_p.to('atm').magnitude)

nK_barrel_mean = ds_p_max_p['nK_m3'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_max_p['nK_m3'].std('run').sel(phi=0.8, method='nearest')
plt.gca().errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='x', label=f'Max ({max_p_str})', capsize=5, )

nK_barrel_mean = ds_p_min_p['nK_m3'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_min_p['nK_m3'].std('run').sel(phi=0.8, method='nearest')
plt.gca().errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='blue', marker='x', label=f'Min ({min_p_str})', capsize=5, )

nK_barrel_mean = ds_p_tophat['nK_m3'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_tophat['nK_m3'].std('run').sel(phi=0.8, method='nearest')
plt.gca().errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='x', label='Tophat (1 atm)', capsize=5, )


plt.legend()
plt.yscale('log')

plt.ylim(5e21, 2e23)

plt.xlim(0,10)

plt.savefig(pjoin('output','figures','nK_barrel_cases.png'))
# %%
