
#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

from mhdpy.analysis import mws
from mhdpy.plot import dropna
from mhdpy.xr_utils import WeightedMeanAccessor
from mhdpy.analysis import absem
from mhdpy.coords.ct import downselect_acq_time

#%%


from mhdpy.fileio.ct import load_df_cuttimes

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)

#%%

tc = '53x'

ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
ds_lecroy = ds_lecroy.xr_utils.stack_run()


ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

ds_lecroy = ds_lecroy.drop(0,'kwt')

ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

ds_fit = ds_lecroy.mean('mnum')
da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))
#%%

ds_p['decay'].mean('run').plot(marker='o')

# ds_p['decay'].plot()k


#%%
# Load cfd K+


from mhdpy.pyvista_utils import CFDDatasetAccessor

fp = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)
ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)
ds_cfd = ds_cfd.interp(kwt=ds_lecroy.coords['kwt']).dropna('kwt', how='all')
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()
ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()


goldi_pos = Quantity(178,'mm')
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')


ds_cfd['Yeq_K+'].plot()

#%%

ds_cfd['T']

#%%

ds_cfd['rho_number'].pint.to('particle/ml').plot()

#%%

ds_cfd[['N2','O2','CO2','H2O']].to_array('var').plot(hue='var')
#%%
ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()

from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted

ds_krb = gen_ds_krb(ds_cfd['T'], ds_cfd['rho_number'])

ds_krb['O2_C'] = calc_krbO2_weighted(ds_cfd)


ds_krb.to_array('species').plot(hue='species', marker='o')

plt.yscale('log')

#%%[markdown]

# # Compare predicted and CFD species density

# predicted density from decay constant and hayhurst reaction rates

#%%

# Factor of 2 only occurs for K+

ds_species_decay = 1/(ds_p['decay'].pint.quantify('us')*ds_krb)

ds_species_decay = ds_species_decay.pint.to('particle/ml')

#%%

fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)


ds_species_decay.to_array('species').mean('run').plot(hue='species', marker='o', ax=axes[0])

axes[0].set_title('Predicted density $k_{r,species}$ and $\\tau_{mws}$')

ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O']].to_array('species').plot(hue='species', ax=axes[1], marker='o')

axes[1].set_title('CFD species density')

plt.yscale('log')

# %$[markdown]

# # calculate predicted decay constant from CFD species density

#%%

ds_species_cfd = ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O', 'Yeq_KOH', 'Yeq_K']]
ds_species_cfd = ds_species_cfd.rename({'Yeq_K+': 'K+', 'Yeq_OH': 'OH'})

from pi_paper_utils.kinetics import calc_krm

ds_krm = calc_krm(ds_krb, ds_species_cfd)

ds_tau = (1/ds_krm).pint.to('us')


#%%
for species in ds_tau.data_vars:
    ds_tau[species].plot(marker='o', label="CFD: {}".format(species))

ds_p['decay'].mean('run').plot(marker='o', label='Experiment')

plt.legend()
plt.yscale('log')
plt.ylabel('Decay constant (us)')
# %%

ds_tau

#%%

ds_tau['O2_C'].plot()

ds_p['decay'].mean('run').plot(marker='o', label='experiment')

plt.yscale('log')
# %%

tau_exp = ds_p['decay'].mean('run').pint.quantify('us')
k_eff = 1/(tau_exp*ds_cfd['rho_number']**2)

k_eff = k_eff.pint.to('cm^6/particle^2/s')

k_eff.plot(marker='o')

plt.yscale('log')
# %%
