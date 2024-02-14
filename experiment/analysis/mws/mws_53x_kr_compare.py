
#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

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

ds_lecroy = xr.load_dataset(pjoin(DIR_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))
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

fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

ds_cfd = xr.load_dataset(fp)

ds_cfd = ds_cfd.interp(kwt=ds_lecroy.coords['kwt']).dropna('kwt', how='all')

ds_cfd = ds_cfd.cfd.convert_species_rho()

goldi_pos = ds_cfd['x'].min().item() + 0.18
ds_cfd = ds_cfd.sel(x = goldi_pos, method='nearest')


ds_cfd['Yeq_K+'].plot()

#%%

ds_cfd['T']

#%%

ds_cfd['rho'].pint.to('1/cm^3').plot()

#%%

ds_cfd[['N2','O2','CO2','H2O']].to_array('var').plot(hue='var')
#%%

# combining k4 from 
# 1. Axford, S.D.T., and Hayhurst, A.N. (1997). Mass spectrometric sampling of negative ions from flames of hydrogen and oxygen: the kinetics of electron attachment and detachment in hot mixtures of H2O, O2, OH and HO2. Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences 452, 1007–1033. 10.1098/rspa.1996.0051.


k4_species = {
    'N2' : Quantity(1e-32, 'cm^6/s'),
    'H2O': Quantity(8e-30, 'cm^6/s'),
    'CO2': Quantity(1e-31, 'cm^6/s'), 
    'O2' : Quantity(2e-30, 'cm^6/s'),
}

weighted_avg = 0

for species, k4 in k4_species.items():
    weighted_avg += k4*ds_cfd[species]

weighted_avg

#%%

rxn_rates = {
    'K+':  Quantity(4e-24, 'K*cm^6/s')*(1/ds_cfd['T'])*ds_cfd['rho'],
    'OH': Quantity(3e-31, 'cm^6/s')*ds_cfd['rho'],
    # 'O2': Quantity(5e-31, 'cm^6/s')*ds_cfd['rho'],
    'O2': weighted_avg,
    # 'O2': Quantity(6e-34, 'cm^6/s')*ds_cfd['rho'],
    'H2O': Quantity(1.6e-6, 'cm^3/s')*np.exp(-(Quantity(36060, 'K')/ds_cfd['T'])),
}

ds_kr = xr.Dataset(rxn_rates).pint.to('cm^3/s')

ds_kr.to_array('species').plot(hue='species', marker='o')

plt.yscale('log')

#%%[markdown]

# # Compare predicted and CFD species density

# predicted density from decay constant and hayhurst reaction rates

#%%

# Factor of 2 only occurs for K+

ds_species_decay = 1/(ds_p['decay'].pint.quantify('us')*ds_kr)

ds_species_decay = ds_species_decay.pint.to('1/cm^3')

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

ds_species_cfd = ds_cfd[['Yeq_K+', 'Yeq_OH', 'O2', 'H2O']]

ds_species_cfd = ds_species_cfd.pint.to('1/cm^3')

ds_species_cfd = ds_species_cfd.rename({'Yeq_K+': 'K+', 'Yeq_OH': 'OH'})

ds_species_cfd

#%%

ds_tau = 1/(ds_species_cfd*ds_kr)

ds_tau = ds_tau.pint.to('us')

ds_tau

#%%
for species in ds_tau.data_vars:
    ds_tau[species].plot(marker='o', label="CFD: {}".format(species))

ds_p['decay'].mean('run').plot(marker='o', label='Experiment')

plt.legend()
plt.yscale('log')
plt.ylabel('Decay constant (us)')
# %%

ds_tau['O2'].plot()

ds_p['decay'].mean('run').plot(marker='o', label='experiment')

plt.yscale('log')
# %%

tau_exp = ds_p['decay'].mean('run').pint.quantify('us')
k_eff = 1/(tau_exp*ds_cfd['rho']**2)

k_eff = k_eff.pint.to('cm^6/s')

k_eff.plot(marker='o')

plt.yscale('log')
# %%
