
#%%[markdown]

# # 53x (seedramp) analysis

#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdpy.fileio.ct import load_df_cuttimes


#%%

fp_ct_seedramp = pjoin(REPO_DIR, 'experiment','metadata','ct_testcase_kwt.csv')
df_cuttimes_seedtcs = load_df_cuttimes(fp_ct_seedramp)
ds_lecroy = ppu.fileio.load_lecroy('53x', avg_mnum=True, df_ct_downselect=df_cuttimes_seedtcs)

# ds_lecroy = downselect_acq_time(ds_lecroy, df_cuttimes_seedtcs)

ds_fit = ds_lecroy
da_fit = ds_fit.mws.calc_mag_phase_AS()['AS']

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)
#%%

ds_p['decay'].mean('run').plot(marker='o')

# ds_p['decay'].plot()k


#%%
# Load cfd K+


ds_cfd = ppu.fileio.load_cfd_centerline(kwt_interp=ds_p.coords['kwt'].values)

ds_cfd = ds_cfd.sel(phi=0.8).sel(offset=0)

goldi_pos = Quantity(180,'mm')
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

ds_p

#%%
# Factor of 2 only occurs for K+

ds_species_decay = 1/(ds_p['decay'].pint.quantify('us')*ds_krb)

ds_species_decay = ds_species_decay.pint.to('particle/ml')

ds_species_decay

#%%

fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)


da_plot = ds_species_decay.to_array('species').mean('run')

da_plot

da_plot.plot(hue='species', marker='o', ax=axes[0])

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

ds_cfd['O2']

#%%

tau_exp = ds_p['decay'].mean('run').pint.quantify('us')
k_eff = 1/(tau_exp*ds_cfd['rho_number']*ds_cfd['O2'])

k_eff = k_eff.pint.to('cm^6/particle^2/s')
k_eff.attrs['long_name'] = 'Effective bimolecular rate'

k_eff.plot(marker='o')


#%%

fig, axes = plt.subplots(3, 1, figsize=(5,7), sharex=True, sharey=False)

# Plot decay
ds_p['decay'].mean('run', keep_attrs=True).plot(ax=axes[0], marker='o')
axes[0].set_title('MWS Time constant')
axes[0].set_ylabel('Decay constant (us)')   
axes[0].set_yscale('log')
axes[0].set_xlabel('')

# Plot ds_cfd['rho_number']
ds_cfd['rho_number'].pint.to('particle/ml').plot(ax=axes[1], label='all species', marker='o')
ds_cfd['O2'].pint.to('particle/ml').plot(ax=axes[1], label='O2', marker='o')
axes[1].set_title('CFD Particle Density')
axes[1].set_yscale('log')
axes[1].set_xlabel('')

axes[1].legend()

# Calculate and plot k_eff

k_eff.plot(marker='o', ax=axes[2])
axes[2].set_title('Effective Bimolecular Rate')
axes[2].set_yscale('log')


plt.xscale('log')

plt.xlabel('Nominal K wt%')
plt.savefig(pjoin(DIR_FIG_OUT, 'krb_eff_53x.png'))

plt.yscale('log')
# %%
