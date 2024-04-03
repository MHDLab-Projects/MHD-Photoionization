
#%%[markdown]

# # 53x (seedramp) O2 kinetics comparison

# TODO: only the final plot is being used, clean up and move to analysis


#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils #Sets matplotlib style
from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted, calc_krm

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()

#%%

krm_eff = ds_p_stats['mws_fit_krm_mean']
krm_eff.plot()

mean_krm = krm_eff.pint.to('1/s').mean().item()

plt.text(0.05, 0.85, 'Average: {:.3e}'.format(mean_krm), transform=plt.gca().transAxes)

print(mean_krm)

#%%

ds_krb = gen_ds_krb(ds_species_cfd['T'], ds_species_cfd['rho_number'])

ds_krb['O2_S'] = calc_krbO2_weighted(ds_species_cfd)

ds_krb.to_array('species').plot(hue='species', marker='o')

plt.yscale('log')

#%%[markdown]

# # Compare predicted and CFD species density

# predicted density from exponential time constant and hayhurst reaction rates


#%%
# Factor of 2 only occurs for K+

ds_species_decay = 1/(ds_p_stats['mws_fit_decay_exp_mean'].pint.quantify('us')*ds_krb)

ds_species_decay = ds_species_decay.pint.to('particle/ml')

ds_species_decay

#%%

fig, axes = plt.subplots(1, 2, figsize=(10,5), sharex=True, sharey=True)


da_plot = ds_species_decay.to_array('species')

da_plot

da_plot.plot(hue='species', marker='o', ax=axes[0])

axes[0].set_title('Predicted density $k_{r,species}$ and $\\tau_{mws}$')

ds_species_cfd[['K+', 'OH', 'O2', 'H2O']].to_array('species').plot(hue='species', ax=axes[1], marker='o')

axes[1].set_title('CFD species density')

plt.yscale('log')

# %$[markdown]

# # calculate predicted exponential time constant from CFD species density

#%%
ds_krm = calc_krm(ds_krb, ds_species_cfd)

ds_tau = (1/ds_krm).pint.to('us')

#%%
for species in ds_tau.data_vars:
    ds_tau[species].plot(marker='o', label="CFD: {}".format(species))

ds_p_stats['mws_fit_decay_exp_mean'].plot(marker='o', color='black', label='Experiment')

plt.legend()
plt.yscale('log')
plt.ylabel('Decay constant (us)')
# %%

ds_tau

#%%

ds_tau['O2_S'].plot()

ds_p_stats['mws_fit_decay_exp_mean'].plot(marker='o', label='experiment')

plt.yscale('log')
# %%

ds_species_cfd['O2']

#%%

tau_exp = ds_p_stats['mws_fit_decay_exp_mean'].pint.quantify('us')

k_eff_bm = 1/(tau_exp*ds_species_cfd['O2'])

k_eff_bm = k_eff_bm.pint.to('ml/particle/s')
k_eff_bm.attrs['long_name'] = '$k_{r,bm,eff}$'

k_eff_bm.plot(marker='o')


k_eff_tm = 1/(tau_exp*ds_species_cfd['rho_number']*ds_species_cfd['O2'])

k_eff_tm = k_eff_tm.pint.to('ml^2/particle^2/s')
k_eff_tm.attrs['long_name'] = '$k_{r,tm,eff}$'

k_eff_tm.plot(marker='o')


#%%

fig, axes = plt.subplots(4, 1, figsize=(5,8), sharex=True, sharey=False)

# Plot decay
ds_p_stats['mws_fit_decay_exp_mean'].plot(ax=axes[0], marker='o')
axes[0].set_title('MWS Time constant')
axes[0].set_ylabel('Decay constant (us)')   
axes[0].set_yscale('log')
axes[0].set_xlabel('')

# Plot ds_species_cfd['rho_number']
ds_species_cfd['rho_number'].pint.to('particle/ml').plot(ax=axes[1], label='all species', marker='o')
ds_species_cfd['O2'].pint.to('particle/ml').plot(ax=axes[1], label='O2', marker='o')
axes[1].set_title('CFD Particle Density')
axes[1].set_yscale('log')
axes[1].set_xlabel('')

axes[1].legend()

# Calculate and plot k_eff_bm
k_eff_bm.plot(marker='o', ax=axes[2])
axes[2].set_title('Effective Bimolecular Rate')
axes[2].set_yscale('log')

k_eff_avg = k_eff_bm.mean('kwt').item()
axes[2].text(0.0, 0.85, 'Average: {:.3e}'.format(k_eff_avg), transform=axes[2].transAxes)
axes[2].set_xlabel('')
axes[2].set_ylabel('Effective Rate\n(ml/particle/s)')

# Calculate and plot k_eff_tm
k_eff_tm.plot(marker='o', ax=axes[3])
axes[3].set_title('Effective Termolecular Rate')
axes[3].set_yscale('log')
axes[3].set_ylabel('Effective Rate\n(ml^2/particle^2/s)')

k_eff_avg = k_eff_tm.mean('kwt').item()
axes[3].text(0.0, 0.85, 'Average: {:.3e}'.format(k_eff_avg), transform=axes[3].transAxes)


plt.xscale('log')

plt.xlabel('Nominal K wt%')
# plt.savefig(pjoin(DIR_FIG_OUT, 'krb_eff_53x.png'))

plt.yscale('log')
# %%

# plot on only bimolecular rate

fig, ax = plt.subplots(1, 1, figsize=(5,2))

k_eff_bm.plot(marker='o', ax=ax)

plt.ylabel('$k_{r,b,eff}$\n(ml/particle/s)')

plt.xscale('log')
plt.xlim(0.05, )
plt.title('')

kr_eff_avg = k_eff_bm.mean('kwt').item()

ax.text(0.0, 0.85, 'Average: {:.3e}'.format(kr_eff_avg), transform=ax.transAxes)


# plt.savefig(pjoin(DIR_FIG_OUT, 'krb_eff_bm_53x.png'))


#%%


# plot on only bimolecular rate

fig, ax = plt.subplots(1, 1, figsize=(5,2))

ln = k_eff_bm.plot(marker='x', ax=ax, label='Experiment')
g = ds_krb[['O2_A', 'O2_G', 'O2_S']].to_array('rxn').plot(hue='rxn', marker='o')

lns = ln + g
labs = ['Experiment', 'O2_A', 'O2_G', 'O2_S']

leg = plt.legend(lns, labs)
leg.set_bbox_to_anchor((1,0.8))

plt.ylabel('$k_{r,b}$\n(ml/particle/s)')
plt.xlabel('Nominal K wt%')

plt.xscale('log')
plt.yscale('log')
plt.xlim(0.05, )
plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'krb_eff_bm_53x_O2.png'))

#%%

labs