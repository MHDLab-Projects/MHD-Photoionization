#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from pi_paper_utils.constants import *

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')


ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(offset=0)
ds_cfd

goldi_pos = Quantity(180, 'mm')

#%%
# Calculate krb

ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()

from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted

ds_krb = gen_ds_krb(ds_cfd['T'], ds_cfd['rho_number'])

ds_krb['O2_S'] = calc_krbO2_weighted(ds_cfd)

ds_krb

#%%

g = ds_krb[['K+', 'OH', 'O2_A', 'H2O']].to_array('var').plot(hue='var', col='phi', row='kwt')

plt.yscale('log')

plt.ylim(1e-14,)

for ax in g.axes.flatten():
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')

#%%

# show species densities on similar plot 

g = ds_cfd[[ppu.CFD_KOH_SPECIES_NAME, 'OH', 'O2', 'H2O']].to_array('var').plot(hue='var', col='phi', row='kwt')

plt.yscale('log')

plt.ylim(1e12,)


for ax in g.axes.flatten():
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')

#%%
plt.figure()
# show species densities on similar plot 
da_plot = ds_cfd[[CFD_Kp_SPECIES_NAME, 'OH', 'O2', 'H2O']].to_array('var').sel(kwt=1, method='nearest')


g = da_plot.plot(hue='var', col='phi')

plt.yscale('log')

plt.ylim(1e12,)


for ax in g.axes.flatten():
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')


plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_species_pos.png'))

#%%


#%%

from pi_paper_utils.kinetics import calc_krm

#TODO: hack until consistent naming is implemented
ds_krm_Yeq = calc_krm(ds_krb, ds_cfd.rename({"Yeq_K+": "K+"})).rename({"K+": "Yeq_K+"})
ds_krm_base = calc_krm(ds_krb, ds_cfd.rename({"Kp": "K+"}))

ds_krm = ds_krm_Yeq.combine_first(ds_krm_base)

ds_tau = (1/ds_krm).pint.to('us')

ds_tau

#%%


tau_species_list = ['K+', 'Yeq_K+', 'OH', 'H2O', 'O2_A', 'O2_G']
cfd_species_list = ['Kp', 'Yeq_K+', 'OH', 'H2O', 'O2']

ds_cfd[cfd_species_list].sel(kwt=1,method='nearest').to_array('var').plot(hue='var', col='phi')

plt.yscale('log')

plt.ylim(1e9,)

#%%

da_plot = ds_tau[tau_species_list].to_array('var').sel(kwt=1, method='nearest')

g= da_plot.plot(hue='var', col='phi')

plt.yscale('log')

plt.ylim(0.01,1e4)

g.axes[0,0].set_ylabel("Tau [$\mu s$]")

for ax in g.axes.flatten():
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')



#%%

da_plot = ds_tau[['K+', 'O2_A']].to_array('var').sel(kwt=1, method='nearest')

da_plot = da_plot.sel(x=slice(100,200))

g= da_plot.plot(row='var', hue='phi', figsize=(3,5), sharey=False)

for ax in g.axes.flatten():
    # ax.set_yscale('log')
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')
    ax.set_ylabel("Tau [$\mu s$]")

g.axes[0,0].set_ylim(0,40)
g.axes[1,0].set_ylim(0,1)

plt.savefig(pjoin(DIR_FIG_OUT, 'tau_cfd_pos_phi_compare.png'))

#%%
fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

# First plot
ds_cfd[cfd_species_list].sel(kwt=1, method='nearest').sel(phi=0.8, method='nearest').to_array('var').plot(ax=axes[0], hue='var')
axes[0].set_yscale('log')
axes[0].get_legend().set_bbox_to_anchor([0, 0, 1.3, 1])
axes[0].set_ylabel('Species Concentration [#/ml]')
axes[0].axvline(goldi_pos.magnitude, color='k', linestyle='--')

# Second plot
ds_tau[tau_species_list].to_array('var').sel(kwt=1, method='nearest').sel(phi=0.8, method='nearest').plot(ax=axes[1], hue='var')
axes[1].set_yscale('log')
axes[1].set_ylim(0.01, 1e4)
axes[1].set_ylabel('Tau [$\mu s$]')
# axes[1].axhline(tau_observed.magnitude, color='gray', linestyle='--')
axes[1].axvline(goldi_pos.magnitude, color='k', linestyle='--')
axes[1].get_legend().set_bbox_to_anchor([0, 0, 1.3, 1])

plt.tight_layout()

plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_species_tau_pos.png'))


#%%

fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

ds_plot = ds_cfd.sel(kwt=1, method='nearest').pint.dequantify()
ds_plot = ds_plot.sel(phi=0.8, method='nearest')

da_plot = ds_plot[['Kp', 'em','OHm', 'O2m']].to_array('var')

g = da_plot.plot(hue='var', ax=axes[0])

axes[0].set_yscale('log')

da_plot = ds_plot['sigma_e']

g = da_plot.plot(hue='var', ax=axes[1])

axes[1].set_yscale('log')
axes[1].set_ylim(1e-15,1e2)


plt.savefig(pjoin(DIR_FIG_OUT, 'cfd_species_charged_sigma_pos.png'))

#%%

ds_plot

#%%

g = ds_tau[['K+','Yeq_K+','O2_A','O2_G', 'OH', 'H2O']].to_array('var').plot(hue='var', col='phi', row='kwt')

plt.yscale('log')

plt.ylim(0.1,1e4)


for ax in g.axes[:,0]:
    ax.set_ylabel("Tau [$\mu s$]")
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')

#%%

ds_cfd_gold = ds_cfd.sel(x = goldi_pos, method='nearest')

#%%
