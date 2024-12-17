#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils #Sets matplotlib style
from pi_paper_utils.kinetics import gen_ds_krb, calc_krbO2_weighted, calc_krm

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_p_stats = xr.open_dataset(pjoin(data_directory, '53x_ds_p_stats.cdf')).pint.quantify()
ds_species_cfd = xr.open_dataset(pjoin(data_directory, '53x_ds_species_cfd.cdf')).pint.quantify()
ds_tau = xr.open_dataset(pjoin(data_directory, 'ds_tau.cdf')).pint.quantify()
# %%
tau_exp = ds_p_stats['mws_fit_decay_exp_mean'].pint.quantify('us')

k_eff_bm = 1/(tau_exp*ds_species_cfd['O2'])

k_eff_bm = k_eff_bm.pint.to('ml/particle/s')
k_eff_bm.attrs['long_name'] = '$k_{r,bm,eff}$'


ds_krb = gen_ds_krb(ds_species_cfd['T'], ds_species_cfd['rho_number'])

ds_krb['O2_S'] = calc_krbO2_weighted(ds_species_cfd)


#%%
# plot on only bimolecular rate

fig, ax = plt.subplots(1, 1, figsize=(5,2))

ln = k_eff_bm.plot(marker='x', ax=ax, label='Experiment')
g = ds_krb[['O2_A', 'O2_G', 'O2_S']].to_array('rxn').pint.dequantify().plot(hue='rxn', marker='o')

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
# %%
