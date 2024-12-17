#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from pi_paper_utils.constants import *

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

# ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
# ds_absem = xr.open_dataset(pjoin(data_directory, 'ds_pos_alpha_fit.cdf')).xr_utils.stack_run()


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

ds_cfd[['Kp','Yeq_K+']].to_array('var').plot(hue='var', col='phi', row='kwt')

plt.yscale('log')

#%%

from pi_paper_utils.kinetics import calc_krm

# ds_krm = calc_krm(ds_krb, ds_cfd.rename({"Yeq_K+": "K+"}))
ds_krm = calc_krm(ds_krb, ds_cfd.rename({"Kp": "K+"}))

ds_tau = (1/ds_krm).pint.to('us')

#%%

da_plot = ds_tau[['K+','O2_A','O2_G', 'OH', 'H2O']].to_array('var').sel(kwt=1, method='nearest')

g= da_plot.plot(hue='var', col='phi')

plt.yscale('log')

plt.ylim(0.01,1e4)

g.axes[0,0].set_ylabel("Tau [$\mu s$]")

for ax in g.axes.flatten():
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')


plt.savefig(pjoin(DIR_FIG_OUT, 'krm_cfd_pos.png'))

#%%

g = ds_tau[['K+','O2_A','O2_G', 'OH', 'H2O']].to_array('var').plot(hue='var', col='phi', row='kwt')

plt.yscale('log')

plt.ylim(0.1,1e4)


for ax in g.axes[:,0]:
    ax.set_ylabel("Tau [$\mu s$]")
    ax.axvline(goldi_pos.magnitude, color='k', linestyle='--')

#%%

ds_cfd_gold = ds_cfd.sel(x = goldi_pos, method='nearest')

#%%
