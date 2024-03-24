# %%

from mhdpy.analysis.standard_import import *

from pi_paper_utils import abscs, noneq
import xyzpy
from dotenv import load_dotenv
load_dotenv()

REPO_DIR = os.getenv('REPO_DIR')


cantera_data_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset', 'output')
PI_modeling_dataset_dir = os.path.join(REPO_DIR, 'modeling', 'viability', 'dataset','output')
if not os.path.exists('output'): os.mkdir('output')

# %%
ds_TP_species = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})
ds_TP_params = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_params.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_TP_species_rho = xr.open_dataset(os.path.join(cantera_data_dir, 'ds_TP_species_rho.cdf'))#.sel({'phi': 0.7, 'Kwt': 0.001})

ds_NE = xr.open_dataset(pjoin(REPO_DIR, 'modeling', 'viability', 'dataset', 'output', 'ds_NE.cdf')).squeeze()
# %%
sig_bk = ds_TP_params['sigma'].sel(T=3000, phi=0.8, Kwt=0.01, P=1e5, method='nearest')
# sig_bk = ds_NE['sigma'].sel(T=3000, method='nearest') #This was used previously, incorrectly
sigma_bl = ds_NE['sigma'].sel(phi=0.8, eta='perf', rxn='mm_sum', Kwt=0.01, P=1e5).squeeze()

combos_l_b = {'l_b' :  np.linspace(0,1,100)}

r = xyzpy.Runner(noneq.calc_dsigma_tot, constants = {'sigma_bk' :sig_bk, 'sigma_bl': sigma_bl }, var_names = ['sigma_tot'])
da_dsigma_tot = r.run_combos(combos_l_b).squeeze()
da_dsigma_tot.name = 'enhancement factor'
da_dsigma_tot.attrs = {}


# gamma_bl = ds_NE['gamma']*da_dsigma_tot

# gamma_bl.name = 'gamma_bl'
# gamma_bl.to_netcdf('output/gamma_bl.cdf')

#%%

sigma_bl.plot()
# sig_bk
# plt.yscale('log')

#%%

# seldict = dict(phi=0.8, eta='perf', rxn='mm_sum', Kwt=0.01, P=1e5)

da_plot = da_dsigma_tot.sel(T=[1500,2000,2500,2700,3000], method='nearest')

da_plot.plot(hue='T')
plt.gca().get_legend().set_bbox_to_anchor((1, 1))

plt.yscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, 'bl_enhance.png'))
# %%


da_dsigma_tot
# %%
