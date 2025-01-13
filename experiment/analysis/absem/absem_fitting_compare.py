#%%[markdown]

# # Absorption emission fitting algorithm comparison


#%%

from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdlab.analysis.absem.fit_prep import interp_alpha
from mhdlab.analysis.absem.fitting import gen_model_alpha_blurred 
from mhdlab.analysis import absem

from mhdlab.analysis.absem.fitting import alpha_2peak
from mhdlab.xr_utils import fit_da_lmfit_global
from mhdlab.xr_utils import fit_da_lmfit

dss_p = []
dss_p_stderr = []

model, pars = gen_model_alpha_blurred(assert_xs_equal_spacing=False)

# %%

ds_absem = ppu.fileio.load_absem('53x')

#%%

ds_sel = ds_absem.sel(mp='barrel').sel(run=('2023-05-24', 1))
ds_sel = ds_sel.dropna('wavelength', how='all')
ds_sel = ds_sel.sel(wavelength=slice(750,790))

#%%

ds_sel['alpha'].plot(hue='mnum', row='kwt')

plt.ylim(-0.1,1.1)

# add horizontal line at 0 
plt.axhline(0, color='k', linestyle='--', linewidth=0.5)

#%%[markdown]

# ## perfrom average before alpha calculation

#%%

ds_sel_2 = ds_sel.mean('mnum').absem.calc_alpha()

ds_absem_fit, ds_p, ds_p_stderr = ds_sel_2['alpha'].absem.perform_fit(model, pars, method='iterative')

dss_p.append(ds_p.assign_coords(method='preavg'))
dss_p_stderr.append(ds_p.assign_coords(method='preavg'))

#%%

da = ds_absem_fit.to_array('var')

da.plot(hue='var', row='kwt')

#%%[markdown]

# ## Global Least squares reduction. 

# All measurement numbers are concatenated into one array, and fitted in one optimization process

#%%

# da_alpha = interp_alpha(da_sel)
da_sel = ds_sel['alpha'].dropna('mnum',how='all')

ds_absem_fit, ds_p, ds_p_stderr = da_sel.absem.perform_fit(model, pars, method='global')

dss_p.append(ds_p.assign_coords(method='global'))
dss_p_stderr.append(ds_p_stderr.assign_coords(method='global'))


#%%

da = ds_absem_fit.mean('mnum').to_array('var')

da.plot(hue='var', row='kwt')

#%%[markdown]

# ## Fitting each individual spectrum

#%%
da_sel = ds_sel['alpha'].dropna('mnum',how = 'all')


ds_absem_fit, ds_p, ds_p_stderr = da_sel.absem.perform_fit(model, pars, method='iterative')

dss_p.append(ds_p.mean('mnum').assign_coords(method='individual'))
dss_p_stderr.append(ds_p.std('mnum').assign_coords(method='individual'))
#%%

da = ds_absem_fit.to_array('var').sel(kwt=1, method='nearest')
ds_p_sel = ds_p.sel(kwt=1, method='nearest')
ds_p_stderr_sel = ds_p_stderr.sel(kwt=1, method='nearest')
da = da.dropna('mnum','all')

#%%

g = da.plot(hue='var', row='mnum')

for i , ax in enumerate(g.axes.flatten()):
    # ax.set_ylim(-0.1,1.1)
    # ax.set_xlim(750,790)
    mnum = da.coords['mnum'].values[i]
    nK_m3 = ds_p_sel.sel(mnum=mnum)['nK_m3'].item()
    nK_m3_stderr = ds_p_stderr_sel.sel(mnum=mnum)['nK_m3'].item()

    ax.set_title('mnum = {}, nK_m3 = {:.2e} +/- {:.2e}'.format(mnum, nK_m3, nK_m3_stderr))

#%%

plt.errorbar(ds_p_sel.coords['mnum'] , ds_p_sel['nK_m3'], ds_p_stderr_sel['nK_m3'], fmt='o')

#%%
nK_mean = ds_p_sel['nK_m3'].mean('mnum').item()
nK_std = ds_p_sel['nK_m3'].std('mnum').item()

print('nK_mean = {:.2e}, nK_std = {:.2e}'.format(nK_mean, nK_std))


#%%[markdown]

# Comparison between the three methods

#%%

ds_p = xr.concat(dss_p, 'method')
ds_p_stderr = xr.concat(dss_p_stderr, 'method')

#%%

ds_p['nK_m3'].plot(marker='o', hue='method')

plt.xscale('log')
plt.yscale('log')

#%%

plot.xr_errorbar(ds_p['nK_m3'], ds_p_stderr['nK_m3'], huedim='method')

#%%

ds_p_stderr['nK_m3'].plot(marker='o', hue='method')

plt.yscale('log')