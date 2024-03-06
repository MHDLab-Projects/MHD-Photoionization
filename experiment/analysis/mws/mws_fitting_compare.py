
#%%[markdown]

# # Microwave Scattering fit comparison

#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.analysis.mws.fitting import gen_model_dnedt

dss_p = []

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x')

#Make a general mws method? 
counts = ds_lecroy['AS'].isel(time=0).count('mnum')
ds_lecroy = ds_lecroy.where(counts > 100)


# %%

da_sel = ds_lecroy['AS']#.sel(kwt=1, method='nearest')
da_fit = da_sel.mean('mnum')


da_fit.plot(hue='run_plot', x='time', row='kwt')

plt.yscale('log')

#%%[markdown]

# # pipe3

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_3
ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_3(da_fit, take_log=False)

ds_p['decay'] = 1/ds_p['krm']


ds_p2 = xr.merge([
    ds_p['decay'].rename('mean'),
    # ds_p_stderr['decay'].rename('std')
    ])

dss_p.append(ds_p2.assign_coords(method='pipe_3'))


#%%

ds_mws_fit.mean('run').to_array('var').plot(hue='var', row='kwt')

plt.yscale('log')

plt.xlim(-1,30)

plt.ylim(1e-3,)

#%%

dne = ds_p['dne'].pint.quantify('particle/um**3').pint.to('particle/cm**3')

dne.mean('run').pint.magnitude

#%%[markdown]

#%%[markdown]

# ## Pipe 2: 

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_mws_2

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_mws_2(da_fit, take_log=False)

# ds_mws_fit = xr.merge([ds_mws_fit, da_fit.rename('AS_all')])

#%%

ds_mws_fit.mean('run').to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

#%%
ds_p2 = xr.merge([
    ds_p['tau'].rename('mean'),
    # ds_p_stderr['tau'].rename('std')
    ])

dss_p.append(ds_p2.assign_coords(method='pipe_2'))


#%%[markdown]

# ## exponential fit

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp


da_fit = da_sel.mean('mnum')
# da_fit = da_sel.copy()

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit, method='iterative', fit_timewindow=slice(Quantity(5, 'us'),Quantity(15, 'us')))



#%%

ds_mws_fit.mean('run').to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

#%%

ds_p

#%%

ds_p2 = xr.merge([
    ds_p['decay'].rename('mean'),
    ds_p['decay'].where(False).fillna(0).rename('std') # kr std not defined yet for exponential fit...
    ])

dss_p.append(ds_p2.assign_coords(method='exp'))

#%%

ds_p = xr.concat(dss_p, dim='method')

ds_p

#%%

# Define your custom plotting function
def custom_plot(x, mean, std, **kwargs):
    x = [str(t) for t in x]
    plt.errorbar(x, mean, yerr=std, capsize=4, marker='o')


fg = xr.plot.FacetGrid(ds_p, row='run', col='method')

# Apply the custom function to each facet

fg.map(custom_plot, 'kwt', 'mean', 'std')


#%%

# ds_p_plot = ds_p.unstack('run').stack

fig, axes = plt.subplots(len(ds_p.indexes['run']), figsize=(5,15))

for i, (method, ds) in enumerate(ds_p.groupby('run')):

    ax=axes[i]

    ds['mean'].plot(hue='method', ax=ax, marker='o')

# %%

from mhdpy.plot.common import xr_errorbar

ds_plot = ds_p.mean('run')

# xr_errorbar(ds_plot['mean'], ds_plot['std'], huedim='method')
ds_plot['mean'].plot(hue='method', marker='o')

plt.ylim(1,20)

plt.yscale('log')
plt.xscale('log')

# %%



da_plot = ds_p.mean('run').drop('pipe_2', 'method')['mean']


da_plot.plot(hue='method', marker='o')

plt.ylim(3,15)

plt.yscale('log')
plt.xscale('log')

plt.legend(['Diff. Eq.', 'Exponential'], title='Method')

plt.ylabel('Decay constant (us)')
plt.xlabel("Kwt [%]")

plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'mws_fitting_compare.png'))
