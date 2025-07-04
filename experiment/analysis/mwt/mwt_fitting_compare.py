
#%%[markdown]

# # Microwave Transmission fit comparison

#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

from mhdlab.analysis import mwt
from mhdlab.analysis.mwt.fitting import gen_model_dnedt

dss_p = []

#%%

ds_lecroy = ppu.fileio.load_lecroy('53x', avg_mnum=True, AS_calc='relative')




# %%

da_sel = ds_lecroy['dAS_rel']#.sel(kwt=1, method='nearest')
da_fit = da_sel
da_fit.plot(hue='run_plot', x='time', row='kwt')

# plt.yscale('log')

#%%
#%%[markdown]

# # pipe3
#%%

from mhdlab.analysis.mwt.fitting import pipe_fit_mwt_3

pipe_fit_mwt_3.fit_prep_kwargs.update({'interpolate_time':Quantity(0.01, 'us'), 'coarsen_amount':10})

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_3(da_fit)

ds_p['decay'] = 1/ds_p['krm']


ds_p2 = xr.merge([
    ds_p['decay'].rename('mean'),
    # ds_p_stderr['decay'].rename('std')
    ])

dss_p.append(ds_p2.assign_coords(method='pipe_3'))

#%%

ds_p['decay'].mean('run').plot()

plt.yscale('log')
plt.ylim(1, 1e2)

#%%

ds_mwt_fit['AS_all'].mean('run').plot(hue='kwt', marker='o')

plt.xlim(-1,1)


#%%

ds_mwt_fit.mean('run').to_array('var').plot(hue='var', row='kwt')

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

from mhdlab.analysis.mwt.fitting import pipe_fit_mwt_2

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_mwt_2(da_fit)

# ds_mwt_fit = xr.merge([ds_mwt_fit, da_fit.rename('AS_all')])
#%%
ne0=Quantity(2e12, 'cm**-3').magnitude
ds_p['decay'] = 1/(2*ne0*ds_p['krb'])

#%%

ds_p

#%%

ds_mwt_fit.mean('run').to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

#%%
ds_p2 = xr.merge([
    ds_p['decay'].rename('mean'),
    # ds_p_stderr['tau'].rename('std')
    ])

dss_p.append(ds_p2.assign_coords(method='pipe_2'))


#%%[markdown]

# ## exponential fit

#%%

from mhdlab.analysis.mwt.fitting import pipe_fit_exp


da_fit = da_sel
# da_fit = da_sel.copy()

ds_mwt_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)



#%%

ds_mwt_fit.mean('run').to_array('var').plot(row='kwt', hue='var')

plt.yscale('log')

plt.ylim(1e-4, )
plt.xlim(-1,40)

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

from mhdlab.plot import xr_errorbar

ds_plot = ds_p.mean('run')

# xr_errorbar(ds_plot['mean'], ds_plot['std'], huedim='method')
ds_plot['mean'].plot(hue='method', marker='o')

plt.ylim(1,20)

plt.yscale('log')
plt.xscale('log')

# %%



da_plot = ds_p.mean('run').drop('pipe_2', 'method')['mean']


da_plot.plot(hue='method', marker='o')

plt.ylim(3,35)

plt.yscale('log')
plt.xscale('log')

plt.legend(['Diff. Eq.', 'Exponential'], title='Method')

plt.ylabel('Decay constant (us)')
plt.xlabel("Kwt [%]")

plt.title('')

plt.savefig(pjoin(DIR_FIG_OUT, 'mwt_fitting_compare.png'))

# %%
