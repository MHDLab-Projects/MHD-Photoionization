#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_absem_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_absem.cdf')).xr_utils.stack_run()
ds_absem_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_absem.cdf')).xr_utils.stack_run()

goldi_pos = Quantity(180, 'mm')


# %%
ds_cfd_beam = ppu.fileio.load_cfd_beam(kwt_interp=ds_absem_53x.coords['kwt'].values)

ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0)

da_cfd_beam = ds_cfd_beam['Yeq_K']
da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

da_cfd_beam
# %%


da_cfd_beam.sel(motor=goldi_pos,method='nearest').plot(col='mp', hue='kwt')

plt.yscale('log')

plt.ylim(1e-7,)
plt.xlim(3,8)
#%%

da_cfd_beam.sel(kwt=1, method='nearest').plot(col='mp', hue='motor')


plt.yscale('log')

plt.ylim(1e-7,)
plt.xlim(3,8)
# %%[markdown]

# # 53x

#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_num_1, pipe_fit_alpha_2

dss_p = []

da_cfd_beam_sel = da_cfd_beam.sel(motor=goldi_pos,method='nearest')

ds_fit = ds_absem_53x#.isel(run=-1)
ds_fit = ds_fit.sel(kwt= da_cfd_beam.kwt.values) #TODO: downselecting as we don't have cfd for all kwt. Remove once we do

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
dss_p.append(ds_p_absem.assign_coords(method='num_1'))

pipe_fit_alpha_num_1.fit_prep_kwargs['led_off_norm_cutoff'] = None

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
dss_p.append(ds_p_absem.assign_coords(method='num_1_full'))

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_2(ds_fit)
dss_p.append(ds_p_absem.assign_coords(method='orig_pipe2'))


# ds_fit_alpha.to_array('var').mean('run').plot(row='kwt', col='mp', hue='var', yscale='log', ylim=(1e-3,1.1), figsize=(10,10))

#%%

ds_p = xr.concat(dss_p, 'method')
ds_p['nK_m3'].mean('run').plot.line(x='kwt', col='mp', hue='method', marker='o')

plt.yscale('log')

#%%[markdown]

# # Position


#%%


dss_p = []

da_cfd_beam_sel = da_cfd_beam.sel(kwt=1, method='nearest')

ds_fit = ds_absem_pos.isel(run=-1)

#TODO: motor coordinates not exactly the same, why? 
da_cfd_beam_sel = da_cfd_beam_sel.sel(motor=ds_fit.coords['motor'].values, method='nearest')
da_cfd_beam_sel = da_cfd_beam_sel.assign_coords(motor=ds_fit.coords['motor'].values) 

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
dss_p.append(ds_p_absem.assign_coords(method='num_1'))

pipe_fit_alpha_num_1.fit_prep_kwargs['led_off_norm_cutoff'] = None

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
dss_p.append(ds_p_absem.assign_coords(method='num_1_full'))

ds_fit_alpha, ds_p_absem, ds_p_stderr_absem = pipe_fit_alpha_2(ds_fit)
dss_p.append(ds_p_absem.assign_coords(method='orig_pipe2'))
# %%

ds_p = xr.concat(dss_p, 'method')
# ds_p = ds_p.mean('run')
g = ds_p['nK_m3'].plot.line(x='motor', col='mp', row='phi', hue='method', marker='o')

dropna(g)

plt.yscale('log')
# %%


ds_fit