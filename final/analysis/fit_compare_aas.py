#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from pi_paper_utils.constants import *

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')

ds_aas_53x = xr.open_dataset(pjoin(data_directory, '53x_ds_aas.cdf')).xr_utils.stack_run()
# ds_aas_53x = ds_aas_53x.isel(run=[-1])
# ds_aas_53x = ds_aas_53x.sel(kwt= [0.1,1], method='nearest')
ds_aas_pos = xr.open_dataset(pjoin(data_directory, 'ds_pos_aas.cdf')).xr_utils.stack_run()
# ds_aas_pos = ds_aas_pos.isel(run=[-1])
# ds_aas_pos = ds_aas_pos.sel(motor=[50, 180 ], method='nearest')

goldi_pos = Quantity(180, 'mm')

#%%

ds_aas_pos


# %%
ds_cfd_beam = ppu.fileio.load_cfd_beam(kwt_interp=ds_aas_53x.coords['kwt'].values)

ds_cfd_beam = ds_cfd_beam.sel(phi=0.8).sel(offset=0)

da_cfd_beam = ds_cfd_beam[CFD_K_SPECIES_NAME]
da_cfd_beam = da_cfd_beam/da_cfd_beam.max('dist')

da_cfd_beam
# %%


#%%

from mhdlab.analysis.aas.fitting import AlphaFitPipeline2, AlphaFitPipelineNum1

def perform_fit_sequence(ds_fit):
    pipe_fit_alpha_2 = AlphaFitPipeline2()
    pipe_fit_alpha_num_1 = AlphaFitPipelineNum1()

    dss_p = []
    dss_fit_alpha = []

    ds_fit_alpha, ds_p_aas, ds_p_stderr_aas = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
    dss_p.append(ds_p_aas.assign_coords(method='num_1'))
    dss_fit_alpha.append(ds_fit_alpha.assign_coords(method='num_1'))

    pipe_fit_alpha_num_1.fit_prep_kwargs['led_off_norm_cutoff'] = None

    ds_fit_alpha, ds_p_aas, ds_p_stderr_aas = pipe_fit_alpha_num_1(ds_fit, perform_fit_kwargs={'nK_profile':da_cfd_beam_sel})
    dss_p.append(ds_p_aas.assign_coords(method='num_1_full'))
    dss_fit_alpha.append(ds_fit_alpha.assign_coords(method='num_1_full'))

    ds_fit_alpha, ds_p_aas, ds_p_stderr_aas = pipe_fit_alpha_2(ds_fit)
    dss_p.append(ds_p_aas.assign_coords(method='orig_pipe2'))
    dss_fit_alpha.append(ds_fit_alpha.assign_coords(method='orig_pipe2'))

    ds_p = xr.concat(dss_p, 'method')
    ds_fit_alpha = xr.concat(dss_fit_alpha, 'method')

    return ds_fit_alpha, ds_p
# %%[markdown]

# # 53x

#%%


dss_p = []

da_cfd_beam_sel = da_cfd_beam.sel(motor=goldi_pos,method='nearest')

ds_fit = ds_aas_53x#.isel(run=-1)
ds_fit = ds_fit.sel(kwt= da_cfd_beam.kwt.values) #TODO: downselecting as we don't have cfd for all kwt. Remove once we do

ds_fit_alpha, ds_p = perform_fit_sequence(ds_fit)

# ds_fit_alpha.to_array('var').mean('run').plot(row='kwt', col='mp', hue='var', yscale='log', ylim=(1e-3,1.1), figsize=(10,10))

#%%

ds_p['nK_m3'].mean('run').plot.line(x='kwt', col='mp', hue='method', marker='o')

plt.yscale('log')

#%%

da_plot = ds_fit_alpha.to_array('var').mean('run').sel(mp='barrel')

da_plot

da_plot.plot(row='kwt', col='method', hue='var', yscale='log', ylim=(1e-3,1.1), figsize=(10,10))

plt.xlim(760,775)


#%%


da_plot = ds_fit_alpha.to_array('var').mean('run').sel(mp='mw_horns')

da_plot

da_plot.plot(row='kwt', col='method', hue='var', yscale='log', ylim=(1e-3,1.1), figsize=(10,10))

plt.xlim(760,775)

#%%[markdown]

# # Position


#%%



da_cfd_beam_sel = da_cfd_beam.sel(kwt=1, method='nearest')

ds_fit = ds_aas_pos

#TODO: motor coordinates not exactly the same, why? 
da_cfd_beam_sel = da_cfd_beam_sel.sel(motor=ds_fit.coords['motor'].values, method='nearest')
da_cfd_beam_sel = da_cfd_beam_sel.assign_coords(motor=ds_fit.coords['motor'].values) 


ds_fit_alpha, ds_p = perform_fit_sequence(ds_fit)

# %%

# ds_p = ds_p.mean('run')
g = ds_p['nK_m3'].mean('run').plot.line(x='motor', col='mp', row='phi', hue='method', marker='o')

dropna(g)

plt.yscale('log')
# %%

da_plot = ds_fit_alpha.to_array('var').sel(mp='mw_horns').mean('run').sel(phi=0.8, method='nearest').dropna('motor',how='all')

da_plot.plot(row='motor', col='method', hue='var', yscale='log', ylim=(1e-3,1.1), figsize=(10,10))

plt.xlim(760,775)
# %%

da_plot = ds_fit_alpha.to_array('var').sel(mp='barrel').mean('run')

da_plot = ds_fit_alpha.to_array('var').sel(mp='barrel').mean('run').sel(phi=0.8, method='nearest').dropna('motor',how='all')
plt.xlim(760,775)