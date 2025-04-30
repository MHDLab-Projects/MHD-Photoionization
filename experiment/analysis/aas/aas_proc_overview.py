# %% [markdown]
# # Atomic Absorption Data Processing Overview

# %%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

# enable dark mode
# plt.style.use('dark_background')

datestr = '2023-05-18'


# %% [markdown]
# ## Initial Munged Data
# 
# This data is the data after the first standardized post experiment data processing. 

# %% [markdown]
# Data is taken continuously with LED (AKA Shutter) continously switching. 

# %%
fp_1 = pjoin(REPO_DIR, 'experiment', 'data','munged',datestr, 'Munged','Spectral','aas.tdms')

ds = mhdlab.fileio.spectral.load_multindexed_spectral(fp_1)

ds

# %% [markdown]
# The Raw Data is divided by the integration time to give counts/ms. Other time-based metadata is dropped, but it is first checked that there is only one value. 

# %%
ds = mhdlab.fileio.spectral.load_aas(fp_1, convert_to_PT=False)

ds

# %%
ds['counts'].mean('time').plot(marker='o')

# %% [markdown]
# To reduce the size of the data, the datapoint frequency is reduced away from the peaks. This can be changed

# %%
ds['counts'].mean('time').plot(marker='o')
plt.xlim(760,780)

# %% [markdown]
# Time plots over experiment

# %%
ds['counts'].mean('wavelength').plot()


# %% [markdown]
# Zoom in showing LED switching

# %%

# ds_sel = ds.isel(time=slice(1000,1100))
ds_sel = ds.isel(time=slice(20000,20100))
tw = slice(ds_sel.coords['time'].values[0], ds_sel.coords['time'].values[-1])

fig, axes = plt.subplots(2,1)

ds_sel['counts'].mean('wavelength').plot(ax=axes[0], marker='o')
ds_sel['led'].plot(ax=axes[1], marker='o')


# %% [markdown]
# ## Additional Processing

# For datasets with multiplexer switching, this coordinate is added 
#%%

import json
from mhdlab.fileio import TFxr

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

has_multiplexer = settings['has_multiplexer']

ds_aas = ds

# Start processing 
from mhdlab.analysis.aas import calc_alpha_simple
from mhdlab.coords import reduce_acq_group, get_value_switches, downselect_num_acq
from mhdlab.coords.spectral import assign_multiplexer_coord

data_folder = os.path.join(REPO_DIR, 'experiment', 'data', 'munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

# Determine LED switching events
switches = get_value_switches(ds_aas.coords['led'].values, switch_to_vals=['led_off','led_on'])
ds_aas = ds_aas.assign_coords(led_switch_num=('time', switches))

if has_multiplexer:
    ds_mp = assign_multiplexer_coord(
    ds_aas,
    mp_sent=dsst['multiplexer_send']['Position'].pint.dequantify(),
    mp_receive=dsst['multiplexer_receive']['Position'].pint.dequantify(),
    mp_port_names={1:'barrel', 2:'mw_horns'}
    )

    # Now we remove data when the multiplexer was switching, kept to allow for accurate determination of switching events
    ds_aas = ds_mp.where(ds_mp['mp'] != 'switch').dropna('time', how='all')
else:
    ds_aas = ds_aas.assign_coords(mp = ('time', ['barrel']*len(ds_aas.coords['time']) ))

ds_aas = ds_aas.groupby('led_switch_num').apply(downselect_num_acq, num_acq=10)
ds_aas = ds_aas.dropna('time', how='all')


# %%

ds_aas
#%%


# ds_sel = ds.isel(time=slice(1000,1100))
ds_sel = ds_aas.sel(time=tw)

fig, axes = plt.subplots(4,1, figsize=(5,10), sharex=True)

ds_sel['counts'].mean('wavelength').plot(ax=axes[0], marker='o')
ds_sel['led'].plot(ax=axes[1], marker='o')
ds_sel['mp'].plot(ax=axes[2], marker='o')
ds_sel['led_switch_num'].plot(ax=axes[3], marker='o')

#%%[markdown]

# For now to reduce the size we group by each led switch and mp switch and average (averaging the 10 rapid-fire acquisitions). This should be able to be changed. 
# We also 'interpolate' (using nearest value) the led off data to the led on data such that the time data matches. 

#%%

from mhdlab.xr_utils import interp_ds_to_var

# Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
#TODO: remove this averaging, should only perform one average. But need to revisit data pipeline to avoid too many large files. 
acq_groups = ['led_switch_num','led','time','mp']
ds_acq_group = ds_aas.set_index(acq_group=acq_groups)

ds_reduce = reduce_acq_group(ds_acq_group)
ds_reduce = ds_reduce.reset_coords('led_switch_num', drop=True)

acq_groups.remove('led_switch_num')
ds_alpha = ds_reduce.set_index(temp=acq_groups).unstack('temp')
ds_alpha = ds_alpha['counts_mean']

ds_alpha = ds_alpha.to_dataset('led')
ds_alpha = interp_ds_to_var(ds_alpha, 'led_on')

ds_alpha
# %%

da = ds_alpha.sel(time=tw).to_array('var')

da.plot(hue='var', col='mp', row='time')
# %%

#%%[markdown]

# ## Calculate Alpha

# # Just pulling in the main dataset for now.

#%%


ds_aas = ppu.fileio.load_aas('53x')

#%%


from mhdlab.analysis.aas.fitting import pipe_fit_alpha_2

ds_fit = ds_aas.mean('mnum')

ds_aas_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

#%%

from mhdlab.analysis.aas.fit_prep import pipe_fit_prep_alpha_2

ds_fit_beta_off = ds_fit.aas.remove_beta_offset(beta_offset_wls=slice(750,755))

# da_fit_prep =  pipe_fit_prep_alpha_2(ds_fit, led_off_norm_cutoff=0.8)
#%%


#%%

#%%

ds_aas['alpha_fit'] = ds_aas_fit['alpha_fit']
ds_aas['alpha_beta_off'] = ds_fit_beta_off['alpha']
ds_aas['alpha_red'] = ds_aas_fit['alpha_red']


# %%

ds_sel = ds_aas.sel(run=('2023-05-24',1)).sel(kwt=1, method='nearest')
ds_sel = ds_sel.mean('mnum')
ds_sel



#%%

ds_sel2 = ds_sel.sel(mp='barrel')

nK_sel = ds_p['nK_m3'].sel(run=('2023-05-24',1)).sel(kwt=1, method='nearest').sel(mp='barrel')


fig, axes = plt.subplots(4, sharex=True, figsize=(4,9))

ds_sel2[['led_on','led_off']].to_array('var').plot(ax=axes[0], hue='var')

ds_sel2[['diff','calib']].to_array('var').plot(ax=axes[1], hue='var')

ds_sel2['alpha'].plot(ax=axes[2], label='Raw')
ds_sel2['alpha_beta_off'].plot(ax=axes[2], label='Offset corr.')
axes[2].legend(loc='upper right')

dropna(axes[2])

ds_sel2['alpha_red'].plot(ax=axes[3], label='Data')
ds_sel2['alpha_fit'].plot(ax=axes[3], label='Fit')

plt.text
labels = ['A)', 'B)', 'C)', 'D)']

for ax, label in zip(axes, labels):
    ax.set_title('')
    ax.text(-0.2, 1.1, label, transform=ax.transAxes, va='top', ha='right')

axes[0].set_ylabel('Counts')
axes[0].set_xlabel('')
axes[0].legend(['Led On','Led Off'])

axes[1].set_ylabel('Counts')
axes[1].set_xlabel('')
axes[1].legend(['On - Off','Calib.'])

axes[2].set_ylabel('Absorption')
axes[2].set_xlabel('')

axes[3].set_ylabel('Absorption')
axes[3].legend()
axes[3].text(750,0.5, '$n_K$ = {:.2e}'.format(nK_sel.values), fontsize=10)

axes[-1].set_xlabel('Wavelength (nm)')

plt.savefig(pjoin(DIR_FIG_OUT, 'aas_proc_overview.png'), bbox_inches='tight')

#%%

