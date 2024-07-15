# %% [markdown]
# # Absorption Emission Data Processing Overview

# %%
from mhdpy.analysis.standard_import import *
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
fp_1 = pjoin(REPO_DIR, 'experiment', 'data','munged',datestr, 'Munged','Spectral','absem.tdms')

ds = mhdpy.fileio.spectral.load_multindexed_spectral(fp_1)

ds

# %% [markdown]
# The Raw Data is divided by the integration time to give counts/ms. Other time-based metadata is dropped, but it is first checked that there is only one value. 

# %%
ds = mhdpy.fileio.spectral.load_absem(fp_1, convert_to_PT=False)

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
from mhdpy.fileio import TFxr

with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
    settings = json.load(f)[datestr]

has_multiplexer = settings['has_multiplexer']

ds_absem = ds

# Start processing 
from mhdpy.analysis.absem import calc_alpha_simple
from mhdpy.coords import reduce_acq_group, get_value_switches, downselect_num_acq
from mhdpy.coords.spectral import assign_multiplexer_coord

data_folder = os.path.join(REPO_DIR, 'experiment', 'data', 'munged',datestr)

dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

# Determine LED switching events
switches = get_value_switches(ds_absem.coords['led'].values, switch_to_vals=['led_off','led_on'])
ds_absem = ds_absem.assign_coords(led_switch_num=('time', switches))

if has_multiplexer:
    ds_mp = assign_multiplexer_coord(
    ds_absem,
    mp_sent=dsst['multiplexer_send']['Position'].pint.dequantify(),
    mp_receive=dsst['multiplexer_receive']['Position'].pint.dequantify(),
    mp_port_names={1:'barrel', 2:'mw_horns'}
    )

    # Now we remove data when the multiplexer was switching, kept to allow for accurate determination of switching events
    ds_absem = ds_mp.where(ds_mp['mp'] != 'switch').dropna('time', how='all')
else:
    ds_absem = ds_absem.assign_coords(mp = ('time', ['barrel']*len(ds_absem.coords['time']) ))

ds_absem = ds_absem.groupby('led_switch_num').apply(downselect_num_acq, num_acq=10)
ds_absem = ds_absem.dropna('time', how='all')


# %%

ds_absem
#%%


# ds_sel = ds.isel(time=slice(1000,1100))
ds_sel = ds_absem.sel(time=tw)

fig, axes = plt.subplots(4,1, figsize=(5,10), sharex=True)

ds_sel['counts'].mean('wavelength').plot(ax=axes[0], marker='o')
ds_sel['led'].plot(ax=axes[1], marker='o')
ds_sel['mp'].plot(ax=axes[2], marker='o')
ds_sel['led_switch_num'].plot(ax=axes[3], marker='o')

#%%[markdown]

# For now to reduce the size we group by each led switch and mp switch and average (averaging the 10 rapid-fire acquisitions). This should be able to be changed. 
# We also 'interpolate' (using nearest value) the led off data to the led on data such that the time data matches. 

#%%

from mhdpy.xr_utils import interp_ds_to_var

# Perform grouping operations over switching groups, to obtain one led off and on for each switch. 
#TODO: remove this averaging, should only perform one average. But need to revisit data pipeline to avoid too many large files. 
acq_groups = ['led_switch_num','led','time','mp']
ds_acq_group = ds_absem.set_index(acq_group=acq_groups)

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


ds_absem = ppu.fileio.load_absem('53x')

#%%


from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit = ds_absem.mean('mnum')

ds_absem_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

#%%

from mhdpy.analysis.absem.fit_prep import pipe_fit_prep_alpha_2

ds_fit_beta_off = ds_fit.absem.remove_beta_offset(beta_offset_wls=slice(750,755))

# da_fit_prep =  pipe_fit_prep_alpha_2(ds_fit, led_off_norm_cutoff=0.8)
#%%


#%%

#%%

ds_absem['alpha_fit'] = ds_absem_fit['alpha_fit']
ds_absem['alpha_beta_off'] = ds_fit_beta_off['alpha']
ds_absem['alpha_red'] = ds_absem_fit['alpha_red']


# %%

ds_sel = ds_absem.sel(run=('2023-05-24',1)).sel(kwt=1, method='nearest')
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

for ax in axes:
    ax.set_title('')

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


plt.savefig(pjoin(DIR_FIG_OUT, 'absem_proc_overview.png'), bbox_inches='tight')

#%%

