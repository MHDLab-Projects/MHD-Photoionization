#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-05-18'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# %%

munged_dir = pjoin(REPO_DIR, 'experiment', 'data', 'munged', datestr)
spe_dir = pjoin(munged_dir, 'spe')

das = {}
names = ['21', '63', '126']

# %%

for name in names:

    fp_in = pjoin(spe_dir, 'PI_topcam_{}.cdf'.format(name))

    ds = xr.load_dataset(fp_in)
    da = ds['counts']

    das[name] = da

#%%

import pandas as pd

df = pd.DataFrame(index=das.keys())

for name, da in das.items():
    # coord_dict[name]['gatedelay'].plot(label=name)
    coords = da.coords
    gatedelay = coords['gatedelay']

    start = gatedelay[0]
    step = gatedelay[1] - gatedelay[0]
    num = len(gatedelay)

    df.loc[name, 'start'] = start
    df.loc[name, 'step'] = step
    df.loc[name, 'num'] = num

df

#%%

for name, da in das.items():

    es_time = da['estime']
    
    es_time.plot(label=name, marker='o')

plt.legend()


#%%

da = das['21']

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da.mean('estime').plot(col='gatedelay', col_wrap=4)

#%%

plt.figure(figsize=(5,4))

da_plot = da.mean('estime').sel(gatedelay=[750, 800, 850])
da_plot.coords['gatedelay'].attrs['units'] = 'ns'

da_plot.plot(row='gatedelay')

#%%

# Create a figure and axes with the desired size
fig, axes = plt.subplots(nrows=3, figsize=(4, 8), sharex=True, sharey=True)

# Select the gatedelays
gatedelays = [750, 800, 850]

# Plot each gatedelay on a different axis
for ax, gatedelay in zip(axes, gatedelays):
    da_plot = da.mean('estime').sel(gatedelay=gatedelay)
    da_plot.coords['gatedelay'].attrs['units'] = 'ns'
    da_plot.plot(ax=ax)

axes[0].set_xlabel('')
axes[1].set_xlabel('')


plt.savefig(pjoin(DIR_FIG_OUT, 'spe_gatedelay_img.png'), bbox_inches='tight', dpi=300)

#%%

da = das['63']

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da_sel = da.mean('estime').sel(x=slice(650,780), y=slice(450,575))

da_sel.mean('gatedelay').plot(robust=True)

#%%

da_sel.mean(['x','y']).plot(marker='o')

#%%

da = das['126']

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%


da_sel = da.mean('estime').sel(gatedelay=slice(790,820))

da_sel.plot(col='gatedelay', col_wrap=4)