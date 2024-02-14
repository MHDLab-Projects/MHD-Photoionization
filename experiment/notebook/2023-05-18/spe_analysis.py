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
fig, axes = plt.subplots(nrows=2, figsize=(12,4), sharex=True, sharey=True)

# Select the gatedelays
gatedelays = [800, 850]

# Plot each gatedelay on a different axis
for ax, gatedelay in zip(axes, gatedelays):
    da_plot = da.mean('estime').sel(gatedelay=gatedelay)
    da_plot.coords['gatedelay'].attrs['units'] = 'ns'
    da_plot.plot(ax=ax)

axes[0].set_xlabel('')
axes[1].set_xlabel('')


plt.savefig(pjoin(DIR_FIG_OUT, 'spe_gatedelay_img.png'), bbox_inches='tight', dpi=300)

#%%
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()

fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

dst_coords

#%%

from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

dst_coords['motor'].sel(time=slice(da.coords['estime'].min(), da.coords['estime'].max())).plot()

# assign_signal(da, dst_coords['motor'])
#%%

tw_536_pos_1 = slice(Timestamp('2023-05-18 22:23:34.172281344'), Timestamp('2023-05-18 22:30:20.482833664'), None)
tw_536_pos_2 = slice(Timestamp('2023-05-18 22:45:57.737082368'), Timestamp('2023-05-18 22:52:50.292869120'), None)
# tw_536_pos_both = slice(Timestamp('2023-05-18 22:21:30.395011328'), Timestamp('2023-05-18 22:57:37.093530624'), None)

da_536_pos =  xr.concat([
    da.sel(estime=tw_536_pos_1),
    da.sel(estime=tw_536_pos_2),
],
    dim='estime'
)


da_536_pos = assign_signal(da_536_pos, dst_coords['motor'], timeindex='estime')

da_536_pos = unstack_multindexed_acq_dim(da_536_pos, 'estime')['counts']

da_536_pos

#%%

da_536_pos.mean('mnum').sel(gatedelay=800).plot(col='motor', col_wrap=4, robust=True)


#%%

tw_536_fullpow = slice(Timestamp('2023-05-18 23:09:06.592949504'), Timestamp('2023-05-18 23:10:35.851570944'), None)

da_536_fullpow = da.sel(estime=tw_536_fullpow).mean('estime')

da_536_fullpow.sel(gatedelay=[800,850]).plot(col='gatedelay', robust=True)



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