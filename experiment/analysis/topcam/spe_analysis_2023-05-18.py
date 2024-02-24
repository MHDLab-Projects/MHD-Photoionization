#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

from calib_utils import calibrate_da_pimax_simple

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

da = das['21'].copy()

da = calibrate_da_pimax_simple(da)

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
fig, axes = plt.subplots(nrows=2, figsize=(8,4), sharex=True, sharey=True)

# Select the gatedelays
gatedelays = [800, 850]

# Plot each gatedelay on a different axis
for ax, gatedelay in zip(axes, gatedelays):
    da_plot = da.mean('estime').sel(gatedelay=gatedelay)
    da_plot.coords['gatedelay'].attrs['units'] = 'ns'
    da_plot.plot(ax=ax)

    ax.set_xlim(-10,220)

axes[0].set_xlabel('')
axes[0].set_ylabel('y (mm)')
axes[1].set_xlabel('x (mm)')
axes[1].set_ylabel('y (mm)')


plt.savefig(pjoin(DIR_FIG_OUT, 'spe_gatedelay_img_2023-05-18.png'), bbox_inches='tight', dpi=300)

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

#%%[markdown]

# # perspective fix and compare to cfd

#%%


import pickle
from calib_utils import pipe_transform_projective, CORRECTED_IMG_SCALE, LOC_BARREL_EXIT, LOC_BARREL_CENTERLINE

fp_transform = pjoin(DIR_DATA_OUT, 'tform_projective.pkl')
with open(fp_transform, 'rb') as f:
    tform = pickle.load(f)

tform

gatedelay_on = slice(797,801)
gatedelay_off = slice(700,793)

#%%
fp_cfd_2d = pjoin(REPO_DIR, 'modeling','cfd','analysis','output', 'mdot0130_phi080_K100.csv')


df_cfd = pd.read_csv(fp_cfd_2d, index_col=[0,1])

ds_cfd = df_cfd.to_xarray()

ds_cfd['pos_x'] = (ds_cfd['pos_x'] * 1000) - 208
ds_cfd['pos_y'] = ds_cfd['pos_y'] * 1000

# ds_cfd = ds_cfd.sel(pos_x=slice(200,))
# ds_cfd = ds_cfd.sel(pos_y=slice())

ds_cfd

#%%

ds_cfd['Yeq_KOH'].plot(x = 'pos_x')



#%%[markdown]

# ## 126 filesize. taken at 526_pos_2

#%%

da = das['126']


# Timewindow for 2023-05-18 526_pos_2
tw = slice(Timestamp('2023-05-18 22:46:08.092976128'), Timestamp('2023-05-18 22:52:35.143241984'), None)

da = da.sel(estime=tw)

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%


da_sel = da.mean('estime').sel(gatedelay=slice(790,820))

da_sel.plot(col='gatedelay', col_wrap=4)
# %%

#%%
#TODO: not have to do this
ds_sel  =  da_sel.to_dataset()

downsel_range = {'x': slice(100, 1024), 'y': slice(420, 620)}
da_sel_tf = pipe_transform_projective(ds_sel, tform, downsel_range=downsel_range)

#%%

da_sel_tf.plot(col='gatedelay', col_wrap=4)
# %%

da_sel_tf_on = da_sel_tf.sel(gatedelay=gatedelay_on).mean('gatedelay')
da_sel_tf_on.name = 'las_on'
da_sel_tf_off = da_sel_tf.sel(gatedelay=gatedelay_off).mean('gatedelay')
da_sel_tf_off.name = 'las_off'

ds = xr.merge([da_sel_tf_on, da_sel_tf_off])


ds.to_array('var').plot(col='var')

#%%

ds['las_off'].plot()

plt.xlim(0,150)

#%%

ds_cfd['p'].plot(x = 'pos_x', alpha=0.5)

plt.xlim(0,150)
#%%
fig, axes = plt.subplots(2,1, figsize=(10,4), sharex=True)

ds['las_off'].sel(y=0).plot(ax=axes[0])

axes[0].set_ylabel('Top Cam Intensity no Laser')

ds_cfd['p'].sel(pos_y=0, method='nearest').plot(color='r', ax=axes[1])

axes[1].set_ylabel('CFD pressure [Pa]')

plt.xlim(0,100)

plt.xlabel('x (mm)')

# da_sel_tf.sel(gatedelay=[790, 798]).plot(col='gatedelay')
# %%

diff = ds['las_on'] - ds['las_off']

diff.plot()

#%%


# %%

ds_cfd['Yeq_KOH'].plot(x = 'pos_x')
diff.plot(alpha=0.5)

plt.xlabel('x (mm)')
plt.ylabel('y (mm)')

plt.title('CFD vs. Topcam, 2023-05-18, 526_pos_2')

#%%[markdown]


# ## 21 filesize. taken at 526_pos_1

# %%

da = das['21']

# Timewindow for 2023-05-18 526_pos_1
tw = slice(Timestamp('2023-05-18 22:23:34.258497280'), Timestamp('2023-05-18 22:30:20.695414272'), None)
da = da.sel(estime=tw)

# da

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da_sel = da.mean('estime')

#TODO: not have to do this
ds_sel  =  da_sel.to_dataset()

downsel_range = {'x': slice(100, 1024), 'y': slice(420, 620)}
da_sel_tf = pipe_transform_projective(ds_sel, tform, downsel_range=downsel_range)

#%%

da_sel_tf.plot(col='gatedelay', col_wrap=4)
# %%

da_sel_tf_on = da_sel_tf.sel(gatedelay=gatedelay_on).mean('gatedelay')
da_sel_tf_on.name = 'las_on'
da_sel_tf_off = da_sel_tf.sel(gatedelay=gatedelay_off).mean('gatedelay')
da_sel_tf_off.name = 'las_off'

ds = xr.merge([da_sel_tf_on, da_sel_tf_off])


ds.to_array('var').plot(col='var')

#%%
diff = ds['las_on'] - ds['las_off']

diff.plot()

#%%

ds_cfd['Yeq_KOH'].plot(x = 'pos_x')
diff.plot(alpha=0.5)

plt.xlabel('x (mm)')
plt.ylabel('y (mm)')

plt.title('CFD vs. Topcam, 2023-05-18, 536_pos_1')

plt.ylim(-30,30)
plt.xlim(0,280)

# %%

#%%


ds['las_off'].plot()


plt.ylim(-30,30)

plt.xlim(0,280)
#%%

ds_cfd['Yeq_K'].plot(x = 'pos_x')

plt.ylim(-30,30)
plt.xlim(0,280)


#%%

ds_cfd['p'].plot(x = 'pos_x', alpha=0.5)

plt.xlim(0,150)
#%%

fig, axes = plt.subplots(2,1, figsize=(10,4), sharex=True)

ds['las_off'].sel(y=0).plot(ax=axes[0])

axes[0].set_ylabel('Top Cam Intensity no Laser')

ds_cfd['p'].sel(pos_y=0, method='nearest').plot(color='r', ax=axes[1])

axes[1].set_ylabel('CFD pressure [Pa]')

plt.xlim(0,100)

plt.xlabel('x (mm)')

for ax in axes:
    ax.set_title('')
# %%
