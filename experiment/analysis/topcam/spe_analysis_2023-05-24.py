#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils
create_standard_folders()

from calib_utils import calibrate_da_pimax_simple

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# %%

munged_dir = pjoin(REPO_DIR, 'experiment', 'data', 'munged', datestr)
spe_dir = pjoin(munged_dir, 'spe')

das = {}
names = ['4', '84']

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

da = das['4']

# da.mean('estime').plot(col='gatedelay', col_wrap=4)
da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da.mean('estime').plot(col='gatedelay')


#%%

da = das['84'].copy()

da.mean('estime').mean('gatedelay').plot(robust=True)
#%%
# Timewindow for these files is entire seedramp
da['estime'].plot()

#%%

# downselect to 1% kwt time window
tw = slice(Timestamp('2023-05-24 22:25:55.952985600'), Timestamp('2023-05-24 22:28:12.304348160'), None) 

da_sel = da.sel(estime=tw).sel(gatedelay=slice(780,810))

da_sel = da_sel.mean('estime')

da_sel

#%%
import pickle
from calib_utils import pipe_transform_projective, CORRECTED_IMG_SCALE, LOC_BARREL_EXIT, LOC_BARREL_CENTERLINE

fp_transform = pjoin(DIR_DATA_OUT, 'tform_projective.pkl')
with open(fp_transform, 'rb') as f:
    tform = pickle.load(f)

da_sel_2 = da_sel.to_dataset()
downsel_range = {'x': slice(100, 1024), 'y': slice(420, 620)}
da_sel_tf = pipe_transform_projective(da_sel_2, tform, downsel_range=downsel_range)

da_sel_tf = da_sel_tf.sel(y=slice(-30,30))

#%%
da_sel_tf.mean('gatedelay').plot()


#%%

beam_yslice = slice(-21,21)
beam_xslice = slice(160,190)

da_sel = da_sel_tf.sel(x=beam_xslice, y=beam_yslice)

da_sel.mean('gatedelay').plot()

plt.xlabel('x (mm)')
plt.ylabel('y (mm)')


#%%

da_sel.mean(['x', 'y']).plot(marker='o')
# %%

#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors


# Define the overload threshold
overload_threshold = 1100

# Create a colormap that has a different color for overloaded values
cmap = plt.get_cmap('viridis')
colors_under_overload = cmap(np.linspace(0, 1, overload_threshold))
overload_color = np.array([[1, 0, 0, 1]])  # Red color for overload
new_colors = np.vstack((colors_under_overload, overload_color))
new_cmap = colors.LinearSegmentedColormap.from_list('new_cmap', new_colors)

# Create a normalization that maps the minimum value in the data to the start of the colormap,
# and values above the overload threshold to the end of the colormap
norm = colors.Normalize(vmin=np.min(da_sel.values), vmax=overload_threshold)


da_sel = da_sel_tf.sel(gatedelay=slice(780,802))

g = da_sel.plot(
    col='gatedelay', 
    col_wrap=4, 
    figsize=(8,4), 
    cmap = new_cmap,
    norm=norm,
    # robust=True
    # vmax=1200
    )  

for ax in g.axes[:,0]:
    ax.set_ylabel('y (mm)')
for ax in g.axes[-1,:]:
    ax.set_xlabel('x (mm)')

# remove grid
for ax in g.axes.flatten():
    ax.grid(False)
    ax.xaxis.grid(False, which='both')
    ax.yaxis.grid(False, which='both')


plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_img_gatedelay.png'))

#%%

offset = da_sel_tf.sel(gatedelay=780).sel(y=slice(25,30)).mean('x').mean('y').item()

offset


#%%
fig, axes = plt.subplots(ncols=2, figsize=(8,2), sharex=True, sharey=True)

# Select the gatedelays

gatedelays = [790, 810]

# Plot each gatedelay on a different axis

for ax, gatedelay in zip(axes, gatedelays):
    da_plot = da_sel_tf.sel(gatedelay=gatedelay) - offset
    da_plot.attrs['long_name'] = 'Counts'
    da_plot.coords['gatedelay'].attrs['units'] = 'ns'
    da_plot.plot(ax=ax, vmin=0)

    ax.set_xlim(-10,220)


plt.savefig(pjoin(DIR_FIG_OUT, 'spe_gatedelay_img_2023-05-24.png'), bbox_inches='tight', dpi=300)
# set the 

# %%

da_sel2 = da_sel.sel(x=beam_xslice, y=beam_yslice)

da_sel2.mean('gatedelay').plot()

plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_laserspot_zoom.png'))

#%%

da_sel2.mean(['x', 'y']).mean('estime').plot(marker='o')

plt.ylabel('Counts')
plt.xlabel('Gate Delay (ns)')

plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_timedecay.png'))


# %%
