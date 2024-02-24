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

da = calibrate_da_pimax_simple(da)

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

#%%
# Timewindow for these files is entire seedramp
da['estime'].plot()

#%%

beam_yslice = slice(-21,21)
beam_xslice = slice(160,195)

da_sel = da.mean('estime').sel(x=beam_xslice, y=beam_yslice)

da_sel.mean('gatedelay').plot(robust=True)

plt.xlabel('x (mm)')
plt.ylabel('y (mm)')


#%%

da_sel.mean(['x', 'y']).plot(marker='o')
# %%
# downselect to 1% kwt time window
tw = slice(Timestamp('2023-05-24 22:25:55.952985600'), Timestamp('2023-05-24 22:28:12.304348160'), None) 

da_sel = da.sel(estime=tw).sel(gatedelay=slice(780,802))

da_sel

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


g = da_sel.mean('estime').plot(
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

# %%

da_sel2 = da_sel.sel(x=beam_xslice, y=beam_yslice)

da_sel2.mean('estime').mean('gatedelay').plot()

plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_laserspot_zoom.png'))

#%%

da_sel2.mean(['x', 'y']).mean('estime').plot(marker='o')

plt.ylabel('Counts')
plt.xlabel('Gate Delay (ns)')

plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_timedecay.png'))


# %%
