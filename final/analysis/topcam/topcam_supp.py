#%%
from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu
import matplotlib.ticker as ticker

create_standard_folders()

# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

# da_sel_tf.to_netcdf(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))

DIR_DATA_OUT = pjoin(REPO_DIR, 'final', 'dataset', 'topcam', 'output')
da_sel_tf = xr.load_dataarray(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))

beam_yslice = slice(-21,21)
beam_xslice = slice(160,190)

#%%
from mhdlab.pyvista_utils import CFDDatasetAccessor
fp_cfd_2d = pjoin(REPO_DIR, 'final','dataset','output', 'mdot0130_phi080_K100.csv')

df_cfd = pd.read_csv(fp_cfd_2d, index_col=[0,1])

ds_cfd = df_cfd.to_xarray()

ds_cfd['pos_x'] = (ds_cfd['pos_x'] * 1000) - 208
ds_cfd['pos_y'] = ds_cfd['pos_y'] * 1000


ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd

#%%

da_sel = da_sel_tf.sel(x=beam_xslice, y=beam_yslice)

da_sel.mean('gatedelay').plot()

plt.xlabel('x [mm]')
plt.ylabel('y [mm]')


#%%
plt.figure()

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

for ax in g.axs[:,0]:
    ax.set_ylabel('y [mm]')
for ax in g.axs[-1,:]:
    ax.set_xlabel('x [mm]')

# remove grid
for ax in g.axs.flatten():
    ax.grid(False)
    ax.xaxis.grid(False, which='both')
    ax.yaxis.grid(False, which='both')


plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_img_gatedelay.png'))




#%%
plt.figure()

da_sel2 = da_sel_tf.sel(x=beam_xslice, y=beam_yslice)

da_sel2.attrs['long_name'] = 'Counts'

da_sel2.mean('gatedelay', keep_attrs=True).plot(vmin=0)


plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_laserspot_zoom.png'))

