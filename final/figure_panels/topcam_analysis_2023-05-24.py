#%%
from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu
create_standard_folders()

# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

# da_sel_tf.to_netcdf(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))

DIR_DATA_OUT = pjoin(REPO_DIR, 'final', 'dataset', 'topcam', 'output')
da_sel_tf = xr.load_dataarray(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))


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


#%%

from matplotlib_scalebar.scalebar import ScaleBar

plt.rcParams.update({'font.size': 8})

fig, axes = plt.subplots(ncols=2, nrows=3, figsize = (5,6))

fig.subplots_adjust(wspace= 0.45)




label_color = 'black' #axis label color

labels = ['A)','B)','C)','D)','E)','F)']

# for ax, label in zip(axes.flatten(), labels):
#     ax.text(-0.25, 1.0, label, transform=ax.transAxes, va='bottom', ha='left')
#     # ax.tick_params(axis='y', which = "both", colors='white')
    
#     ax.xaxis.label.set_color(label_color)
#     ax.yaxis.label.set_color(label_color)




for ax, label in zip(axes.flatten(), labels):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - .1, Y + .015, label)
    # ax.tick_params(axis='y', which = "both", colors='white')
    
    ax.xaxis.label.set_color(label_color)
    ax.yaxis.label.set_color(label_color)


# Select the gatedelays

gatedelays = [780, 790]

# Plot each gatedelay on a different axis

ax = axes[0,1]
da_plot = da_sel_tf.sel(gatedelay=gatedelays[0])
da_plot.plot(ax=ax, vmin=0)
axes[0,1].set_title('ICCD before [counts]')


ax = axes[1,1]
da_plot = da_sel_tf.sel(gatedelay=gatedelays[1]) 
da_plot.plot(ax=ax, vmin=0)
axes[1,1].set_title('ICCD during [counts]')

ax = axes[0,0]
ds_cfd[ppu.CFD_K_SPECIES_NAME].plot(ax=ax, vmin=0, x='pos_x')
axes[0,0].set_title(r'K $\mathrm{[\#/cm^3]}$')

ax = axes[1,0]
ds_cfd[ppu.CFD_KOH_SPECIES_NAME].plot(ax=ax, vmin=0, x='pos_x', vmax=2e16)
axes[1,0].set_title(r'KOH $\mathrm{[\#/cm^3]}$')


ax = axes[2,0]
ds_cfd['T'].plot(ax=ax, x='pos_x', vmax=3200)
axes[2,0].set_title('T [K]')



for ax in axes.flatten()[:-1]:
    ax.set_ylabel('')
    ax.set_xlabel('')
    # ax.set_xticks([])
    # ax.set_yticks([])

    ax.set_xlim(0,250)
    ax.set_ylim(-25,25)

    # get the colorbar of each axis

    if ax.collections:
        cbar = ax.collections[0].colorbar
        cbar.set_label('')

axes[0,0].set_xticks([])
axes[0,1].set_xticks([])
axes[1,0].set_xticks([])

axes[0,0].set_ylabel('y [mm]')
axes[1,0].set_ylabel('y [mm]')
axes[2,0].set_ylabel('y [mm]')
axes[2,0].set_xlabel('x [mm]')

# draw beam slice box on axes [1,1]

axes[1,1].add_patch(plt.Rectangle((160,-21), 30, 42, fill=False, edgecolor='r', lw=1, linestyle='--', alpha=0.5))


ax_lpf = axes[2,1]
ax_lpf_twin = ax_lpf.twinx()

# axes[2, 1].cla()
# axes[2, 1].axis('off')
# plt.savefig(pjoin(DIR_FIG_OUT, 'spe_cfd_compare_2ns_2023-05-24.png'), bbox_inches='tight', dpi=300)
# set the 


da_sel2 = da_sel_tf.sel(x=beam_xslice, y=beam_yslice)

da_sel2.attrs['long_name'] = 'Counts'


fp_laser_profile = pjoin(REPO_DIR, 'experiment', 'notebook', '2018-11-20', 'output', 'laser_profile_1.csv')

df_laser_profile = pd.read_csv(fp_laser_profile, index_col=0)

# df_laser_profile.plot()


ln_cam = da_sel2.mean(['x', 'y']).plot(marker='o', ax = ax_lpf)

laser_timeshift = 30



# for ax in axes.flatten():

#     ax.set_anchor('W')

df_laser_profile_shift = df_laser_profile.copy()
df_laser_profile_shift.index = df_laser_profile_shift.index + laser_timeshift
df_laser_profile_shift.plot(ax=ax_lpf_twin, color='r')
ln_laser = ax_lpf_twin.get_lines()

ax_lpf_twin.get_legend().remove()

lns = ln_cam + ln_laser
labs = ['Top Cam', 'Laser Profile']
ax_lpf.legend(lns, labs, loc=0)

ax_lpf.set_xlabel('Gate Delay [ns]')
ax_lpf.set_ylabel('ICCD [counts]')
ax_lpf_twin.set_ylabel('Laser Profile \n[counts]')

ax_lpf_twin.grid(False)

ax_lpf.set_xlim(780,830)

fig.align_ylabels(axes[:, 1])

# plt.savefig(pjoin(DIR_FIG_OUT, '536_iccd_timedecay.png'))

# ax_lpf.set_yscale('log')
# ax_lpf_twin.set_yscale('log')

# %%