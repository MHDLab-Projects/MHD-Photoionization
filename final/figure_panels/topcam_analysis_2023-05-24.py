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

fp_laser_profile = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'laser_profile_1.csv')
laser_timeshift = -80
df_laser_profile = pd.read_csv(fp_laser_profile, index_col=0)
df_laser_profile_shift = df_laser_profile.copy()
df_laser_profile_shift.index = df_laser_profile_shift.index + laser_timeshift

df_laser_profile_shift = df_laser_profile_shift.loc[780:]


#%%

def scale_cbar_ticks(cbar, scale_factor, new_label = None, rounding = None):
    ''' scales tick labels of cbar by scale_factor (new labels = old labels * scale_factor) 
        and renames label to new_label if given'''

    cbar_ticks = cbar.get_ticks()
    if rounding is not None:
        new_ticks = [round(x, rounding) for x in cbar_ticks * scale_factor]
    else: 
        new_ticks = cbar_ticks * scale_factor

    cbar.set_ticklabels(new_ticks)

    if new_label is not None:
        cbar.set_label(new_label)


gatedelays = [780, 790]
label_color = 'black' #axis label color

labels = ['A)','B)','C)','D)','E)','F)']



#%%
from matplotlib.colors import Normalize

#this cell generates the full figure
#Use gridspec to precisely align axes

plt.rcParams.update({'font.size': 8})


norm = Normalize(vmin = 0, vmax = 1E16)

shared_colorbar = False

cbar_location = "right"


fig = plt.figure(figsize = (7.5, 7))
# fig.set_figwidth(7.5)

if shared_colorbar:
    height_ratios = [10,1,10,10]
    spec = fig.add_gridspec(4, 2, height_ratios = height_ratios)

else: 
    height_ratios = [10,10,10]
    spec = fig.add_gridspec(3, 2, height_ratios = height_ratios)



if shared_colorbar:

    ax = fig.add_subplot(spec[0, 0])
    plot_obj = ds_cfd[ppu.CFD_K_SPECIES_NAME].plot(ax=ax, x='pos_x', add_colorbar = False, norm = norm)

    ax.set_title(r'K (CFD)')
    ax.set_xlabel("x position [mm]")


    ax2 = fig.add_subplot(spec[0, 1])
    plot_obj = ds_cfd[ppu.CFD_KOH_SPECIES_NAME].plot(ax=ax2, x='pos_x', add_colorbar = False, norm = norm)

    ax_cb = fig.add_subplot(spec[1, :])
    cbar = plt.colorbar(plot_obj,  cax = ax_cb, location = "bottom")

    #Scale ticks to make units #*10^16/cm^3
    tick_scale_factor = 1E16
    cbar_ticks = cbar.get_ticks()
    new_ticks = cbar_ticks / tick_scale_factor
    cbar.set_ticklabels(new_ticks)
    ax_cb.set_xlabel(r"Concentration $[10^{16}/m^3]$")


else:     
    ax = fig.add_subplot(spec[0, 0])
    plot_obj = ds_cfd[ppu.CFD_K_SPECIES_NAME].plot(ax=ax, vmin=0, vmax = 1E16, x='pos_x', add_colorbar = False)

    cb = plt.colorbar(plot_obj)
    scale_cbar_ticks(cb, 10**(-16), new_label = r"K $\mathrm{[10^{16}/m^3]}$", rounding = 2)

    ax.set_title(r'K (CFD)')
    ax.set_xlabel("x position [mm]")


    ax2 = fig.add_subplot(spec[0, 1])
    plot_obj = ds_cfd[ppu.CFD_KOH_SPECIES_NAME].plot(ax=ax2, vmin=0, vmax = 2E16, x='pos_x', add_colorbar = False)

    cb2 = plt.colorbar(plot_obj)

    scale_cbar_ticks(cb2, 10**(-16), new_label = r"KOH $\mathrm{[10^{16}/m^3]}$", rounding = 2)

    ax2.set_title(r'KOH (CFD)')
    ax2.sharex(ax)
    ax2.set_xlabel("x position [mm]")






if shared_colorbar:
    ax3 = fig.add_subplot(spec[2,0])
else:
    ax3 = fig.add_subplot(spec[1,0])

da_plot = da_sel_tf.sel(gatedelay=gatedelays[0])
da_plot.plot(ax=ax3, vmin=0, cbar_kwargs = {"location":cbar_location, "label":"counts"})
ax3.set_title('ICCD before laser')

if shared_colorbar:
    ax4 = fig.add_subplot(spec[2,1])
else:
    ax4 = fig.add_subplot(spec[1,1])

da_plot = da_sel_tf.sel(gatedelay=gatedelays[1]) 
da_plot.plot(ax=ax4, vmin=0, cbar_kwargs = {"location": cbar_location, "label":"counts"})
ax4.set_title('ICCD during laser')

# draw beam slice box on ax4
ax4.add_patch(plt.Rectangle((160,-21), 30, 42, fill=False, edgecolor='r', lw=1, linestyle='--', alpha=0.5))

if shared_colorbar:
    ax5 = fig.add_subplot(spec[3,0])
else:
    ax5 = fig.add_subplot(spec[2,0])


ds_cfd['T'].plot(ax=ax5, x='pos_x', vmax=3200, cbar_kwargs = {"location": cbar_location, "label":"T [K]"})
ax5.set_title('Temperature')

#Must adjust layout like this, instead of with a kwarg in plt.figure, then place text. Otherwise the layout changes the text position
# fig.tight_layout()


# ax6 = 



if shared_colorbar:
    ax_lpf = fig.add_subplot(spec[3,1])
else:
    ax_lpf = fig.add_subplot(spec[2,1])

ax_lpf_twin = ax_lpf.twinx()

# axes[2, 1].cla()
# axes[2, 1].axis('off')
# plt.savefig(pjoin(DIR_FIG_OUT, 'spe_cfd_compare_2ns_2023-05-24.png'), bbox_inches='tight', dpi=300)
# set the 


da_sel2 = da_sel_tf.sel(x=beam_xslice, y=beam_yslice)

da_sel2.attrs['long_name'] = 'Counts'



ln_cam = da_sel2.mean(['x', 'y']).plot(marker='o', ax = ax_lpf)

# for ax in axes.flatten():

#     ax.set_anchor('W')

df_laser_profile_shift.plot(ax=ax_lpf_twin, color='r')
ln_laser = ax_lpf_twin.get_lines()

ax_lpf_twin.get_legend().remove()

lns = ln_cam + ln_laser
labs = ['ICCD during laser', 'Laser pulse time profile']
ax_lpf.legend(lns, labs, loc=0)

ax_lpf.set_xlabel('Gate Delay [ns]')
ax_lpf.set_ylabel('Intensity (normalized)')
ax_lpf.set_yticklabels(['', '0','','', '','', '', 1])
ax_lpf.tick_params(axis = 'y', which = 'both', right = True)
# ax_lpf_twin.set_ylabel('Laser Profile \n[counts]')

ax_lpf_twin.grid(False)
ax_lpf_twin.set_yticklabels([])
ax_lpf_twin.tick_params(axis = 'y', which = 'both', right = False)

ax_lpf.set_xlim(780,830)


#Set layout
fig.tight_layout(w_pad = 2, h_pad = 1)


#modify labels and add A, B, C... 
axes = [ax, ax2, ax3, ax4, ax5, ax_lpf]

right_axes = [ax2, ax4]



for ax, label in zip(axes, labels):
    X = ax.get_position().x0
    Y = ax.get_position().y1    
    fig.text(X - .07, Y + .015, label)
    # ax.tick_params(axis='y', which = "both", colors='white')
    
    ax.xaxis.label.set_color(label_color)
    ax.yaxis.label.set_color(label_color)

#set parameters for all axes except for the bottom right one
for ax in axes[:-1]:
    ax.set_ylabel('y position [mm]')
    ax.set_xlabel('x position [mm]')
    # ax.set_xticks([])
    # ax.set_yticks([])

    ax.set_xlim(0,250)
    ax.set_ylim(-25,25)

# for ax in right_axes:
#     ax.set_ylabel('')

#saving this as svg creates a huge file, I think because it's creating a vector object for every pixel of the CCD images
# output dir here will not be make yet in main data pipeline...
output_dir = os.path.join(REPO_DIR, 'final','figures', 'output')
if not os.path.exists(output_dir): os.makedirs(output_dir)
fig.savefig(os.path.join(output_dir, 'Fig2_ICCD_CFD.png'))

# %%
