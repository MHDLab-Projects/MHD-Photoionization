#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu
from mhdlab.plot import xr_errorbar_axes

# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300

data_directory = pjoin(REPO_DIR, 'final', 'dataset', 'output')


ds_lecroy = xr.open_dataset(pjoin(data_directory, 'ds_pos_lecroy.cdf')).xr_utils.stack_run()
ds_lecroy = ds_lecroy.sel(phi=0.8, method='nearest')

ds_lecroy_536 = ppu.fileio.load_lecroy('536_pos', avg_mnum=False, AS_calc='absolute')
ds_lecroy_536 = ds_lecroy_536.mwt.calc_time_stats()

# Position dependent barrel AES
ds_p = xr.open_dataset(pjoin(data_directory, 'ds_p_stats_pos.cdf')).xr_utils.stack_run()
ds_p = ds_p.pint.quantify()
ds_p['nK_mw_horns'] = ds_p['nK_mw_horns'].pint.to('particle/cm^3')

# Position dependent AES spectra
ds_alpha_fit = xr.open_dataset(pjoin(data_directory, 'ds_pos_alpha_fit.cdf')).xr_utils.stack_run()

# Non position dependent barrel AES
ds_p_barrel = xr.open_dataset(pjoin(data_directory, 'ds_p_barrel.cdf')).xr_utils.stack_run()
ds_p_barrel = ds_p_barrel.pint.quantify()
ds_p_barrel['nK_barrel_cfdprofile'] = ds_p_barrel['nK_barrel_cfdprofile'].pint.to('particle/cm^3')

#%%

# This shows the dates of barrel exit data for 0.8 phi
ds_p_barrel.sel(phi=0.8, method='nearest').dropna('run', how='all').coords['run']

#%%

ds_cfd = ppu.fileio.load_cfd_centerline()

ds_cfd = ds_cfd.sel(kwt=1)

ds_cfd['nK_cm3'] = ds_cfd[ppu.CFD_K_SPECIES_NAME].pint.to('particle/cm^3')


#%%

nrows = 3
fig = plt.figure()

fig.set_figheight(nrows * 3)
fig.set_figwidth(5)

# fig.set_figheight(8)

sfigs = fig.subfigures(2,1, height_ratios = [1,2], hspace = 0)

sfigs[0].subplots(ncols = 3, sharey = True)


ds_plot = ds_alpha_fit.sel(run=('2023-05-18', 1)).sel(phi=0.8, method='nearest')

# da_plot.plot(col='phi', hue='var', row='motor')
motor_sel = [80, 130, 180]

ds_plot = ds_plot.sel(motor=motor_sel, method='nearest')

# fig, axes = plt.subplots(1,len(motor_sel), figsize=(5,2), sharey=True)

axes = sfigs[0].get_axes()

for i, motor in enumerate(ds_plot.coords['motor']):
    ax = axes[i]
    ds_plot.sel(motor=motor).to_array('var').plot(hue='var', ax=ax)
    ax.set_title('{:.0f} mm'.format(motor))

    ax.set_xlabel("Wavelength [nm]")

    ax.set_ylim(-0.1,1.1)
    ax.set_xlim(760,775)

    if i != 2:
        ax.get_legend().remove()

    else:
        ax.legend(['Fit', 'Fitted Data', 'Raw Data'])

axes[0].set_ylabel('$\\alpha$')



(ax1, ax2) = sfigs[1].subplots(nrows = 2)

# Phi =0.8
# ax1 = axes[0]

ds_p_sel = ds_p.sel(phi=0.8, method='nearest').pint.dequantify()
ds_cfd_sel = ds_cfd.sel(phi=0.8, method='nearest').pint.dequantify()

phi_val_expt = ds_p_sel.coords['phi'].item()


nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: removing last datapoint. Can this be fixed with calibraiton/processing?
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,220))

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-18 Run 1', ax=ax1)
g = nK_mw_horns.isel(run=1).plot(x='motor', marker='o', label='2023-05-18 Run 2', ax=ax1)

ds_cfd_sel['nK_cm3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax1)
# ds_cfd_sel['nK_cm3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax1)


# Barrel exit
nK_barrel_mean = ds_p_barrel['nK_barrel_cfdprofile'].mean('run').sel(phi=0.8, method='nearest')
nK_barrel_std = ds_p_barrel['nK_barrel_cfdprofile'].std('run').sel(phi=0.8, method='nearest')
ax1.errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='o', label='Barrel Exit', capsize=5, )

ax1.set_title('Equivalence Ratio = {}'.format(phi_val_expt))


# Phi = 0.65

# ax2 = axes[1]

ds_p_sel = ds_p.sel(phi=0.6, method='nearest').pint.dequantify()
ds_cfd_sel = ds_cfd.sel(phi=0.6, method='nearest').pint.dequantify()

phi_val_expt = ds_p_sel.coords['phi'].item()


nK_mw_horns = ds_p_sel['nK_mw_horns'].dropna('run', how='all').dropna('motor', how='all')

#TODO: final datapoint for 516 has no potassium peaks. However, calibration is off and fit gives a large absorption value. Removing is probably best as any absorption is not meaninful anyway. 
nK_mw_horns = nK_mw_horns.sel(motor = slice(50,200)) 

g = nK_mw_horns.isel(run=0).plot(x='motor', marker='o', label='2023-05-24 Run 1', ax=ax2)

ds_cfd_sel['nK_cm3'].sel(offset=0).plot(color='black', label ='CFD centerline', linestyle='-.', ax=ax2)
# ds_cfd_sel['nK_cm3'].sel(offset=3).plot(color='black', label ='CFD 3mm offset', linestyle='--', ax=ax2)

nK_barrel_mean = ds_p_barrel['nK_barrel_cfdprofile'].mean('run').sel(phi=0.6, method='nearest')
nK_barrel_std = ds_p_barrel['nK_barrel_cfdprofile'].std('run').sel(phi=0.6, method='nearest')

assert nK_barrel_std.item().magnitude ==0 # there is only one run for phi=0.6

# ax2.errorbar(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, yerr=nK_barrel_std, color='green', marker='o', label='Barrel Exit', capsize=5, )
ax2.scatter(ppu.AES_BARREL_OFFSET.to('mm').magnitude, nK_barrel_mean, color='green', marker='o', label='Barrel Exit')


ax2.set_title('Equivalence Ratio = {}'.format(phi_val_expt))

for ax in [ax1, ax2]:
    ax.axvline(180, color='gold', linestyle='--')
    ax.set_ylim(1e10,1e17)
    ax.set_xlim(-10,220)
    ax.set_yscale('log')
    ax.legend()
    ax.set_ylabel('$n_K [\\mathrm{cm^{-3}}$]')

ax2.set_xlabel('Stage Position [mm]')
ax1.set_xlabel('')
ax1.set_xticklabels([])


#bounding box position changes for subfigure with multiple subplots, so the top subfigure's text does not end up in the same place as the bottom two. Right now, doing this hacky way of getting the top label to the right place. More robust would be to use figure coordinates for putting text, but that makes it way harder to put labels on multiple subplots.

labels = ['B)', 'C)']

axlist = [ax1, ax2]

for ax, label in zip(axlist, labels):
    ax.text(-0.15, 1.1, label, transform=ax.transAxes, va='bottom', ha='left')

sfigs[0].get_axes()[0].text(-0.15, 2.75, 'A)', transform=ax1.transAxes)

# plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mwt_cfd.png'), dpi=300, bbox_inches='tight')

# # Phi =0.65
# # ax2 = axes[1]

# #%%
# plt.savefig(pjoin(DIR_FIG_OUT, 'pos_alpha_fit.png'), dpi=300, bbox_inches='tight')

#%%

#Concatenate three figures into a single figure

ds_plot = ds_lecroy_536.sel(run=('2023-05-18', 1)).dropna('motor', how='all')

da_plot = ds_plot['AS_abs']

motor_sel = [30, 100, 180]
motor_actual = []

nrows = 3 #to scale plot height by number of rows

fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3)
fig.set_figheight(nrows * 3)

fig.subplots_adjust(hspace = 0.25)


labels = ['A)', 'B)', 'C)']

axlist = [ax1, ax2, ax3]

for ax, label in zip(axlist, labels):
    ax.set_title('')
    ax.text(-0.15, 1.1, label, transform=ax.transAxes, va='top', ha='right')


#Create two subfigures, where bottom subfigure has no vertical space between subplots


# fig = plt.figure()
# fig.set_figheight(nrows * 3)

# sfigs = fig.subfigures(2, 1, height_ratios = [1,2], hspace = 0)
# ax1 = sfigs[0].subplots(1, 1)
# (ax2, ax3) = sfigs[1].subplots(2, 1, sharex = True)

# sfigs[1].subplots_adjust(hspace = 0)



colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Define a list of colors for the traces

for i, motor in enumerate(motor_sel):
    # Calculate mean and standard deviation
    mean = da_plot.sel(motor=motor,method='nearest').mean(dim='mnum', keep_attrs=True)
    std = da_plot.sel(motor=motor,method='nearest').std(dim='mnum', keep_attrs=True)

    motor_val = da_plot.sel(motor=motor, method='nearest').coords['motor'].values
    motor_actual.append(motor_val)
    
    # Plot mean and shaded region for standard deviation
    ax1.plot(mean['time'].pint.magnitude, mean.pint.magnitude, color=colors[i % len(colors)], label='{:.0f} mm'.format(motor_val))
    ax1.fill_between(mean['time'].pint.magnitude, (mean-std).pint.magnitude, (mean+std).pint.magnitude, color=colors[i % len(colors)], alpha=0.2)

ax1.set_xlabel('Time [$\mu s$]')
ax1.set_ylabel('AS')
ax1.legend(bbox_to_anchor=(0.5,0.87))
ax1.set_ylim(-0.05,1.05)



# plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_time_20230518.png'), dpi=300, bbox_inches='tight')


da_plot = ds_lecroy_536['AS_abs'].sel(run=('2023-05-18', 1)).dropna('motor', how='all')

AS_pp = da_plot.sel(time=slice(-50,-1)).mean('mnum').mean('time').pint.dequantify()
AS_pp_fluct = da_plot.sel(time=slice(-50,-1)).std('mnum').mean('time').pint.dequantify()
AS_laser = da_plot.mwt._pulse_max().mean('mnum').pint.dequantify()
AS_laser_fluct = da_plot.mwt._pulse_max().std('mnum').pint.dequantify()

xr_errorbar_axes(AS_pp, AS_pp_fluct, axes=ax2, capsize=5, label='AS Pre Laser')
xr_errorbar_axes(AS_laser, AS_laser_fluct, axes=ax2, capsize=5, label='AS Laser')

for motor_val in motor_actual:
    plt.axvline(motor_val, color='k', linestyle='--', alpha=0.5)

dropna(ax2)


ax2.legend()
ax2.set_xlim(-5,270)
ax2.set_xlabel('Stage Position [mm]')
ax2.set_ylabel('AS')

# plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_stats_20230518.png'))


da_plot = ds_lecroy_536['SFR_abs'].dropna('run', how='all')

g =  da_plot.mean('mnum', keep_attrs=True).plot(hue='run_plot', x='motor', marker='o', ax = ax3)

for motor_val in motor_actual:
    ax3.axvline(motor_val, color='k', linestyle='--', alpha=0.5)



ax3.set_title('')
ax3.set_xlim(-5, 270)
ax3.set_ylim(-5,125)
ax3.set_xlabel('Stage Position [mm]')
ax3.set_ylabel('SFR')

dropna(g)

fig.savefig('Fig4_Position_MWT.svg')




# %%

axlist = [ax1, ax2, ax3]

fig_combined, axs = plt.subplots()

[fig_combined.add_subplot(ax) for ax in axlist]

fig_combined.show()
# %%
