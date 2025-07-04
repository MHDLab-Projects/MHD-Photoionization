#%%


from mhdlab.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
import pi_paper_utils as ppu

# plt.rcParams.update({'font.size': 16})

# %%

tc = '536_pos'

ds_aas = ppu.fileio.load_aas(tc)
ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')

# add 536 photodiode data from 5x6_pos run. Do not have photodiode data for 536_pos run
tc = '5x6_pos'

ds_lecroy_5x6 = ppu.fileio.load_lecroy('5x6_pos', avg_mnum=True, AS_calc='relative')
ds_lecroy_5x6 = ds_lecroy_5x6.sel(phi=0.8, method='nearest')
ds_lecroy_5x6 = ds_lecroy_5x6[['pd1','pd2']]

ds_lecroy = xr.merge([ds_lecroy_5x6, ds_lecroy])

ds_lecroy = ds_lecroy.mwt.calc_time_stats()


# ds_lecroy.to_array('var').mean('mnum').mean('motor').mean('run').sel(time=slice(-1,1)).plot(col='var', sharey=False)

#%%

ds_lecroy['pd1'].isnull().all()

#%%[markdown]

# # aas

#%%

ds_aas.sel(mp='barrel')['alpha'].mean('mnum').mean('wavelength').dropna('run', how='all')

#%%

ds_aas['alpha'].mean('mnum').plot(col='mp', hue='run_plot',row='motor', x='wavelength')
plt.ylim(0,1.1)
plt.xlim(760,780)


#%%

from mhdlab.analysis.aas.fitting import pipe_fit_alpha_1, pipe_fit_alpha_2

ds_fit = ds_aas.mean('mnum')
ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

#%%

ds = ds_alpha_fit[['alpha_red', 'alpha_fit']]
ds = ds.rename({'alpha_red': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')

ds_p.coords['motor'].attrs = dict(long_name="Stage Position", units='mm')
ds_p.coords['run_plot'].attrs = dict(long_name="Run")

#%%

g  = ds_p['nK_m3'].plot(hue='run_plot', x='motor',col='mp', marker='o')

dropna(g)

plt.yscale('log')

#%%

plt.figure(figsize=(6,3))

da_sel = ds_p['nK_m3'].sel(mp='mw_horns')

da_sel = da_sel.dropna('run', 'all')

da_sel.plot(hue='run_plot', x='motor', marker='o')

plt.yscale('log')

#%%[markdown]

# # Lecroy

# %%

g = ds_lecroy['AS_abs'].plot(hue='run_plot', row='motor', x='time')

dropna(g)

#%%

motor_sel = [34.81, 104.8, 178, 226.7]

da_sel = ds_lecroy['AS_abs'].sel(motor=motor_sel, method='nearest')

fig, axes = plt.subplots(4, figsize=(3,12), sharex=True, sharey=True)

for i, motor in enumerate(da_sel.coords['motor'].values):
    ax = axes[i]
    da_sel.sel(motor=motor).plot(hue='run_plot', x='time', ax=ax)
    ax.set_title('Position: {} mm'.format(motor))
    ax.set_xlabel('')
    ax.set_ylabel('AS_abs')

    if i == len(da_sel.coords['motor'].values) - 1:
        ax.set_xlabel('Time [us]')

    else:
        ax.get_legend().remove()

# da_sel.plot(row='motor', hue='run_plot', x='time', figsize=(8,25))

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_AS_sel.png'), dpi=300, bbox_inches='tight')
#%%

da_sel = ds_lecroy['AS_abs'].sel(motor=motor_sel, method='nearest')

da_sel = da_sel.mean('run')

da_sel.plot(hue='motor', x='time', marker='o')

plt.yscale('log')

plt.xlim(-1, 20)


#%%

# plt.figure(figsize=(4,6))

# da_sel = ds_lecroy['AS_abs'].mean('mnum').sel(motor=[50,100,150,180,225], method='nearest')
da_sel = ds_lecroy['mag'].sel(motor=[50,100,150,180,225], method='nearest')

da_sel = da_sel.dropna('run', how='all')

da_sel.plot(col='motor', hue='run_plot', x='time', figsize=(10,3))

# plt.yscale('log')
#%%
ds_lecroy

#%%

nK = ds_p['nK_m3'].sel(mp='mw_horns')


ds_stats_lecroy = ds_lecroy.mwt.calc_time_stats()[['dAS_abs_max', 'mag_pp', 'mag_fluct', 'SFR_abs', 'dpd1_max']]

#%%


#%%

ds = xr.merge([
    nK,
    ds_stats_lecroy,
    ]).sortby('motor').dropna('run', how='all')

ds
#%%

g = ds.pint.dequantify().to_array('var').plot(hue='run_plot', row='var', x='motor', sharey=False)

dropna(g)

#%%


fig, axes = plt.subplots(5, figsize=(4,12), sharex=True)

ds['nK_m3'].plot(hue='run_plot', x='motor', ax=axes[0],marker='o')
ds['dAS_abs_max'].plot(hue='run_plot', x='motor', ax=axes[1], marker='o')
ds['mag_fluct'].plot(hue='run_plot', x='motor', ax=axes[2],marker='o')
ds['SFR_abs'].plot(hue='run_plot', x='motor', ax=axes[3],marker='o')
ds['dpd1_max'].plot(hue='run_plot', x='motor', ax=axes[4],marker='o')

for ax in axes:
    ax.get_legend().remove()
    ax.set_title('')
    ax.set_xlabel('')


axes[3].set_xlabel('Position [mm]')

#%%


fig, axes = plt.subplots(3, figsize=(3,10), sharex=True)

for i, var in enumerate(['dAS_abs_max', 'mag_fluct', 'SFR_abs']):
    ax = axes[i]

    ds_stat = ds.wma.initialize_stat_dataset(var, 'run')

    # plot with confidence interval

    # ds_stat['mean'].plot(ax=ax)
    # ax.fill_between(ds_stat.coords['motor'], ds_stat['mean'] - ds_stat['std'], ds_stat['mean'] + ds_stat['std'], alpha=0.2)

    ax.errorbar(ds_stat.coords['motor'], ds_stat['mean'], yerr=ds_stat['std'], marker='o', capsize=5)

    ax.set_title('')
    ax.set_xlabel('')

    unit_str = '[' + ds[var].attrs['units']+ ']' if 'units' in ds[var].attrs else '' 
    ax.set_ylabel(ds[var].attrs['long_name'] + unit_str)
    
axes[2].set_xlabel('Position [mm]')

ta = axes[0].twinx()

delta_pd1 = ds['dpd1_max'].mean('run')
# delta_pd1 = delta_pd1.pint.quantify('V').pint.to('mV')
delta_pd1.plot(ax=ta, color='red', marker='o')
ta.set_ylabel('Delta PD1 [mV]')
ta.set_title('')

for mp in motor_sel:
    if mp == 178:
        axes[2].axvline(mp, color='gold', linestyle='--')
    else:
        axes[2].axvline(mp, color='gray', linestyle='--')


plt.savefig(pjoin(DIR_FIG_OUT, 'pos_mwt_stats.png'), dpi=300, bbox_inches='tight')


#%%

ds.unstack('run').plot.scatter(x='nK_m3', y='mag_pp')

plt.xscale('log')
# plt.yscale('log')

# %%

from mhdlab.pyvista_utils import CFDDatasetAccessor

fp_cfd = pjoin(os.getenv('REPO_DIR'), 'final', 'dataset', 'output', 'cfd_profiles_centerline.cdf' )

ds_cfd = xr.load_dataset(fp_cfd)
ds_cfd = ds_cfd.cfd.quantify_default()
ds_cfd = ds_cfd.cfd.convert_all_rho_number()

ds_cfd = ds_cfd.sel(kwt=1).sel(phi=0.8).sel(offset=0)

ds_cfd['nK_m3'] = ds_cfd[ppu.CFD_K_SPECIES_NAME].pint.to('particle/m^3')

ds_cfd_norm = ds_cfd/ds_cfd.max()



#%%

ds_cfd_beam = ppu.fileio.load_cfd_beam()

ds_cfd_beam = ds_cfd_beam.sel(offset=0).sel(kwt=1)

xs = ds_cfd_beam.coords['dist'].values
xs = (xs - xs[0])*1000

ds_cfd_beam = ds_cfd_beam.assign_coords(dist = xs)

#TODO: Convert to number density
ds_cfd_beam_norm = ds_cfd_beam/ds_cfd_beam.max()

#%%

ds_cfd_beam

#%%

plt.figure(figsize=(5,3))

da = ds['nK_m3'].dropna('run', how='all')
da = da.sel(motor = slice(50,300))

g = da.plot(hue='run_plot', x='motor', marker='o')
plot.dropna(g)

ax1 = plt.gca()
ax1.get_legend().set_title('Experiment (date, #)')
ax1.get_legend().set_bbox_to_anchor([0,0,0.5,0.4])
ax1.set_yscale('log')
ax1.set_title('')
ax1.set_xlabel('Position [mm]')

ax1.axvline(178, color='gold', linestyle='--')

plt.twinx()

ds_cfd_norm[ppu.CFD_K_SPECIES_NAME].plot(color='black', label ='CFD centerline')
# ds_cfd_beam_norm['K'].plot(color='grey', label = 'beam conv (TODO)', marker='o')

ax2 = plt.gca()
ax2.legend()
ax2.set_ylabel('Normalized $n_{K, CFD}$')
ax2.set_yscale('log')
ax2.set_title('')

plt.xlim(0,310)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mwt_cfd.png'), dpi=300, bbox_inches='tight')

#%%

nK_barrel_mean = ds_p['nK_m3'].sel(mp='barrel').mean('motor').mean('run')
nK_barrel_std = ds_p['nK_m3'].sel(mp='barrel').mean('motor').std('run')


#%%

plt.figure(figsize=(4.5,2))

da = ds['nK_m3'].dropna('run', how='all')
da = da.sel(motor = slice(50,300))

g = da.isel(run=0).plot(x='motor', marker='o', label='2023-05-12 Run 1')
dropna(g)
g = da.isel(run=1).plot(x='motor', marker='o', label='2023-05-12 Run 2')
dropna(g)

ax1 = plt.gca()
ax1.set_yscale('log')
ax1.set_title('')
ax1.set_xlabel('Position [mm]')

# ax1.get_legend().set_title('Experiment (date, #)')
# ax1.get_legend().set_bbox_to_anchor([0,0,0.5,0.4])
ax1.axvline(178, color='gold', linestyle='--')


ds_cfd['nK_m3'].plot(color='black', label ='CFD centerline')

plt.errorbar(0, nK_barrel_mean, yerr=nK_barrel_std, color='red', marker='o', label='Barrel nK avg')

plt.ylim(1e16,3e22)
plt.xlim(-10,250)



ax1.legend()
ax1.set_title('')

ax1.set_xlabel('Position [mm]')

# ax2 = plt.gca()
# ax2.legend()
# ax2.set_ylabel('Normalized $n_{K, CFD}$')
# ax2.set_yscale('log')
# ax2.set_title('')

# plt.xlim(0,310)

plt.savefig(pjoin(DIR_FIG_OUT, 'pos_nK_mwt_cfd.png'), dpi=300, bbox_inches='tight')




# %%

g = ds['dAS_abs_max'].plot(hue='run_plot', x='motor', marker='o')

dropna(g)

leg = plt.gca().get_legend()
leg.set_title('Experiment (date, #)')
leg.set_bbox_to_anchor([0,0,1.5,1])

plt.title('CFD KOH and MWT AS compare')

plt.twinx()

ds_cfd_norm[ppu.CFD_KOH_SPECIES_NAME].plot(color='black', label ='CFD centerline')
# ds_cfd_beam_norm[ppu.CFD_KOH_SPECIES_NAME].plot(color='gray', label = 'beam conv (TODO)')

ax2 = plt.gca()
ax2.legend(bbox_to_anchor=[0,0,1.5,0.3])

plt.xlim(0,310)
plt.title('')

# %%
