#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-05-18'

data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)


dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()

fp = pjoin(data_folder, 'Lecroy', 'ds_silicon_power_init_12kV.cdf')

ds = xr.load_dataset(fp)

ds.coords['time'].attrs = dict(units='s', long_name='Time')   
ds.coords['time'] = ds.coords['time'].pint.quantify('s').pint.to('us')

# %%

from mhdpy.analysis.mws import MwsAccessor

ds = ds.mws.calc_mag_phase_AS()

# %%

ds[['i', 'q']].to_array('var').isel(acq_time=0).plot(hue='var')

#%%

ds['q'].isel(acq_time=[0,1,2]).plot(hue='acq_time', marker='o')

plt.xlim(-1,30)


#%%

ds['phase'].isel(acq_time=0).plot()
# %%

tw_power_sweep = slice(Timestamp('2023-05-18 19:13:17.603790336'), Timestamp('2023-05-18 19:30:55.000467200'), None)

# dsst['lasen_meter1']['Power'].plot()

dsst['filterwheel']['Filter Position'].sel(time=tw_power_sweep).plot(marker='o')

ds['i'].mean('time').plot(marker='o')

#%%

from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

from pi_paper_utils import convert_fw_pos_relpower

da_fw = dsst['filterwheel']['Filter Position'].sel(time=tw_power_sweep)

da_fw = convert_fw_pos_relpower(da_fw)

da_fw = da_fw.round(3)

da_fw
#%%

ds

#%%

ds = assign_signal(ds, da_fw, timeindex='acq_time')

ds = unstack_multindexed_acq_dim(ds, acq_dim_name='acq_time')

ds = ds.mean('mnum')

ds = ds.rename({'Filter Position': 'fw'})

ds

#%%

fig, axes = plt.subplots(3, sharex=True, figsize=(5,5))

ds['i'].plot(hue='fw', ax=axes[0])
axes[0].set_ylabel('I (V)')

ds['q'].plot(hue='fw', ax=axes[1])
axes[1].set_ylabel('Q (V)')

ds['mag'].plot(hue='fw', ax=axes[2])
axes[2].set_ylabel('Mag (V)')

plt.xlim(-1, 200)
for i, ax in enumerate(axes):
    if i != 0:
        ax.get_legend().remove()
    ax.set_xlabel('')

axes[0].get_legend().set_title('Laser Power (relative)')

axes[2].set_xlabel('Time (us)')

plt.savefig(pjoin(DIR_FIG_OUT, 'silicon_power_raw.png'), bbox_inches='tight')

#%%

ds[['i','q']].mws.calc_mag_phase_AS()['AS'].plot(hue='fw')

plt.yscale('log')

plt.xlim(-1, 200)

plt.ylim(1e-5,)

#%%

from mhdpy.analysis.mws.fitting import pipe_fit_exp


# da_fit = da_sel.mean('mnum')

da_fit = ds.copy()['AS'].drop(0, 'fw')

da_fit

#%%

da_fit

#%%


pipe_fit_exp.perform_fit_kwargs.update({'fit_timewindow': slice(Quantity(60, 'us'),Quantity(130, 'us'))})
pipe_fit_exp.fit_prep_kwargs['pre_norm_cutoff'] = None
ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(da_fit)

#%%
# %%

ds_mws_fit.to_array('var').plot(hue='var', row='fw')

plt.yscale('log')

#%%

ds_p['decay'].plot(marker='o')


#%%

da_fit_sel = da_fit.sel(fw=1)

ds_mws_fit, ds_p, ds_p_stderr = pipe_fit_exp(
    da_fit_sel
    )

#%%

plt.figure(figsize=(4,2))

ds_mws_fit[['AS_all','AS_fit']].to_array('var').plot(hue='var')

plt.yscale('log')

plt.xlim(-1, 200)

plt.ylim(1e-4,)

decay_str = f'$\\tau = {ds_p.decay.values:.2f} \\mu s$'
# plt.text(100, 1, decay_str, fontsize=12)
# text centered at top of plot

plt.text(130, 0.2, decay_str, fontsize=14, ha='center', va='bottom')

plt.title('')

plt.legend(['Data','Fit'], loc='lower left')

plt.ylabel('AS')

plt.savefig(pjoin(DIR_FIG_OUT, 'silicon_power_fit.png'), bbox_inches='tight')

# %%
