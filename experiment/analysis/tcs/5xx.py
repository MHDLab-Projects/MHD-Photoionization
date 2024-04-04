#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu

from mhdpy.analysis import mws
from mhdpy.analysis import absem

from mhdpy.plot import dropna

#%%
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dst_coords = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()

#%%

ds_cfd = ppu.fileio.load_cfd_centerline()
ds_cfd = ds_cfd.sel(offset=0)

# %%

tc = '5x6_pos'

ds_absem = ppu.fileio.load_absem(tc)

# ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='relative')
ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')
da_lecroy = ds_lecroy['AS_abs']

#%%

da_max = da_lecroy.mws._pulse_max()

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')

#%%

da_max.plot(hue='phi', marker='o')

plt.twinx()

da_cfd_sel = ds_cfd.sel(kwt=1)['Yeq_KOH']

da_cfd_sel.plot(hue='phi', linestyle='--')

plt.gca().get_legend().set_bbox_to_anchor((1, 0.6))


#%%


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

da_cfd_sel = ds_cfd.sel(kwt=1)['Yeq_KOH']

lines = []
labels = []

for phi in da_max['phi']:
    # line, = ax1.plot(da_max.sel(phi=phi), label=phi.values, marker='o')
    line, = da_max.sel(phi=phi).plot(label=phi.values, marker='o', ax=ax1)
    lines.append(line)
    labels.append(phi.values)

ax1.legend(lines, labels, title="Expt.", loc="upper right")
ax1.set_ylim(-0.05,0.3)

lines = []
labels = []

for phi in da_cfd_sel['phi']:
    line, = ax2.plot(da_cfd_sel.sel(phi=phi), label=phi.values, linestyle='--')
    lines.append(line)
    labels.append(phi.values)

ax2.legend(lines, labels, title="CFD", loc="upper right")
ax2.set_ylim(0.4e16,)

ax1.set_xlim(0,260)

plt.savefig(pjoin(DIR_FIG_OUT, '5x6_pos_mws_KOH.png'), dpi=300)


#%%


da_lecroy_mean = da_lecroy

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)

plt.xlim(-1,)

#%%

da_lecroy_mean = da_lecroy

da_lecroy_mean.plot(col='motor',hue='phi', col_wrap=3)

plt.yscale('log')

plt.xlim(-1,)

plt.savefig(pjoin(DIR_FIG_OUT, '5x6_pos_mws_AStime_phi.png'), dpi=300)


#%%

da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']
# da = da.max('wavelength').count('mnum')
da = da.max('wavelength').mean('mnum')

g = da.plot(hue='phi', marker='o')

dropna(g)

#%%
da = ds_absem.sel(mp='mw_horns').drop('run')['alpha']

da.sel(motor=100, method='nearest', phi=1).dropna('mnum', how='all')


#%%

## Examine time data...
# I believe this is redundant with other time analysis, probably can remove 

ds_absem_time = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'ds_absem.cdf'))

ds_absem_time['diff'] = ds_absem_time['led_on'] - ds_absem_time['led_off']
ds_absem_time['alpha'] = 1 - ds_absem_time['diff']/ds_absem_time['calib']

ds_absem_time = ds_absem_time.set_index(acq=['time','mp']).unstack('acq')

ds_absem_time = ds_absem_time['alpha']

ds_absem_time

#%%
dsst = mhdpy.fileio.TFxr(pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')).as_dsst()

# tw = slice(Timestamp('2023-05-24 19:45:01.091800832'), Timestamp('2023-05-24 20:39:19.309871616'), None)
tw = slice(Timestamp('2023-05-24 20:12:07.301042944'), Timestamp('2023-05-24 20:12:46.490260736'), None)

da = ds_absem_time.sel(time=tw).dropna('time', how='all')
da.mean('wavelength').plot(hue='mp', marker='o')

plt.twinx()

dsst['hvof']['CC_equivalenceRatio'].sel(time=tw).plot()

#%%

# Check assigned coord values

da = ds_absem_time.sel(time=tw).dropna('time', how='all')
da.mean('wavelength').plot(hue='mp', marker='o')

plt.twinx()

dst_coords['phi'].sel(time=tw).plot(marker='o')

#%%


from mhdpy.analysis.absem.fitting import pipe_fit_alpha_1

spectral_reduction_params_fp = os.path.join(REPO_DIR,'experiment','metadata', 'spectral_reduction_params.csv')
spect_red_dict = pd.read_csv(spectral_reduction_params_fp, index_col=0).squeeze().to_dict()

ds_fit = ds_absem.mean('mnum')


ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_1(ds_fit, fit_prep_kwargs={'spect_red_dict': spect_red_dict})
ds_alpha_fit['alpha'] = ds_fit['alpha']


#%%

ds = ds_alpha_fit[['alpha', 'alpha_fit']]
ds = ds.rename({'alpha': 'data', 'alpha_fit':'fit'})
ds = ds.to_array('var')


# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)
# %%
# %%

tc = '5x3_pos'

ds_absem = ppu.fileio.load_absem(tc)

ds_lecroy = ppu.fileio.load_lecroy(tc, avg_mnum=True, AS_calc='absolute')
da_lecroy = ds_lecroy['AS_abs']

#%%


da_lecroy_mean = da_lecroy

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)

plt.xlim(-1,)

#%%

da_lecroy_mean = da_lecroy

da_lecroy_mean.plot(row='motor',hue='phi')

plt.yscale('log')


plt.xlim(-1,)

#%%[markdown]

# # quick partial fix calibraiton data for 5x3_pos. 


# the calibration data appears to have a significant offset on the high wavlength side.
# I looked at the alpha in 53x notebook and it doesn't appear to be as bad as here
# TODO: some method of correcting calibration for this case. beta_offfset_line? or multiply calibration by a linear function of wavelength?

#%%

fp_calib = pjoin(REPO_DIR, 'experiment', 'data','proc_data', 'ds_calib.cdf')

ds_calib = xr.load_dataset(fp_calib)

ds_calib

#%%

da_calib = ds_calib.sel(mp='mw_horns').sel(time=Timestamp('2023-05-24T18:12:49.047832386'))
# da_calib = ds_calib.sel(mp='mw_horns').sel(time=Timestamp('2023-05-24T22:42:25.960157725'))

da_calib

ds_absem['calib'] = da_calib['diff']

ds_absem = ds_absem.absem.calc_alpha()

#%%

ds_sel = ds_absem.mean('run').sel(mp='mw_horns').mean('mnum')

# ds_sel[['diff','calib']].to_array('var').plot(col='motor', row='phi', hue='var')

#%%

ds_sel[['diff', 'calib']].to_array('var').sel(motor=180, method='nearest').plot(hue='var', row='phi')


#%%

ratio = (ds_sel['diff']/ds_sel['calib'])

ratio.sel(motor=50, method='nearest').plot(row='phi')

plt.ylim(0.8,1.2)


#%%
ds_sel['alpha'].sel(motor=50, method='nearest').plot(hue='phi')

plt.ylim(0,1)

plt.savefig(pjoin(DIR_FIG_OUT, '5x3_pos_alpha.png'), dpi=300)


# %%

da_lecroy_mean = da_lecroy

da_lecroy_mean.plot(row='motor',hue='phi', sharey=False)



#%%

da_max = da_lecroy_mean.mws._pulse_max()

da_max.plot()
#%%

da_max.plot(hue='phi', marker='o')

#%%

da_max

#%%


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), sharex=True)

da_cfd_sel = ds_cfd.sel(kwt=0.1)['Yeq_KOH']

lines = []
labels = []

for phi in da_max['phi']:
    # line, = ax1.plot(da_max.sel(phi=phi), label=phi.values, marker='o')
    line, = da_max.sel(phi=phi).plot(label=phi.values, marker='o', ax=ax1)
    lines.append(line)
    labels.append(phi.values)

ax1.legend(lines, labels, title="Expt.", loc="upper right")
ax1.set_ylim(-0.01,0.1)

lines = []
labels = []

for phi in da_cfd_sel['phi']:
    line, = ax2.plot(da_cfd_sel.sel(phi=phi), label=phi.values, linestyle='--')
    lines.append(line)
    labels.append(phi.values)

ax2.legend(lines, labels, title="CFD", loc="lower right")
ax2.set_ylabel("CFD KOH [#/cm^3]")

ax1.set_xlim(0,260)
ax2.set_ylim(0.4e15)

plt.savefig(pjoin(DIR_FIG_OUT, '5x3_pos_mws_KOH.png'), dpi=300)



#%%

from mhdpy.analysis.absem.fitting import pipe_fit_alpha_2

ds_fit = ds_absem.isel(run=0) #only one run
ds_fit = ds_fit.sel(wavelength=slice(750,790)).mean('mnum')

ds_alpha_fit, ds_p, ds_p_stderr = pipe_fit_alpha_2(ds_fit)

#%%

ds_alpha_fit['alpha_red'].sel(mp='mw_horns').plot(hue='phi', row='motor')


# %%

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

g = da.plot(hue='phi', marker='o')


dropna(g)

plt.yscale('log')
# %%

da_count = ds_absem['alpha'].isel(wavelength=0).count('mnum')

da_count.sel(mp='mw_horns').plot(hue='phi', marker='o')

# %%



da_cfd_sel = ds_cfd.sel(kwt=0.1)['Yeq_K'].pint.to('particle/m^3')

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')


for phi in da['phi']:
    da.sel(phi=phi).plot(label=phi.values, marker='o')


for phi in da_cfd_sel['phi']:
    da_cfd_sel.sel(phi=phi).plot(label=phi.values, linestyle='--')


plt.legend()
plt.ylim(1e17,)
plt.xlim(0,300)

plt.yscale('log')
#%%

plt.figure(figsize=(6,4))

import matplotlib.lines as mlines

da_cfd_sel = ds_cfd.sel(kwt=0.1)['Yeq_K'].pint.to('particle/m^3')

da = ds_p['nK_m3'].sel(mp='mw_horns').drop('run')

lines = []
labels = []

for phi in da['phi']:
    line, = da.sel(phi=phi).plot(label=phi.values, marker='o')
    lines.append(line)
    labels.append(phi.values)

legend1 = plt.legend(lines, labels, title="Expt.", loc="lower right",)

lines = []
labels = []

for phi in da_cfd_sel['phi']:
    line, = da_cfd_sel.sel(phi=phi).plot(label=phi.values, linestyle='--')
    lines.append(line)
    labels.append(phi.values)

legend2 = plt.legend(lines, labels, title="CFD", loc="lower left")

plt.gca().add_artist(legend1)

plt.ylim(1e17,)
plt.xlim(0,360)

plt.title('')

plt.yscale('log')

plt.savefig(pjoin(DIR_FIG_OUT, '5x3_pos_nK_m3.png'), dpi=300)
# %%
