"""
Analysis of 2D parameter sweeps

#TODO: Add other dates/2d sweeps, improve plots, and convert to use main function like analsis_1d.py
"""
#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_LECROY_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

figure_out_dir = pjoin(DIR_DATA_OUT, 'mwt_2d', datestr)

dsst = mhdlab.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)
# #%%

tc_dict = {
'2023-04-07':[
    '536_pos_power'
],
'2023-05-24':[
    '5x6_pos',
    '5x3_pos'
],
}

tcs = tc_dict[datestr]

hudim_dict = {'2023-04-07':'power', '2023-05-24':'phi'}

huedim = hudim_dict[datestr]


coldim = 'motor'

from mhdlab.analysis import mwt

mag_0 = ppu.fileio.load_mwt_T0()

dss = {}

for tc in tcs:

    fp_in = pjoin(DIR_LECROY_PROC_DATA, '{}.cdf'.format(tc))

    ds_in = xr.load_dataset(fp_in)

    ds_in = ds_in.mwt.calc_AS_abs(mag_0=mag_0).mwt.calc_time_stats()
    ds_in = ds_in.sel(date=datestr).sel(run_num=1)
    ds_in = ds_in[['i','q','mag','phase','AS_abs', 'dAS_abs']] #Plotted variables

    tc_dim = [dim for dim in ds_in.dims if dim not in ['time','mnum']][0]
    mnum_counts = ds_in['i'].mean('time').groupby(tc_dim).count('mnum')

    ds = ds_in.mean('mnum', keep_attrs=True)

    tc_fig_dir = pjoin(figure_out_dir, tc)
    if not os.path.exists(tc_fig_dir): os.makedirs(tc_fig_dir)

    for fn in os.listdir(tc_fig_dir):
        os.remove(os.path.join(tc_fig_dir, fn))

    ds['count'] = mnum_counts
    ds['count'].name ='Acquisitions'


    dss[tc] = ds

# %%

count_cutoff = 20

figsize = (8, 11)
for tc in dss:
    ds = dss[tc]


    plt.figure()

    # mnum_counts.plot(marker='o')
    ds['count'].plot(hue=huedim, marker='o')


    plt.axhline(count_cutoff, linestyle='--')
    # plt.ylabel('Acquisitions')

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'acquisition_counts.png'))

#%%

for tc in dss:
    dss[tc] = dss[tc].where(dss[tc]['count'] > count_cutoff)
    dss[tc] = dss[tc].dropna(huedim, 'all')
    dss[tc] = dss[tc].dropna(coldim, 'all')

#%%


figsize = (8, 11)
for tc in dss:
    ds = dss[tc]
    ds = ds.drop('count')

    col_vals = ds.coords[coldim].values

    mosaic = []
    for var in list(ds.data_vars):
        mosaic_row = []
        for col_val in col_vals:
            mosaic_str = "{}_{}".format(var, col_val)
            mosaic_row.append(mosaic_str)

        mosaic.append(mosaic_row)


    fig, axes = plt.subplot_mosaic(mosaic, figsize=figsize, sharex=True)

    for var in list(ds.data_vars):
        for col_val in col_vals:
            ds_sel = ds[var].sel({coldim:col_val})

            mosaic_str = "{}_{}".format(var, col_val)

            ax = axes[mosaic_str]
            
            ds_sel.plot(ax=ax, hue=huedim)

            if var != 'i': ax.get_legend().remove()
            if var != 'i': ax.set_title('')

    tc_fig_dir = pjoin(figure_out_dir, tc)
    plt.savefig(pjoin(tc_fig_dir, 'raw_signals_hue.png'))
# %%
