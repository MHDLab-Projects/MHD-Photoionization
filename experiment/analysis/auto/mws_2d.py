"""
Analysis of 2D parameter sweeps

#TODO: Add other dates/2d sweeps, improve plots, and convert to use main function like analsis_1d.py
"""
#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()



datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

figure_out_dir = pjoin(DIR_DATA_OUT, 'mws_2d', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
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

from mhdpy.analysis.mws import calc_mag_phase_AS

dss = {}

for tc in tcs:

    fp_in = pjoin(DIR_PROC_DATA, '{}.cdf'.format(tc))

    ds_in = xr.load_dataset(fp_in)
    ds_in = ds_in.sel(date=datestr).sel(run_num=1)

    ds_in = ds_in.drop('kwt')

    ds = calc_mag_phase_AS(ds_in).drop('mag_pp')


    tc_dim = [dim for dim in ds.dims if dim not in ['time','mnum']][0]
    mnum_counts = ds['i'].mean('time').groupby(tc_dim).count('mnum')

    ds = ds.mean('mnum', keep_attrs=True)

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
