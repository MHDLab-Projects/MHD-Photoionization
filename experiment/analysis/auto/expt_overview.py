#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

#have time axes show up in PST
plt.rcParams['timezone'] = 'US/Pacific'
import matplotlib.dates as mdates

# %%

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst(convert_to_PT=False)

fp_dst_coords = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst(convert_to_PT=False)['coords']

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'ct_sequence.csv')
df_ct = mhdpy.fileio.load_df_cuttimes(fp_cuttimes)

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws)
# df_exptw = pd.read_csv(fp_expt_tws)

df_ct['date'] = df_ct['Start Time'].apply(lambda x: x.date())
df_exptw['date'] = df_exptw['Start Time'].apply(lambda x: x.date())
df_exptw = df_exptw.set_index('date')


#%%

coord_keys = [
    ('hvof', 'CC_K_massFrac_in'),
    ('filterwheel', 'Power_Relative'),
    ('motor', 'Motor C Relative'),
    ('hvof','CC_equivalenceRatio')
]

das = []

for c in coord_keys:
    da = dsst[c[0]][c[1]]
    das.append(da)

ds_orig = xr.merge(das)

da = ds_orig['CC_K_massFrac_in']
ds_orig['CC_K_massFrac_in'] = da.where((da >= 0) & (da <= 1)).dropna('time', how='all')

da = ds_orig['CC_equivalenceRatio']
ds_orig['CC_equivalenceRatio'] = da.where((da >= 0.2) & (da <= 1.5)).dropna('time', how='all')

# Needed to get the legend colors to line up correctly
ds_orig['Power_Relative'] = ds_orig['Power_Relative'].ffill('time')
ds_orig['Motor C Relative'] = ds_orig['Motor C Relative'].ffill('time')

ds_orig

#%%


# %%

from mhdpy.plot import simple_ds_plot, tc_plot

simple_ds_plot(ds_orig)


#%%
import textwrap

plt.rcParams.update({'font.size': 14})

for date, df in df_ct.groupby('date'):

    timewindow = slice(
        Timestamp(df_exptw.loc[date]['Start Time']), 
        Timestamp(df_exptw.loc[date]['Stop Time'])
        )

    ds_plot = ds_orig.sel(time=timewindow)
    # TODO: something in here is messing with the order of the values
    fig = tc_plot(ds_plot, df, legend_axes=0, figsize=(8,10), legend_column='Event')

    # fig.axes[1].get_legend().set_bbox_to_anchor((1, 1))
    # fig.axes[1].get_legend.set(bbox_to_anchor=(0.5, 1.5), loc="lower center", bbox_transform=fig.transFigure)   
    fig.axes[0].get_legend().set_bbox_to_anchor([0.65,0.8], transform=fig.transFigure)
    fig.subplots_adjust(top=0.8, left=0.2)  # Adjust the top spacing

    fig.axes[0].set_yscale('log')
    fig.axes[1].set_yscale('log')
    fig.axes[1].set_ylim(0.001,2)
    fig.axes[3].set_ylim(0.5,1.5)

    fig.suptitle(date, fontdict={'fontsize': 16})

    for ax in fig.axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_ylabel(textwrap.fill(ax.get_ylabel(), 20))
        ax.set_xlabel('')

    fig.axes[-1].set_xlabel('Time (PST)')

    # plt.tight_layout()

    plt.savefig(pjoin(DIR_FIG_OUT, 'tc_plot_{}.png'.format(date)))

    # break

# %%


# %%

output_dir = pjoin(DIR_FIG_OUT, 'other_signals')
os.makedirs(output_dir, exist_ok=True)

coord_keys = [
    ('hvof', 'CC_fuel_flow_in'),
    ('hvof', 'CC_o2_flow_in'),
    ('hvof','CC_P'),
    ('calorimetry', 'CC_water_flow_in'),
    ('calorimetry', 'CC_water_T_in'),
    ('calorimetry', 'CC_water_T_out'),
]

das = []

for c in coord_keys:
    da = dsst[c[0]][c[1]]
    das.append(da)

ds_orig = xr.merge(das)

#%%

#TODO: why are these missing?
ds_orig['CC_o2_flow_in'].attrs['long_name'] = 'O2 Flow Rate'
ds_orig['CC_water_flow_in'].attrs['long_name'] = 'Cooling Water Flow Rate'

#%%


plt.rcParams.update({'font.size': 12})

for date, df in df_ct.groupby('date'):

    timewindow = slice(
        Timestamp(df_exptw.loc[date]['Start Time']), 
        Timestamp(df_exptw.loc[date]['Stop Time'])
        )

    ds_plot = ds_orig.sel(time=timewindow).pint.dequantify()

    # TODO: something in here is messing with the order of the values
    fig = tc_plot(ds_plot, df, legend_axes=None, figsize=(8,12), legend_column='Event')

    # fig.axes[0].set_yscale('log')
    # fig.axes[3].set_ylim(0.5,1.5)

    # fig.axes[-1].get_legend().remove()

    fig.suptitle(date, fontdict={'fontsize': 16})

    for ax in fig.axes:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_ylabel(textwrap.fill(ax.get_ylabel(), 20))

        ax.set_xlabel('')


    fig.axes[-1].set_xlabel('Time (PST)')
    plt.tight_layout()

    plt.savefig(pjoin(output_dir, 'tc_plot_{}.png'.format(date)))


# %%
