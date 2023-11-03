#%%

from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

# %%

fp_dsst = pjoin(DIR_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.io.TFxr(fp_dsst).as_dsst()

fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.io.TFxr(fp_dst_coords).as_dsst()['coords']

fp_cuttimes = pjoin(DIR_PROC_DATA, 'cuttimes.csv')
df_ct = mhdpy.io.load_df_cuttimes(fp_cuttimes)

fp_expt_tws = r'C:\Users\aspitarl\Git\MHDLab\video-conversion\scripts\experiment_timewindows.csv'
df_exptw = mhdpy.io.load_df_cuttimes(fp_expt_tws)
# df_exptw = pd.read_csv(fp_expt_tws)

df_ct['date'] = df_ct['Start Time'].apply(lambda x: x.date())
df_exptw['date'] = df_exptw['Start Time'].apply(lambda x: x.date())
df_exptw = df_exptw.set_index('date')


#%%

coord_keys = [
    ('hvof', 'CC_K_massFrac_in'),
    ('filterwheel', 'Filter Position'),
    ('motor', 'Motor C Relative'),
    ('hvof','CC_equivalenceRatio')
]

das = []

for c in coord_keys:
    da = dsst[c[0]][c[1]]
    das.append(da)

ds_orig = xr.merge(das)

da = ds_orig['CC_K_massFrac_in']
ds_orig['CC_K_massFrac_in'] = da.where((da >= 0) & (da <= 1)).dropna('time','all')

da = ds_orig['CC_equivalenceRatio']
ds_orig['CC_equivalenceRatio'] = da.where((da >= 0.2) & (da <= 1.5)).dropna('time','all')


ds_orig

# %%

from mhdpy.plot.common import simple_ds_plot, tc_plot

simple_ds_plot(ds_orig)


#%%

plt.rcParams.update({'font.size': 12})

for date, df in df_ct.groupby('date'):

    timewindow = slice(
        Timestamp(df_exptw.loc[date]['Start Time']), 
        Timestamp(df_exptw.loc[date]['Stop Time'])
        )

    ds_plot = ds_orig.sel(time=timewindow)

    # Generate da_ct from df_cuttimes with date column
    
    df['tw'] = df.apply(lambda x: slice(x['Start Time'], x['Stop Time']), axis=1)
    da_ct = xr.DataArray(df.set_index(['Event'])['tw'])

    # TODO: something in here is messing with the order of the values
    fig = tc_plot(ds_plot, da_ct, legend_axes=1)

    fig.axes[0].set_yscale('log')
    fig.axes[3].set_ylim(0.5,1.5)

    fig.suptitle(date, fontdict={'size': 16})

    plt.tight_layout()

    plt.savefig(pjoin(DIR_FIG_OUT, 'tc_plot_{}.png'.format(date)))

# %%


da_ct