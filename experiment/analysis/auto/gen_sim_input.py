#%%
from mhdpy.analysis.standard_import import *
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst()


fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

dst_coords

#%%

# Combustion chamber pressure sensor was faulty at the beginning of 04-07-2023. This removes those datapoints
# Note. the starting time is just by eye here. 

tw_CC_P_valid = slice(pd.Timestamp('2023-04-07 19:45:00'), pd.Timestamp('2024-05-24 23:59:59'))
dsst['hvof']['CC_P'] = dsst['hvof']['CC_P'].sel(time=tw_CC_P_valid)


# %%

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list

fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'cuttimes_kwt.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes).sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

df_cuttimes['date'] = df_cuttimes['Start Time'].dt.date

cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

# Convert to local time and drop all timezone information. 
# https://stackoverflow.com/questions/61802080/excelwriter-valueerror-excel-does-not-support-datetime-with-timezone-when-savin


# df_cuttimes['Start Time'] = df_cuttimes['Start Time'].dt.tz_localize('UTC').dt.tz_convert('US/Pacific').dt.tz_localize(None)
# df_cuttimes['Stop Time'] = df_cuttimes['Stop Time'].dt.tz_localize('UTC').dt.tz_convert('US/Pacific').dt.tz_localize(None)

#%%

tc_names = []

for idx, row in df_cuttimes.iterrows():
    timeslice = slice(row['Start Time'], row['Stop Time'])
    da = dst_coords.sel(time=timeslice)['kwt'].dropna('time')
    val_list = list(set(da.values))

    assert len(val_list) == 1, f'Not all values are the same for {idx}'

    val = round(val_list[0]*1, 2)

    tc_names.append("{}_{}".format(row['date'], val))

df_cuttimes.index = tc_names

df_cuttimes


df_cuttimes = df_cuttimes.reset_index()

# df_cuttimes.groupby()
for idx, row in df_cuttimes.groupby('index'):

    mnums = range(1, len(row)+1)

    df_cuttimes.loc[row.index, 'index'] = [f'{idx}_{mnum}' for mnum in mnums]

    pass


df_cuttimes = df_cuttimes.set_index('index')



#%%

#TODO: Rename these in lab config
long_name_dict = {
'hvof': {
'CC_o2_flow_in': 'Oxygen Flow'
},
'calorimetry':{
'CC_water_flow_in': "JP Cooling Water Flow",
'CC_heatTransfer': "CC Wall Heat Transfer",
}
}

for group in long_name_dict:
    group_long_name_dict = long_name_dict[group]
    for var in group_long_name_dict:
        dsst[group][var].attrs['long_name'] = group_long_name_dict[var]


dsst['calorimetry']['CC_water_T_in'] = dsst['calorimetry']['CC_water_T_in'].pint.to('K')
dsst['calorimetry']['CC_water_T_out'] = dsst['calorimetry']['CC_water_T_out'].pint.to('K')
dsst['calorimetry']['CC_water_dT'] = dsst['calorimetry']['CC_water_dT'].pint.to('K')

first_vars = ['CC_total_flow_in','CC_equivalenceRatio','CC_K_massFrac_in']
dsst['hvof'] = dsst['hvof'][[*first_vars, *[var for var in dsst['hvof'].data_vars if var not in first_vars]]]

# dsst['syringe'] = dsst['syringe'][['syringe_em_flow_out', 'syringe_fuel_flow_out', 'syringe_water_flow_out', 'syringe_K2CO3_flow_out', 'syringe_K_flow_out']]

#%%

key_sel = ['hvof', 'calorimetry']
sheet_names = {'hvof': "HVOF Process Inputs", 'calorimetry': "Calorimetry"}

dsst_stats = {key: dsst[key] for key in key_sel}

da_ct = xr.DataArray(cuttimes, coords = {'tc': df_cuttimes.index.values}, dims = ['tc'])

from mhdpy.coords.ct import assign_tc_general
from mhdpy.xr_utils import calc_stats

#%%


def extract_em_recipe_table(recipe_em):

    col_select=['em_M_tween80', 'em_rho', 'em_M_surf', 'em_M_brine', 'em_M_total',
         'em_f_fuel', 'em_f_K2CO3', 'f_fuel', 'f_K2CO3', 'f_water', 'f_surf', 'rho_em',
         'percent_K_to_K2CO3', 'f_water']
        
    col_select = [col for col in col_select if col in recipe_em]

    recipe_em = recipe_em.loc[col_select]

    recipe_split = recipe_em.str.replace(' / ', '/').str.split(' ', expand=True)
    recipe_split.columns = ['val', 'unit']
    recipe_split = recipe_split.T

    return recipe_split

dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input_all.xlsx'), engine='xlsxwriter') as writer:
    for key in dsst_stats:
        ds = dsst_stats[key]
        ds = assign_tc_general(ds, da_ct)

        # ds_out = ds.mean('time', keep_attrs=True)#.to_dataframe()
        ds_out = calc_stats(ds)

        df_out = ds_out.pint.dequantify().to_dataframe()

        units = [ds_out[var].pint.units for var in ds_out.data_vars]
        long_names = [ds_out[var].attrs['long_name'] if 'long_name' in ds_out[var].attrs else '' for var in ds_out.data_vars]


        df_out.columns = [df_out.columns, units, long_names]
        df_out.index.names = ['Statistic', 'Test Case']
        df_out.columns.names = ['Name','Unit','Long Name']

        df_out.to_excel(writer, sheet_name=sheet_names[key])

    recipe_em = pd.Series(dsst['hvof']['CC_K_massFrac_in'].attrs)
    recipe_split = extract_em_recipe_table(recipe_em)
    recipe_split.to_excel(writer, sheet_name='Emulsion Recipe')

    df_cuttimes.to_excel(writer, sheet_name='Experiment Time Windows')
    
# %%


from mhdpy.fileio.tdms import tdms2ds
from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim

fp_dst_coords = pjoin(DIR_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdpy.fileio.TFxr(fp_dst_coords).as_dsst()['coords']

dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input_mean.xlsx'), engine='xlsxwriter') as writer:
    for key in dsst_stats:
        ds = dsst_stats[key]

        dss = []

        for idx, row in df_cuttimes.iterrows():
            timeslice = slice(row['Start Time'], row['Stop Time'])

            ds_sel = ds.sel(time=timeslice)
            ds_sel = assign_signal(ds_sel, dst_coords['kwt'].round(2), timeindex='time')

            dss.append(ds_sel)

        ds = xr.concat(dss, dim='time')

        ds = unstack_multindexed_acq_dim(ds)

        ds_out = calc_stats(ds, stat_dim='mnum')

        df_out = ds_out.pint.dequantify().to_dataframe()

        units = [ds_out[var].pint.units for var in ds_out.data_vars]
        long_names = [ds_out[var].attrs['long_name'] if 'long_name' in ds_out[var].attrs else '' for var in ds_out.data_vars]


        df_out.columns = [df_out.columns, units, long_names]
        df_out.index.names = ['Statistic', 'Test Case']
        df_out.columns.names = ['Name','Unit','Long Name']

        df_out.to_excel(writer, sheet_name=sheet_names[key])

    recipe_em = pd.Series(dsst['hvof']['CC_K_massFrac_in'].attrs)
    recipe_split = extract_em_recipe_table(recipe_em)
    recipe_split.to_excel(writer, sheet_name='Emulsion Recipe')

    df_cuttimes.to_excel(writer, sheet_name='Experiment Time Windows')



#%%

df_ct = df_cuttimes.copy()

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'experiment_timewindows.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws)
# df_exptw = pd.read_csv(fp_expt_tws)

df_ct['date'] = df_ct['Start Time'].apply(lambda x: x.date())
df_exptw['date'] = df_exptw['Start Time'].apply(lambda x: x.date())
df_exptw = df_exptw.set_index('date').sort_index()

df_ct

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
ds_orig['CC_K_massFrac_in'] = da.where((da >= 0) & (da <= 1)).dropna('time', how='all')


#%%

plt.rcParams.update({'font.size': 16})

fig, axes = plt.subplots(len(df_exptw), 1, figsize=(10, 5*len(df_exptw)))

# for i, date, df in df_ct.groupby('date'):
for i, date in enumerate(df_exptw.index):



    tw_expt = df_exptw.loc[date]
    tw_expt = slice(tw_expt['Start Time'], tw_expt['Stop Time'])

    dsst['hvof']['CC_K_massFrac_in'].sel(time=tw_expt).plot(ax=axes[i])

    df_ct_date = df_ct[df_ct['date'] == date]

    for idx, row in df_ct_date.iterrows():
        axes[i].axvspan(row['Start Time'], row['Stop Time'], color='green', alpha=0.3)
    
    axes[i].set_title(date)

    if i != len(df_exptw)-1:
        axes[i].set_xlabel('')
        axes[i].set_xticklabels([])

#%%

df_exptw