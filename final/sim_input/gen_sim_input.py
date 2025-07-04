#%%
from mhdlab.analysis.standard_import import *
create_standard_folders()

#TODO: move plotting to a separate file
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
plt.rcParams.update({'font.size': 8, 'timezone': 'America/Los_Angeles'})

from mhdlab.xr_utils import calc_stats
from mhdlab.coords import assign_signal, unstack_multindexed_acq_dim

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdlab.fileio.TFxr(fp_dsst).as_dsst(convert_to_PT=False)

fp_dst_coords = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdlab.fileio.TFxr(fp_dst_coords).as_dsst(convert_to_PT=False)['coords']

dst_coords
#%%

dsst.keys()

#%%

# Combustion chamber pressure sensor was faulty at the beginning of 04-07-2023. This removes those datapoints
# Note. the starting time is just by eye here. 

tw_CC_P_valid = slice(pd.Timestamp('2023-04-07 19:45:00'), pd.Timestamp('2024-05-24 23:59:59'))
dsst['hvof']['CC_P'] = dsst['hvof']['CC_P'].sel(time=tw_CC_P_valid)

# %%

from mhdlab.fileio.ct import load_df_cuttimes

# Test case cuttimes
fp_cuttimes = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_testcase_kwt.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes).sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

# Sequence cuttimes. Downselect to only 53x

#TODO: rename the file to sequence timewindows
fp_sequence_tw = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_sequence.csv')

df_sequence_tw = load_df_cuttimes(fp_sequence_tw)
# add tc and run_num columns. tc is everything before the last underscore in 'Event' column
df_sequence_tw['tc'] = df_sequence_tw['Event'].str.extract(r'(.*)_\d$')
df_sequence_tw['run_num'] = df_sequence_tw['Event'].str.extract(r'.*_(\d)$').astype(int)

df_sequence_tw = df_sequence_tw.where(df_sequence_tw['tc'].isin(['53x', '516_pos'])).dropna()

# Experiment timewindows

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws)
df_exptw = df_exptw.set_index('date').sort_index()


#%%

# Filter the cuttimes to only include the ones that are in the sequence timewindows


def filter_cuttimes(df_ct_testcase, df_ct_sequence, mode='both'):
    """
    Filter the testcase cuttimes to only those that are within the sequence timewindows 
    df_ct_testcase: dataframe of cuttimes for the test case
    df_ct_sequence: dataframe of cuttimes for the sequences, should be longer than df_ct_testcase

    mode: 'both' or 'either'
    'both' will only include cuttimes that are entirely within the sequence timewindow
    'either' will include cuttimes that are partially within the sequence timewindow
    """


    df_ct_testcase_filtered = pd.DataFrame()

    for _, row in df_ct_sequence.iterrows():
        start_time, stop_time = row['Start Time'], row['Stop Time']

        if mode == 'both':
            df_temp = df_ct_testcase[(df_ct_testcase['Start Time'] >= start_time) & (df_ct_testcase['Stop Time'] <= stop_time)]
        elif mode == 'either':
            df_temp = df_ct_testcase[((df_ct_testcase['Start Time'] >= start_time) & (df_ct_testcase['Start Time'] <= stop_time)) | 
                                  ((df_ct_testcase['Stop Time'] >= start_time) & (df_ct_testcase['Stop Time'] <= stop_time))]
        else:
            raise ValueError("Invalid mode. Choose either 'both' or 'either'.")

        df_ct_testcase_filtered = pd.concat([df_ct_testcase_filtered, df_temp])

    return df_ct_testcase_filtered


df_cuttimes = filter_cuttimes(df_cuttimes, df_sequence_tw, mode='either')
df_cuttimes

#%%

# Manually add 516 pos 1 from df_sequence_tw to df_cuttimes. Add as the final row in df_cuttimes

df_516_pos1 = df_sequence_tw[df_sequence_tw['tc'] == '516_pos'].iloc[0]
df_516_pos1

df_516_pos1 = df_516_pos1.to_frame().T
df_516_pos1['Event'] = '516_pos1'

df_cuttimes.loc[len(df_cuttimes)] = df_516_pos1.iloc[0]

# df_cuttimes = df_cuttimes

#%%

"Form the test case names"

kwts = []
phis = []

for idx, row in df_cuttimes.iterrows():
    timeslice = slice(row['Start Time'], row['Stop Time'])
    da = dst_coords.sel(time=timeslice)['kwt'].dropna('time')
    val_list = list(set(da.values))

    assert len(val_list) == 1, f'Not all values are the same for {idx}'

    val = round(val_list[0]*1, 2)

    kwts.append(val)

    da = dst_coords.sel(time=timeslice)['phi'].dropna('time')
    val_list = list(set(da.values))

    assert len(val_list) == 1, f'Not all values are the same for {idx}'

    val = round(val_list[0]*1, 2)

    phis.append(val)




df_cuttimes['kwt'] = kwts
df_cuttimes['phi'] = phis

# groupby kwt and date columns and add a measurement number column

df_cuttimes['mnum'] = df_cuttimes.groupby(['kwt', 'phi', 'date']).cumcount()+1
df_cuttimes['repeat'] = df_cuttimes['date'].astype(str) + '_' + df_cuttimes['mnum'].astype(str)
df_cuttimes['tc'] = df_cuttimes['kwt'].astype(str) + '_' + df_cuttimes['phi'].astype(str) + '_' + df_cuttimes['repeat'].astype(str)


#TODO: can we have kwt as it's own multindex level? Mainly for display in output
# df_cuttimes = df_cuttimes.set_index(['kwt', 'repeat'])

df_cuttimes = df_cuttimes.set_index('tc')
df_cuttimes = df_cuttimes.sort_index()

cuttimes = df_cuttimes.ct.slice_list()


df_cuttimes

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
dsst['other_temperature'] = xr.merge([
    dsst['o2']['JP_o2_T_out'].pint.to('K'),
    dsst['tc1']['ambient_T'].pint.to('K')
])

#%%

key_sel = ['hvof', 'calorimetry', 'other_temperature']
sheet_names = {'hvof': "HVOF Process Inputs", 'calorimetry': "Calorimetry", 'other_temperature': "Other Temperatures"}

dsst_stats = {key: dsst[key] for key in key_sel}



#%%

# Plot the experiment coordinates

coord_keys = [
    ('hvof', 'CC_K_massFrac_in'),
    ('filterwheel', 'Filter Position'),
    ('motor', 'Motor C Relative'),
    ('hvof','CC_equivalenceRatio')
]

das = [dsst[c[0]][c[1]] for c in coord_keys]
ds_orig = xr.merge(das)

da = ds_orig['CC_K_massFrac_in']
ds_orig['CC_K_massFrac_in'] = da.where((da >= 0) & (da <= 1)).dropna('time', how='all')

plt.rcParams.update({'font.size': 16})

fig, axes = plt.subplots(len(df_exptw), 1, figsize=(10, 5*len(df_exptw)), sharex=False)

for i, date in enumerate(df_exptw.index):

    tw_expt = df_exptw.loc[date]
    tw_expt = slice(tw_expt['Start Time'], tw_expt['Stop Time'])

    dsst['hvof']['CC_K_massFrac_in'].sel(time=tw_expt).plot(ax=axes[i])

    df_cuttimes_date = df_cuttimes[df_cuttimes['date'] == date]

    for idx, row in df_cuttimes_date.iterrows():
        axes[i].axvspan(row['Start Time'], row['Stop Time'], color='green', alpha=0.3)

    df_sequence_date = df_sequence_tw[df_sequence_tw['date'] == date]

    for idx, row in df_sequence_date.iterrows():
        axes[i].axvspan(row['Start Time'], row['Stop Time'], color='gray', alpha=0.3)
    
    axes[i].set_title(date)

    if i != len(df_exptw)-1:
        axes[i].set_xlabel('')
        axes[i].set_xticklabels([])


plt.savefig(pjoin(DIR_FIG_OUT, 'sim_input_timewindows.png'), dpi=300, bbox_inches='tight')

#%%

# Make a plot of each individual datarray in dsst with the timewindows. 
downselect_keys = {
    'hvof': ['CC_total_flow_in', 'CC_fuel_flow_in', 'CC_o2_flow_in', 'CC_em_flow_in', 'CC_P'],
    'calorimetry': ['CC_water_flow_in', 'CC_water_T_in', 'CC_water_T_out', 'CC_heatTransfer'],
    'other_temperature': ['JP_o2_T_out', 'ambient_T']
}

output_dir = pjoin(DIR_FIG_OUT, 'individual_signals')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for i, date in enumerate(df_exptw.index):
    tw_expt = df_exptw.loc[date]
    tw_expt = slice(tw_expt['Start Time'], tw_expt['Stop Time'])

    df_cuttimes_sel = df_cuttimes[df_cuttimes['date'] == date]
    df_sequence_tw_sel = df_sequence_tw[df_sequence_tw['date'] == date]

    for key in dsst_stats:

        ds = dsst_stats[key].sel(time=tw_expt)

        ds = ds[downselect_keys[key]]

        for var in ds.data_vars:
            da = ds[var]

            fig, ax = plt.subplots()

            da.plot(ax=ax)

            for idx, row in df_cuttimes_sel.iterrows():
                ax.axvspan(row['Start Time'], row['Stop Time'], color='green', alpha=0.3)

            for idx, row in df_sequence_tw_sel.iterrows():
                ax.axvspan(row['Start Time'], row['Stop Time'], color='gray', alpha=0.3)

            ax.set_title(var)



            plt.savefig(pjoin(output_dir, f'{date}_{key}_{var}.png'), dpi=300, bbox_inches='tight')

            # break

        

#%%

def extract_params(ds_hvof):
    #TODO: name this dict to be more general
    recipe_em = pd.Series(ds_hvof['CC_K_massFrac_in'].attrs)

    col_select=['em_M_tween80', 'em_rho', 'em_M_surf', 'em_M_brine', 'em_M_total',
         'em_f_fuel', 'em_f_K2CO3', 'f_fuel', 'f_K2CO3', 'f_water', 'f_surf', 'rho_em',
         'percent_K_to_K2CO3', 'f_water']
        
    col_select = [col for col in col_select if col in recipe_em]

    recipe_em = recipe_em.loc[col_select]

    recipe_em['calorific value'] = ds_hvof['CC_heatInput'].attrs['calorific value']
    recipe_em['CH_weight_ratio'] = ds_hvof['CC_equivalenceRatio'].attrs['CHweightratio']
    recipe_em['N_Carbon'] = ds_hvof['CC_equivalenceRatio'].attrs['N_Carbon']

    #TODO: hacky way to separate out value and units. 
    recipe_split = recipe_em.str.replace(' / ', '/').str.replace(' * ', '*', regex=False).str.split(' ', expand=True)
    recipe_split.columns = ['val', 'unit']
    recipe_split = recipe_split.T

    return recipe_split

dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input_all.xlsx'), engine='xlsxwriter') as writer:
    for key in dsst_stats:
        ds = dsst_stats[key]

        # Assign a time signal of the test case names, then calculate stats
        tc_time_signal = df_cuttimes.reset_index().ct.column_to_coord_signal('tc', ds.coords['time'])
        ds = ds.sel(time=tc_time_signal.time)
        ds = assign_signal(ds, tc_time_signal, 'time')
        ds = unstack_multindexed_acq_dim(ds, )

        ds_out = calc_stats(ds, stat_dim='mnum')

        df_out = ds_out.pint.dequantify().to_dataframe()

        units = [ds_out[var].pint.units for var in ds_out.data_vars]
        long_names = [ds_out[var].attrs['long_name'] if 'long_name' in ds_out[var].attrs else '' for var in ds_out.data_vars]


        df_out.columns = [df_out.columns, units, long_names]
        df_out.index.names = ['Statistic', 'Test Case']
        df_out.columns.names = ['Name','Unit','Long Name']

        df_out.to_excel(writer, sheet_name=sheet_names[key])

    params_table = extract_params(dsst['hvof'])
    params_table.to_excel(writer, sheet_name='Parameters')

    # Convert the start and stop time columns to PST 

    df_cuttimes_out = df_cuttimes.copy(deep=True)

    df_cuttimes_out['Start Time'] = df_cuttimes_out['Start Time'].dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.tz_localize(None)
    df_cuttimes_out['Stop Time'] = df_cuttimes_out['Stop Time'].dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.tz_localize(None)

    df_cuttimes_out.to_excel(writer, sheet_name='Experiment Time Windows')

#%%
    
df_out
    
# %%

import subprocess
# Get the current commit hash
commit_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip().decode('utf-8')




from mhdlab.fileio.tdms import tdms2ds
from mhdlab.coords import assign_signal, unstack_multindexed_acq_dim

fp_dst_coords = pjoin(DIR_EXPT_PROC_DATA, 'dst_coords.tdms')
dst_coords = mhdlab.fileio.TFxr(fp_dst_coords).as_dsst(convert_to_PT=False)['coords']

dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input_mean.xlsx'), engine='xlsxwriter') as writer:
    for key in dsst_stats:
        ds = dsst_stats[key]

        dss = []

        for idx, row in df_cuttimes.iterrows():
            timeslice = slice(row['Start Time'], row['Stop Time'])

            ds_sel = ds.sel(time=timeslice)
            ds_sel = assign_signal(ds_sel, dst_coords['kwt'].round(3), timeindex='time')
            ds_sel = assign_signal(ds_sel, dst_coords['phi'].round(3), timeindex='time')

            dss.append(ds_sel)

        ds = xr.concat(dss, dim='time')

        ds = unstack_multindexed_acq_dim(ds)

        ds_out = calc_stats(ds, stat_dim='mnum')

        df_out = ds_out.pint.dequantify().to_dataframe()

        units = [ds_out[var].pint.units for var in ds_out.data_vars]
        long_names = [ds_out[var].attrs['long_name'] if 'long_name' in ds_out[var].attrs else '' for var in ds_out.data_vars]


        df_out.columns = [df_out.columns, units, long_names]
        df_out.index.names = ['Statistic', 'K wt', 'Phi']
        df_out.columns.names = ['Name','Unit','Long Name']

        df_out.to_excel(writer, sheet_name=sheet_names[key])

    params_table = extract_params(dsst['hvof'])

    # Add the commit hash to the params_table DataFrame
    params_table['Commit Hash'] = commit_hash

    params_table.to_excel(writer, sheet_name='Parameters')

