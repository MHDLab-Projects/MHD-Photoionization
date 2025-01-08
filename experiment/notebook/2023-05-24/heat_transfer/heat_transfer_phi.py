# %% [markdown]
# # Simulation input calc and heat transfer fix for 2023-05-24 experiment
# 
# This notebook corrects the incorrect CC water flow scale discovered on 2023-06-07. see Equipment lab notebook page and DAQ channel Teams post. 
# 
# The simulation input file is calculated for the seed ramp (0.1, 0.31, 0.1 kwt %, 0.8 equivalence ratio, and 13 g/s total flow), which includes the fixed heat transfer values. 

# %%

from mhdpy.analysis.standard_import import *
create_standard_folders()
from mhdpy.plot import xr_errorbar, xr_errorbar_axes

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)

# %% [markdown]
# # Test Case analysis
# 
# select out time windows corresponding to different seed kwt. 

# %%

from mhdpy.fileio.ct import load_df_cuttimes

fp_cuttimes_phi = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_testcase_phi.csv')
df_cuttimes = load_df_cuttimes(fp_cuttimes_phi)

start_shift = np.timedelta64(1, 'm')
stop_shift = np.timedelta64(10, 's')

df_cuttimes['Start Time'] = df_cuttimes['Start Time'] + start_shift
df_cuttimes['Stop Time'] = df_cuttimes['Stop Time'] - stop_shift

cuttimes = df_cuttimes.ct.slice_list()

#%%

df_cuttimes = df_cuttimes.ct.add_data_stat(dsst['hvof']['CC_K_massFrac_in']*100, round_value=2)
df_cuttimes = df_cuttimes.ct.add_data_stat(dsst['hvof']['CC_equivalenceRatio'], round_value=1)

df_cuttimes = df_cuttimes.rename({'CC_K_massFrac_in': 'kwt', 'CC_equivalenceRatio': 'phi'}, axis=1)
df_cuttimes


#%%

from mhdpy.coords import assign_signal, unstack_multindexed_acq_dim
from mhdpy.coords.ct import downselect_acq_time

da_ht  = dsst['calorimetry', 'CC_heatTransfer']

da_ht = downselect_acq_time(da_ht, df_cuttimes, 'time')

da_signal = df_cuttimes.ct.column_to_coord_signal('kwt', da_ht.coords['time'])
da_ht = assign_signal(da_ht, da_signal, 'time')

da_signal = df_cuttimes.ct.column_to_coord_signal('phi', da_ht.coords['time'])
da_ht = assign_signal(da_ht, da_signal, 'time')

da_ht = unstack_multindexed_acq_dim(da_ht)['CC_heatTransfer']

da_ht

#%%

# region = get_region(da_ct, buffer=0.2)
region = df_cuttimes.ct.timewindow(buffer=0.2)

# %%

ds_plot = xr.merge([
    dsst['calorimetry'][['CC_heatTransfer', 'CC_water_dT']],
    dsst['hvof'][['CC_total_flow_in', 'CC_K_massFrac_in', 'CC_equivalenceRatio']]
])

from mhdpy.plot import simple_ds_plot, tc_plot

df_cuttimes['plot_label'] = df_cuttimes['kwt'].astype(str) + ' kwt, ' + df_cuttimes['phi'].astype(str) + ' phi'
fig =  tc_plot(ds_plot.sel(time=region), df_cuttimes, legend_axes=3, legend_column='plot_label')

# # fig.axes[0].get_legend().remove()
# # fig.axes[3].legend()

plt.savefig(pjoin(DIR_FIG_OUT,'phi_heattransfer_tcplot.png'))

#%%
da_ht_mean = da_ht.mean('mnum').pint.dequantify()
da_ht_std = da_ht.std('mnum').pint.dequantify()

g = da_ht_mean.plot(hue='kwt', marker='o')

#%%

#%%

xr_errorbar(da_ht_mean, da_ht_std, huedim='kwt')

plt.savefig(pjoin(DIR_FIG_OUT, 'phi_heattransfer.png'))

# %% [markdown]
# # Make simulation input file

# %%

#TODO: combine with 04-07 Microwave transmission, move to mhdpy 


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


key_sel = ['hvof', 'calorimetry']
sheet_names = {'hvof': "HVOF Process Inputs", 'calorimetry': "Calorimetry"}

dsst_stats = {key: dsst[key] for key in key_sel}


from mhdpy.xr_utils import calc_stats

#TODO: automate output name

dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input_2023-05-24_phi.xlsx'), engine='xlsxwriter') as writer:
    for key in dsst_stats:
        ds = dsst_stats[key]

        # Assign a time signal of the test case names, then calculate stats
        tc_time_signal = df_cuttimes.reset_index().ct.column_to_coord_signal('index', ds.coords['time'])
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

    recipe_em = pd.Series(dsst['hvof']['CC_K_massFrac_in'].attrs)
    # recipe_em = recipe_em.to_frame().T


    # This is to remove the attributes associated with CC_K_massFrac_in signal. Is there a better way to 'carry' the recipe?
    col_select=['em_M_tween80', 'em_M_span80', 'em_M_surf', 'em_M_brine','em_M_water','em_M_fuel', 'em_M_total',
       'f_fuel', 'f_K2CO3', 'rho_em',
       'percent_K_to_K2CO3', 'f_water']
    
    recipe_em = recipe_em.loc[col_select]


    recipe_split = recipe_em.str.replace(' / ', '/').str.split(' ', expand=True)
    recipe_split.columns = ['val', 'unit']
    recipe_split = recipe_split.T
    


    recipe_split.to_excel(writer, sheet_name='Emulsion Recipe')

    df_cuttimes.to_excel(writer, sheet_name='Experiment Time Windows')
    


# %%
