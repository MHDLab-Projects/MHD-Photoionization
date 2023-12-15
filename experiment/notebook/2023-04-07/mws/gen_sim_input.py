#%%
from mhdpy.analysis.standard_import import *

datestr = '2023-04-07'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
# #%%

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
df_cuttimes = load_df_cuttimes('cuttimes_sim.csv').sort_values('Start Time').reset_index(drop=True)
df_cuttimes = df_cuttimes.set_index('Event')

cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

# Convert to local time and drop all timezone information. 
# https://stackoverflow.com/questions/61802080/excelwriter-valueerror-excel-does-not-support-datetime-with-timezone-when-savin


df_cuttimes['Start Time'] = df_cuttimes['Start Time'].dt.tz_localize('UTC').dt.tz_convert('US/Pacific').dt.tz_localize(None)
df_cuttimes['Stop Time'] = df_cuttimes['Stop Time'].dt.tz_localize('UTC').dt.tz_convert('US/Pacific').dt.tz_localize(None)


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


#%%

key_sel = ['hvof', 'calorimetry']
sheet_names = {'hvof': "HVOF Process Inputs", 'calorimetry': "Calorimetry"}

dsst_stats = {key: dsst[key] for key in key_sel}

da_ct = xr.DataArray(cuttimes, coords = {'tc': df_cuttimes.index.values}, dims = ['tc'])

from mhdpy.analysis.ct import assign_tc_general
from mhdpy.analysis.xr import calc_stats


dss_out = []
with pd.ExcelWriter(pjoin(DIR_DATA_OUT, 'sim_input.xlsx'), engine='xlsxwriter') as writer:
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
    # recipe_em = recipe_em.to_frame().T


    # This is to remove the attributes associated with CC_K_massFrac_in signal. Is there a better way to 'carry' the recipe?
    col_select=['em_M_tween80', 'em_rho', 'em_M_surf', 'em_M_brine', 'em_M_total',
       'em_f_fuel', 'em_f_K2CO3', 'f_fuel', 'f_K2CO3', 'rho_em',
       'percent_K_to_K2CO3', 'f_water']
    
    recipe_em = recipe_em.loc[col_select]


    recipe_split = recipe_em.str.replace(' / ', '/').str.split(' ', expand=True)
    recipe_split.columns = ['val', 'unit']
    recipe_split = recipe_split.T
    


    recipe_split.to_excel(writer, sheet_name='Emulsion Recipe')

    df_cuttimes.to_excel(writer, sheet_name='Experiment Time Windows')
    
#%%


#%%

# Testing that the total flow makes sense, but won't be exactly equal after the average...But did remoev argon during post processing. 

# df_out['CC_fuel_flow_in'] + df_out['CC_o2_flow_in'] + df_out['syringe_em_flow_out']

# df_out['CC_total_flow_in']
# dsst['JP_ar_flow_out']

#%%


# %%

#%%


