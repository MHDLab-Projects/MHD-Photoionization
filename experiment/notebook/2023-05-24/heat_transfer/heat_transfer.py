# %% [markdown]
# # Simulation input calc and heat transfer fix for 2023-05-24 experiment
# 
# This notebook corrects the incorrect CC water flow scale discovered on 2023-06-07. see Equipment lab notebook page and DAQ channel Teams post. 
# 
# The simulation input file is calculated for the seed ramp (0.1, 0.31, 0.1 kwt %, 0.8 equivalence ratio, and 13 g/s total flow), which includes the fixed heat transfer values. 

# %%

from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()


# %% [markdown]
# # Test Case analysis
# 
# select out time windows corresponding to different seed kwt. 

# %%

from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
from mhdpy.coords.ct import gen_da_ct_data, assign_tc_general, get_region


df_cuttimes = load_df_cuttimes('cuttimes_seedramp.csv')

cuttimes = extract_cuttime_list(df_cuttimes)


#TODO: simplify this in mhdpy

da_ct = gen_da_ct_data(cuttimes, dim_dict={
    # 'phi': {'data': dsst['hvof']['CC_equivalenceRatio']},
    'kwt': {'data' : dsst['hvof']['CC_K_massFrac_in']*100},
    # 'tf': {'data': dsst['hvof']['CC_total_flow_in']},


})

da_ct = da_ct.rename(ct='kwt').set_index(kwt='kwt')

da_ct = da_ct.assign_coords(kwt = [round(kwt,3) for kwt in da_ct.coords['kwt'].values])

da_ct.coords['kwt'].attrs['long_name'] = 'K wt%'

region = get_region(da_ct, buffer=0.2)

da_ct


# %%

ds_plot = xr.merge([
    dsst['calorimetry'][['CC_heatTransfer', 'CC_water_dT']],
    dsst['hvof'][['CC_total_flow_in', 'CC_K_massFrac_in', 'CC_equivalenceRatio']]
])

from mhdpy.plot.common import simple_ds_plot, tc_plot


region = slice(Timestamp('2023-05-24 19:51:43.989046272'), Timestamp('2023-05-24 20:47:57.945108480'), None)

fig =  tc_plot(ds_plot.sel(time=region), da_ct, )

plt.savefig(pjoin(DIR_FIG_OUT,'seedramp_heattransfer.png'))

# %% [markdown]
# # Make simulation input file

# %%

#TODO: combine with 04-07 microwave scattering, move to mhdpy 


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

da_ct = xr.DataArray(cuttimes, coords = {'tc': df_cuttimes.index.values}, dims = ['tc'])

from mhdpy.coords.ct import assign_tc_general
from mhdpy.xr_utils import calc_stats


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
    col_select=['em_M_tween80', 'em_M_span80', 'em_M_surf', 'em_M_brine','em_M_water','em_M_fuel', 'em_M_total',
       'f_fuel', 'f_K2CO3', 'rho_em',
       'percent_K_to_K2CO3', 'f_water']
    
    recipe_em = recipe_em.loc[col_select]


    recipe_split = recipe_em.str.replace(' / ', '/').str.split(' ', expand=True)
    recipe_split.columns = ['val', 'unit']
    recipe_split = recipe_split.T
    


    recipe_split.to_excel(writer, sheet_name='Emulsion Recipe')

    df_cuttimes.to_excel(writer, sheet_name='Experiment Time Windows')
    


