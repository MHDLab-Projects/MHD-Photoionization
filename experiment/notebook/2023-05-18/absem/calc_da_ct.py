#%%
from mhdpy.analysis.standard_import import *

from mhdpy.fileio import TFxr
from mhdpy.fileio.path import gen_path_date

data_folder = gen_path_date('2023-05-18')
dsst = TFxr(os.path.join(data_folder,'Processed_Data.tdms')).as_dsst()


#%%


from mhdpy.fileio.ct import load_df_cuttimes, extract_cuttime_list
df_cuttimes = load_df_cuttimes('cuttimes_tcm.csv').sort_values('Start Time').reset_index(drop=True)
cuttimes = extract_cuttime_list(df_cuttimes)
timewindow = slice(cuttimes[0].start, cuttimes[-1].stop)

timewindow

# %%

timewindow = slice(Timestamp('2023-05-18 21:44:30.206302208'), Timestamp('2023-05-18 22:11:50.342727424'), None)

dsst['hvof']['CC_K_massFrac_in'].sel(time=timewindow).plot()
# %%


from mhdpy.analysis.ct import gen_da_ct_data
from mhdpy.analysis.xr import fix_coord_grid_bins

tf = dsst['hvof']['CC_total_flow_in'].pint.to('g/s').pint.dequantify()
kwt = dsst['hvof']['CC_K_massFrac_in'].pint.dequantify()
# kwt = kwt*100

da_ct = gen_da_ct_data(cuttimes, {
    # 'tf': {'data': tf},
    'kwt': {'data': kwt}
    })

bins = [-0.0001, 0.0001, 0.0008, 0.0015, 0.0025 , 0.004, 0.0075, 0.015, 0.025, 0.05]
# da_ct = fix_coord_grid_bins(da_ct,'tf',3, round_value=3)
da_ct = fix_coord_grid_bins(da_ct,'kwt',bins, round_value=5)

# da_ct = da_ct.set_index(ct=['kwt','tf']).unstack('ct')

da_ct.coords['kwt'].attrs = dict(long_name='K Mass Pct.', units='%')
# da_ct.coords['tf'].attrs = dict(long_name='Total Mass Flow', units='g/s')


import pickle
with open('output/da_ct.pickle', 'wb') as f:
    pickle.dump(da_ct, f)

da_ct
# %%
