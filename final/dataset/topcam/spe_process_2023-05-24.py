#%%

from mhdlab.analysis.standard_import import *
import pi_paper_utils
create_standard_folders()

import pickle
import pandas as pd
from pi_paper_utils.spe_calib_utils import pipe_transform_projective

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)

dsst = mhdlab.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst(convert_to_PT=False)
# %%

munged_dir = pjoin(REPO_DIR, 'experiment', 'data', 'munged', datestr)
spe_dir = pjoin(munged_dir, 'spe')

das = {}
names = ['4', '84']

# %%

for name in names:

    fp_in = pjoin(spe_dir, 'PI_topcam_{}.cdf'.format(name))

    ds = xr.load_dataset(fp_in)
    da = ds['counts']

    das[name] = da

#%%


df = pd.DataFrame(index=das.keys())

for name, da in das.items():
    # coord_dict[name]['gatedelay'].plot(label=name)
    coords = da.coords
    gatedelay = coords['gatedelay']

    start = gatedelay[0]
    step = gatedelay[1] - gatedelay[0]
    num = len(gatedelay)

    df.loc[name, 'start'] = start
    df.loc[name, 'step'] = step
    df.loc[name, 'num'] = num

df

#%%

for name, da in das.items():

    es_time = da['estime']
    
    es_time.plot(label=name, marker='o')

plt.legend()


#%%

da = das['84'].copy()

da.mean('estime').mean('gatedelay').plot(robust=True)
#%%
# Timewindow for these files is entire seedramp
da['estime'].plot()

#%%

# downselect to 1% kwt time window
tw = slice(Timestamp('2023-05-24 22:25:55.952985600'), Timestamp('2023-05-24 22:28:12.304348160'), None) 

da_sel = da.sel(estime=tw).sel(gatedelay=slice(780,810))

da_sel = da_sel.mean('estime')

da_sel

#%%

fp_transform = pjoin(DIR_DATA_OUT, 'tform_projective.pkl')
with open(fp_transform, 'rb') as f:
    tform = pickle.load(f)

da_sel_2 = da_sel.to_dataset()
downsel_range = {'x': slice(100, 1024), 'y': slice(420, 620)}
da_sel_tf = pipe_transform_projective(da_sel_2, tform, downsel_range=downsel_range)

da_sel_tf = da_sel_tf.sel(y=slice(-30,30))

#%%

# subtract off offset

offset = da_sel_tf.sel(gatedelay=780).sel(y=slice(25,30)).mean('x').mean('y').item()

offset

da_sel_tf = da_sel_tf - offset

da_sel_tf.attrs['long_name'] = 'Counts'
da_sel_tf.coords['gatedelay'].attrs['units'] = 'ns'

#%%

da_sel_tf.to_netcdf(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))