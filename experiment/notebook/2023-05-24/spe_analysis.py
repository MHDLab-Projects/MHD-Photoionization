#%%

from mhdpy.analysis.standard_import import *
create_standard_folders()

datestr = '2023-05-24'
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')

dsst = mhdpy.fileio.TFxr(pjoin(data_folder, 'Processed_Data.tdms')).as_dsst()
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

import pandas as pd

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

da = das['4']

# da.mean('estime').plot(col='gatedelay', col_wrap=4)
da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da.mean('estime').plot(col='gatedelay')


#%%

da = das['84']

da.mean('estime').mean('gatedelay').plot(robust=True)

#%%

da_sel = da.mean('estime').sel(x=slice(650,780), y=slice(450,575))

da_sel.mean('gatedelay').plot(robust=True)

#%%

da_sel.mean(['x', 'y']).plot(marker='o')