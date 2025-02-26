#%%
from mhdlab.analysis.standard_import import *
from mhdlab.fileio.spe import spe2ds_img, _get_gatedelays, spe2ds_spect
from tqdm import tqdm
from collections import defaultdict

create_standard_folders()
# %%

import dotenv; load_dotenv()
RAW_DATA_DIR = os.getenv('RAW_DATA_FOLDER') #This will throw error if no .env file with REPO_DIR defined in analysis repo. 
folder = pjoin(RAW_DATA_DIR, '2018-11-20', 'PIMAX_2')

# folder = pjoin(RAW_DATA_DIR, '2017-10-20', 'LightfieldRaw')

fns = os.listdir(folder)

dss = []

for fn in fns:

    fp = pjoin(folder, fn)

    ds = spe2ds_img(fp)

    dss.append(ds)

ds = xr.concat(dss, 'estime')
ds = ds.sortby('estime')

#%%

# A couple times in here from 2018-12-04? Just dropping. 
tw = slice(Timestamp('2018-11-19'),Timestamp('2018-11-22'))
ds = ds.sel(estime=tw)


#%%


output_dir = pjoin('munged','2018-11-20')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ds.to_netcdf(pjoin(output_dir, 'ds_beam_timing.cdf'))


# %%
