#%%

from mhdpy.analysis.standard_import import *
from mhdpy.fileio.spe import spe2ds_img, _get_gatedelays

datestr = '2023-05-18'

folder = r'Z:\HVOF Booth\H\2023-05-18\PI_TopCam'

# %%

output_dir = os.path.join('munged',datestr, 'spe')
if not os.path.exists(output_dir): os.makedirs(output_dir)

#%%

fns = os.listdir(folder)#[0:10]

# for fn in fns:
#     _get_gatedelays()

fn_groups = {}

# group all files in folder by their file size on disk, rounded to mega bytes

from collections import defaultdict

fn_grouped = defaultdict(list)

for fn in fns:
    fp = os.path.join(folder, fn)
    size = os.path.getsize(fp)
    size_mb = round(size/1e6)
    fn_grouped[size_mb].append(fn)

fn_grouped

#%%
fn_grouped_count = {k: len(v) for k, v in fn_grouped.items()}

fn_grouped_count = {k: v for k, v in fn_grouped_count.items() if v > 10}

fn_grouped_count

#%%



for size in fn_grouped_count:

    dss = []

    fns = fn_grouped[size][0:5]

    for fn in fns:

        fp = os.path.join(folder, fn)

        ds = spe2ds_img(fp)

        dss.append(ds)

    # Gate delays should be same sinze file sizes are the same, exact asserts this
    ds_out = xr.concat(dss,'estime', join='exact')

    fp_out = os.path.join(output_dir, 'PI_topcam_{}.cdf'.format(size))

    ds_out.to_netcdf(fp_out)


    

#%%

# not any faster and shows xml errors

# import spe_loader as sl

# fps = [pjoin(folder, fn) for fn in fns]

# spe_file = sl.load_from_files(fps)