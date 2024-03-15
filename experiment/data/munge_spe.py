"""
Munging script for princeton instruments camera files. 

The files are grouped by filesize. 
TODO: read file attributes to categorize 
"""


#%%

from mhdpy.analysis.standard_import import *
from mhdpy.fileio.spe import spe2ds_img, _get_gatedelays
from tqdm import tqdm
from collections import defaultdict

#%%

import dotenv; dotenv.load_dotenv()
RAW_DATA_DIR = os.getenv('RAW_DATA_FOLDER') #This will throw error if no .env file with REPO_DIR defined in analysis repo. 


def main(datestr):
    # folder = r'Z:\HVOF Booth\H\{}\PI_TopCam'.format(datestr)
    folder = pjoin(RAW_DATA_DIR, datestr, 'PI_TopCam' )


    output_dir = os.path.join('munged',datestr, 'spe')
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    fns = os.listdir(folder)#[0:10]

    # group all files in folder by their file size on disk, rounded to mega bytes

    fn_grouped = defaultdict(list)

    for fn in fns:
        fp = os.path.join(folder, fn)
        size = os.path.getsize(fp)
        size_mb = round(size/1e6)
        fn_grouped[size_mb].append(fn)

    fn_grouped_count = {k: len(v) for k, v in fn_grouped.items()}

    fn_grouped_count = {k: v for k, v in fn_grouped_count.items() if v > 10}

    #quick hack to skip these as they are repetitive. Can we check that all files are repetitive quickly? Think checking attrs requires loading of whole file. 
    if datestr == '2023-05-18':
        del fn_grouped_count[210] 


    for size in fn_grouped_count:
        print("loading files of size: {}".format(size))

        dss = []

        fns = fn_grouped[size]
        # fns = fns[0:5]

        for fn in tqdm(fns):

            fp = os.path.join(folder, fn)

            #TODO: add the repetitve images? 
            ds = spe2ds_img(fp, roiy=slice(412,612), gatingmode_require='Sequential')


            if ds is None:
                print("Got None for dataset, skipping")
                continue

            dss.append(ds)

        if len(dss) == 0:
            print("Did not get any datasets for file size {}".format(size))
            continue

        # Gate delays should be same sinze file sizes are the same, exact asserts this
        # ds_out = xr.concat(dss,'estime', join='exact')

        # Gate delays can be different if the start time was changed through run
        #TODO: outer should be ok, but maybe interpolate gatedelays before exact join to be more explicit. 
        ds_out = xr.concat(dss, 'estime', join='outer')
        ds_out = ds_out.sortby('estime')


        # #Coarsen to reduce file size. TODO: Remove for final analysis
        # ds_out = ds_out.coarsen(y=4, boundary='trim').mean()
        # ds_out = ds_out.coarsen(x=4, boundary='trim').mean()

        fp_out = os.path.join(output_dir, 'PI_topcam_{}.cdf'.format(size))

        ds_out.to_netcdf(fp_out)


datestrs = ['2023-05-18', '2023-05-24']

for datestr in datestrs:
    main(datestr)

#%%

# Convert calibration images to cdf 

folder_calibration_image = pjoin(RAW_DATA_DIR, '2023-05-18', 'PI_TopCam', 'Manual_Image')
output_dir = pjoin('munged','2023-05-18', 'spe', 'calibration')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for fn in os.listdir(folder_calibration_image):
    fp = pjoin(folder_calibration_image, fn)
    ds = spe2ds_img(fp, gatingmode_require='Repetitive')

    fn_base = os.path.splitext(fn)[0]
    fp_out = pjoin(output_dir, "{}.cdf".format(fn_base))

    ds.to_netcdf(fp_out)
