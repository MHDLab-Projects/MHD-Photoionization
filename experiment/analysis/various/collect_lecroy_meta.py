#%%

#%%

from mhdpy.analysis.standard_import import *
import pi_paper_utils as ppu
import json

datestrs = ['2023-04-07', '2023-05-12', '2023-05-18', '2023-05-24']

for datestr in datestrs:

    with open(pjoin(REPO_DIR, 'experiment', 'metadata', 'settings.json')) as f:
        settings = json.load(f)[datestr]

    process_tcs = settings['mwt_process_tcs']

    # Start processing

    data_folder = os.path.join(REPO_DIR, 'experiment','data','munged',datestr)

    lecroy_munged_folder = pjoin(data_folder, 'Lecroy')

    process_fns = ['ds_{}.cdf'.format(tc) for tc in process_tcs]
    input_fps = [pjoin(lecroy_munged_folder, fn) for fn in process_fns] 


    # input_fps = ['output/ds_516.cdf']
    # fp = input_fps[0]

    if datestr == '2023-05-24':
        vars = ['i','q','pd1','pd2']
    else:
        vars = ['i','q']

    for var in vars:
        ss = []
        for fp in input_fps:
            if not os.path.exists(fp):
                print("Couldn't find: {}".format(fp))
                continue
            ds = xr.load_dataset(fp)

            fn = os.path.basename(fp)
            fn = fn.split('.')[0]

            s = pd.Series(ds[var].attrs)
            s.name = fn
            ss.append(s)

        df_out = pd.concat(ss, axis=1)
        df_out

        output_dir = pjoin(DIR_DATA_OUT, 'lecroy_meta')

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)


        df_out.to_csv(pjoin(output_dir, 'lecroy_meta_{}_{}.csv'.format(datestr,var)))
