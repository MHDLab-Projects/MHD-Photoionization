#%%
import os
import pandas as pd
from pandas import Timestamp

from mhdlab.fileio.ct import load_df_cuttimes

from mhdlab.video import load_df_dav_meta, ffmpeg_davs_to_mp4

import dotenv; dotenv.load_dotenv(override=True)
#%%

#Run gen_dav_meta.py first
fp_df_dav_meta = os.path.join('output', 'flir', 'df_dav_meta.csv')

fp_timewindows = os.path.join(os.getenv('repo_dir'), 'experiment', 'metadata', 'ct_experiment.csv')
df_ct = load_df_cuttimes(fp_timewindows).set_index('Event')

fp_timewindows_extra = os.path.join('ct_movie.csv')
df_ct_add = load_df_cuttimes(fp_timewindows_extra).set_index('Event')

df_ct = pd.concat([df_ct, df_ct_add])

# idxs = ['2023-04-07', '2023-05-12', '2023-05-18', '2023-05-24']
idxs = ['2023-05-18', '2023-05-18_setup']
channels = ['5','6','7','8']

for idx in idxs:

    #TODO: FLIR timestamps are in PT. See issue #13
    ts_start = Timestamp(df_ct.loc[idx]['Start Time'], tz='UTC').tz_convert('US/Pacific').tz_localize(None)
    ts_stop = Timestamp(df_ct.loc[idx]['Stop Time'], tz='UTC').tz_convert('US/Pacific').tz_localize(None)

    output_timewindow = slice(
        ts_start,
        ts_stop
        )

    if 'setup' in idx:
        date = idx.split('_')[0]
    else:
        date = idx

    date = Timestamp(date).date()

    folder_in = os.path.join(os.getenv('RAW_FLIR_FOLDER'), date.strftime("%Y"))

    folder_out = os.path.join('output', 'flir', idx)
    if not os.path.exists(folder_out): os.makedirs(folder_out)

    dfs_meta = []

    for channel in channels:

        df_dav_meta = load_df_dav_meta(fp_df_dav_meta, output_timewindow=output_timewindow, channel_sel=channel)

        #TODO: how to consolidate rows to one per outputfile, describing start/end timestamp for lining up grid video
        # maybe integrate with the for loop in the ffmpeg_davs_to_mp4 fn 
        dfs_meta.append(df_dav_meta)

        fn_out = 'Timelapse_{}_ch{}.mp4'.format(idx, channel)
        fp_out = os.path.join(folder_out, fn_out)

        if os.path.exists(fp_out):
            print(f"{fp_out} already exists, skipping")
            continue

        if len(df_dav_meta) == 0:
            print(f"Got zero length for {fp_out}")
            continue

        ffmpeg_davs_to_mp4(df_dav_meta, fp_out, fps=30, crf=20, framestep=30)


    df_dav_meta_combined = pd.concat(dfs_meta)

    fn_out = 'Timelapse_{}_meta.csv'.format(date.strftime("%Y-%m-%d"))
    fp_out = os.path.join(folder_out, fn_out)

    df_dav_meta_combined.to_csv(fp_out)