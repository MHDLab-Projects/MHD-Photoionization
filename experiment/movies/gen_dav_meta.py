"""
Generates a dataframe of metadata for all .dav files in the FLIR folder. Saved to CSV to be loaded in when generating timelapses.
"""
#%%
import os
import pandas as pd
import dotenv; dotenv.load_dotenv(override=True)

from mhdlab.video import generate_df_dav_meta

#%%

dfs = []
dvr_folder = os.getenv('RAW_FLIR_FOLDER')
date_range = pd.date_range('2023-04-06', '2023-05-25').date

df_dav_meta = generate_df_dav_meta(dvr_folder, date_range=date_range)

df_dav_meta


#%%
# Quick fix for this specific file, which seem t obe the only that throws above error with no start time...

df_dav_meta.loc['DNR708_ch5_main_20230524104412_20230524110000.dav', 'ffprobe_start'] = df_dav_meta.loc['DNR708_ch6_main_20230524104412_20230524110000.dav', 'ffprobe_start']
df_dav_meta.loc['DNR708_ch5_main_20230524104412_20230524110000.dav', 'ffprobe_stop'] = df_dav_meta.loc['DNR708_ch6_main_20230524104412_20230524110000.dav', 'ffprobe_stop']

#%%
df_dav_meta = df_dav_meta.dropna(how='any')

fp_out = os.path.join('output', 'flir', 'df_dav_meta.csv')

if not os.path.exists(os.path.dirname(fp_out)):
    os.makedirs(os.path.dirname(fp_out))

df_dav_meta.to_csv(fp_out)