"""
Script to munge all relevant lecroy dates. This should probably be incorporated in the PostProcessor repository at some point.
"""


from mhdpy.analysis.standard_import import *
create_standard_folders()

from mhdpy.fileio.lecroy import get_matching_filepath_mnums, load_trc_mnum_simple, setupfile_timeoffset
from mhdpy.fileio.lecroy import NoMatchingFilesError, NotAllFilenameBaseMatchingError

# # Setup File to calibrate time difference between lecroy and lab computer. 
# https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
def get_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]



def main(
        lecroy_wfm_dir, 
        channel_dict, 
        time_offset, 
        output_dir, 
        chunk_size=1000,
        coarsen_amount=1000,
        skip_existing=False,
        subfolder_downselect=None
        ):

    channel_list = list(channel_dict.keys())

    middle_regex = '(?:IQ_)?\w+'
    channel_regex_str = "".join(channel_list)

    full_regex = 'C([{}])--_({})_--(\d+).trc'.format(channel_regex_str, middle_regex)


    #TODO: Check for regex matching here (i.e. ignore C3 and C4)
    for dirpath, dirnames, filenames in os.walk(lecroy_wfm_dir):
        if len(filenames) > 3:
            print(dirpath)

            subfolder = os.path.relpath(dirpath, lecroy_wfm_dir)


            # out_fn = subfolder.replace('\\','_') # Old windows method, remove when tested on windows
            out_fn = subfolder.replace(os.path.sep, '_')
            fp_out = os.path.join(output_dir, 'ds_{}.cdf'.format(out_fn))

            if skip_existing:
                if os.path.exists(fp_out):
                    print("{} exists, skipping".format(fp_out))
                    continue

            if subfolder_downselect:
                if subfolder != subfolder_downselect:
                    continue
            

            try:
                fname_middle, mnums = get_matching_filepath_mnums(filenames, full_regex)
            except NoMatchingFilesError as e:
                print("No Matching files at all, skipping")
                continue
            except NotAllFilenameBaseMatchingError as e:
                print("Middle portion of files did not all match in folder: {}".format(dirpath))
                continue

            dss_chunk = []

            for mnum_list in get_chunks(mnums, chunk_size):
                # mnum_list = mnum_list[::500]    
                try:
                    ds = load_trc_mnum_simple(dirpath, fname_middle, mnum_list, channels=channel_list, time_offset=time_offset)
                except NoMatchingFilesError as e:
                    print("No Matching files")
                except AttributeError as e:
                    print(e)
                except FileNotFoundError as e:
                    print(e)
                else: 

                    #TODO: make kwarg once method is dialed in
                    #TODO: move down where first time array is assigned?
                    if coarsen_amount:
                        # ds = ds.coarsen(time=coarsen_amount, boundary='trim').mean()

                        t1 = 0.7e-6
                        t2 = 1.5e-6

                        ds_1 = ds.sel(time=slice(-1, t1))
                        ds_1 = ds_1.coarsen(time=1000, boundary='trim').mean(keep_attrs=True) 
                        ds_2 = ds.sel(time=slice(t1,t2))
                        ds_2 = ds_2.coarsen(time=10, boundary='trim').mean(keep_attrs=True) 
                        ds_3 = ds.sel(time=slice(t2, 1))
                        ds_3 = ds_3.coarsen(time=1000, boundary='trim').mean(keep_attrs=True) 

                        ds = xr.concat([ds_1, ds_2, ds_3], 'time')

                    #TODO: This happens for one folder on 2023-05-12...why? move this into load_trc_mnum_simple?
                    if not all([name in ds.data_vars for name in channel_dict]):
                        print("Did not get all channels in returned dataset...skipping")
                        continue
                    ds = ds.rename(channel_dict)

                    dss_chunk.append(ds)

            if len(dss_chunk):

                for i in range(len(dss_chunk)):
                    dss_chunk[i] = dss_chunk[i].assign_coords(time=dss_chunk[0].coords['time'])

                ds = xr.concat(dss_chunk, 'acq_time')

                ds.to_netcdf(fp_out)

process_date_dict = {
    '2023-04-07': {
        'time_offset_file': 'Setup_2023_04_07_13_48_34_810411--00000.lss',
        'channel_dict': {'1':'i','2':'q'},
    },
    '2023-05-12': {
        'time_offset_file': 'Setup_2023_05_12_13_22_06_740182--00000.lss',
        'channel_dict': {'1':'i','2':'q'},
    },
    '2023-05-18': {
        'time_offset_file': 'Setup_2023_05_18_15_22_31_204405--00000.lss',
        'channel_dict': {'1':'i','2':'q'},
    },
    '2023-05-24': {
        'time_offset_file': 'Setup_2023_05_24_14_22_00_592076--00000.lss',
        'channel_dict': {'1':'i','2':'q','3':'pd1','4':'pd2'},
    },
}

# No resample settings, huge filesize
# output_base_dir = 'munged_100'
# skip_existing=False
# coarsen_amount=100
# subfolder_downselect='pos_1'

# # Resampel settings 
output_base_dir = 'munged'
skip_existing=True
coarsen_amount=100 #TODO: coarsen amount not being used. 
subfolder_downselect=None
chunk_size = 200 

import dotenv; dotenv.load_dotenv()
LECROY_DIR = os.getenv('LECROY_RAW_FOLDER') #This will throw error if no .env file with REPO_DIR defined in analysis repo. 

lecroy_setup_dir = pjoin(LECROY_DIR, 'Setups')
lecroy_wfm_dir = pjoin(LECROY_DIR, 'Waveforms')


from multiprocessing import Pool

def process_date(datestr, date_dict):
    print("Processing Lecroy Data {}".format(datestr))

    lecroy_wfm_date_dir = pjoin(lecroy_wfm_dir, datestr)

    output_dir = os.path.join(output_base_dir, datestr, 'Lecroy')
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    channel_dict = date_dict['channel_dict']
    time_offset = date_dict['time_offset_file']

    timeoffset_1 = setupfile_timeoffset(pjoin(lecroy_setup_dir, time_offset))
    extra_time_offset =  - np.timedelta64(7,'h')
    time_offset = timeoffset_1 + extra_time_offset

    print("Time offset for date {} is {}".format(datestr, time_offset))

    main(
        lecroy_wfm_date_dir, 
        channel_dict, 
        time_offset, 
        output_dir,
        skip_existing=skip_existing,
        coarsen_amount=coarsen_amount,
        subfolder_downselect=subfolder_downselect,
        chunk_size=chunk_size
    )

for datestr, date_dict in process_date_dict.items():
    process_date(datestr, date_dict)

# #TODO: multiprocess hangs, memory error?
# #TODO: multiprocess all folders in a date (>4) vs all dates (4)? Believe can set max threads. 
# #TODO: tried to ask copilot to use tqdm.contrib.concurrent.process_map,  to stop the output progress bars from being overwritten, but it didn't work.``
# # Create a pool of workers
# with Pool() as p:
#     # Map the function to the data
#     p.starmap(process_date, process_date_dict.items())