# %% [markdown]
#  # 2022-10-05 Photoionization experiment with motor position dependence and IQ mixer
#  Laser and mixer setup were moving together.

# %%

from mhdpy.analysis.standard_import import *

datestr = '2023-05-24'
datestr = datestr
data_folder = pjoin(REPO_DIR, 'experiment','data','munged', datestr)
DIR_PROC_DATA = pjoin(REPO_DIR, 'experiment','data', 'proc_data', 'lecroy')
lecroy_setup_folder = r'Z:\Lecroy Oscope\Setups'
lecroy_raw_folder = os.path.join(r'Z:\Lecroy Oscope\Waveforms', datestr)

channel_dict = {
    '1':'i',
    '2':'q',
    '3': 'pd1',
    '4': 'pd2',
}

# %%

# # Setup File to calibrate time difference between lecroy and lab computer. 
from mhdpy.mws_utils.trc import setupfile_timeoffset

timeoffset_1 = setupfile_timeoffset(pjoin(lecroy_setup_folder, r'Setup_2023_05_24_14_22_00_592076--00000.lss'))

extra_time_offset =  - np.timedelta64(7,'h')
time_offset = timeoffset_1 + extra_time_offset

time_offset

#%%

from mhdpy.mws_utils.fileio import get_matching_filepath_mnums, NoMatchingFilesError, load_trc_mnum_simple

filepaths = [pjoin(lecroy_raw_folder,  fn) for fn in os.listdir(lecroy_raw_folder) if os.path.isfile(pjoin(lecroy_raw_folder, fn))]

# https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

chunk_size = 1000

#TODO: Check for regex matching here (i.e. ignore C3 and C4)
for dirpath, dirnames, filenames in os.walk(lecroy_raw_folder):
    if len(filenames) > 3:
        print(dirpath)

        subfolder = os.path.relpath(dirpath, lecroy_raw_folder)

        try:
            mnums = get_matching_filepath_mnums(filenames, channels=['1','2','3','4'])
        except NoMatchingFilesError as e:
            print("No Matching files at all, skipping")
            continue

        dss_chunk = []

        for mnum_list in chunks(mnums, chunk_size):

            # mnum_list = mnum_list[::100]

            try:
                ds = load_trc_mnum_simple(dirpath, mnum_list,time_offset=time_offset)
            except NoMatchingFilesError as e:
                print("No Matching files")
            except AttributeError as e:
                print(e)
            except FileNotFoundError as e:
                print(e)
            else: 

                ds = ds.coarsen(time=1000, boundary='trim').mean()
                ds = ds.rename(channel_dict)

                dss_chunk.append(ds)


        for i in range(len(dss_chunk)):
            dss_chunk[i] = dss_chunk[i].assign_coords(time=dss_chunk[0].coords['time'])

        ds = xr.concat(dss_chunk, 'acq_time')
        output_folder = os.path.join(data_folder, 'Lecroy')

        out_fn = subfolder.replace('\\','_')

        if not os.path.exists(output_folder): os.makedirs(output_folder)
        ds.to_netcdf(os.path.join(output_folder, 'ds_{}.cdf'.format(out_fn)))

