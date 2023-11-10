"""
Script to munge all relevant lecroy dates. This should probably be incorporated in the PostProcessor repository at some point.
"""


from mhdpy.analysis.standard_import import *

from mhdpy.mws_utils.fileio import get_matching_filepath_mnums, NoMatchingFilesError, load_trc_mnum_simple
from mhdpy.mws_utils.trc import setupfile_timeoffset

datestr = '2023-05-24'

lecroy_setup_dir = r'Z:\HVOF Booth\F\Lecroy Oscope\Setups'
lecroy_wfm_dir = os.path.join(r'Z:\HVOF Booth\F\Lecroy Oscope\Waveforms', datestr)

output_dir = os.path.join('munged',datestr, 'Lecroy')
if not os.path.exists(output_dir): os.makedirs(output_dir)

channel_dict = {
    '1':'i',
    '2':'q',
    '3': 'pd1',
    '4': 'pd2',
}

timeoffset_1 = setupfile_timeoffset(pjoin(lecroy_setup_dir, r'Setup_2023_05_24_14_22_00_592076--00000.lss'))
extra_time_offset =  - np.timedelta64(7,'h')
time_offset = timeoffset_1 + extra_time_offset

# # Setup File to calibrate time difference between lecroy and lab computer. 
# https://stackoverflow.com/questions/312443/how-do-i-split-a-list-into-equally-sized-chunks
def get_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def main(lecroy_wfm_dir, channel_dict, time_offset, output_dir, chunk_size=1000):

    #TODO: Check for regex matching here (i.e. ignore C3 and C4)
    for dirpath, dirnames, filenames in os.walk(lecroy_wfm_dir):
        if len(filenames) > 3:
            print(dirpath)

            subfolder = os.path.relpath(dirpath, lecroy_wfm_dir)

            try:
                mnums = get_matching_filepath_mnums(filenames, channels=channel_dict.keys())
            except NoMatchingFilesError as e:
                print("No Matching files at all, skipping")
                continue

            dss_chunk = []

            for mnum_list in get_chunks(mnums, chunk_size):
                mnum_list = mnum_list[::100]    
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

            out_fn = subfolder.replace('\\','_')
            

            ds.to_netcdf(os.path.join(output_dir, 'ds_{}.cdf'.format(out_fn)))

main(lecroy_wfm_dir, channel_dict, time_offset, output_dir)


