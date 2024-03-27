#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
import pi_paper_utils as ppu

raw_folder = r"C:\Users\ASPITARL\OneDrive - National Energy Technology Laboratory\Documents - MHD Lab Team\Data Share\MHD Lab\Instrumentation Room\2024-03-27"


from mhdpy.fileio.lecroy import get_matching_filepath_mnums, load_trc_mnum_simple, setupfile_timeoffset

filenames = os.listdir(raw_folder)

channel_dict= {'1':'i','2':'q'}
channel_list = list(channel_dict.keys())

middle_regex = '(?:IQ_)?\w+'
channel_regex_str = "".join(channel_list)

full_regex = 'C([{}])--_({})_--(\d+).trc'.format(channel_regex_str, middle_regex)



fname_middle, mnums = get_matching_filepath_mnums(filenames, full_regex)
ds = load_trc_mnum_simple(raw_folder, fname_middle, mnums)

ds = ds.rename(channel_dict)[['i','q']]


ds['i'].attrs['units'] = 'V'
ds['q'].attrs['units'] = 'V'
ds.coords['time'].attrs['units'] = 's'
ds = ds.pint.quantify()
ds['time'] = ds['time'].pint.to('microsecond')

#There doesn't appear to be any sort of sorting going on in load_trc function, this should maintain order
ds = ds.assign_coords(mnum = ('acq_time', [int(m) for m in mnums]))
ds = ds.set_index(acq_time='mnum').rename(acq_time='mnum')

ds = ds.mws.calc_mag_phase()

ds=  ds.mws.calc_time_stats()

ds 

#%%
ds['mag']


#%%

mnum_dict = {
    5: 'nothing_1', 
    6: 'metal_middle_1',
    11: 'terminator',
    21: 'nothing_2',
    22: 'metal_middle_2',
}

ds2 = ds.sel(mnum=list(mnum_dict.keys()))

ds2= ds2.assign_coords(mnum = [mnum_dict[m] for m in mnum_dict])

ds2 = ds2.rename(mnum='tc')

ds2

#%%

ds2['mag_pp'].plot()

plt.xticks(rotation=90)


#%%
mag_block = ds2['mag_pp'].sel(tc='terminator')

mag_block

#%%

mag_nothing = ds2['mag_pp'].sel(tc='nothing_2')

mag_nothing


#%%

T_block = mag_block/mag_nothing

AS_block = 1 - T_block

AS_block

# ds.mws.calc_time_stats()['mag_pp'].plot()
# %%
