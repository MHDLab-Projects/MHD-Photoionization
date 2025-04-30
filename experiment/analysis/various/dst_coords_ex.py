"""
Examine process of binning data into coordinates. 

The equivalence ratio coordinates appear to be off at 0.65. Examining why
Appears to be because including all experiment data, not just the test case sequences.
When the syringe shuts off, the torch jumps to 0.68 equivalence ratio.
Going to downselect to test case sequences and rebin the data.
"""
#%%

from mhdlab.analysis.standard_import import *
create_standard_folders()
import re
from collections import defaultdict

from mhdlab.coords import gen_coords_to_assign_1

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

# Absorption emsssion and lecroy
ds_aas = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'ds_aas.cdf'))
ds_aas = ds_aas.set_index(acq=['time','mp']).unstack('acq')
ds_aas = ds_aas.rename(time='acq_time')
ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'ds_lecroy.cdf'))


dsst = mhdlab.fileio.TFxr(pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')).as_dsst(convert_to_PT=False)
# Generate binned coordinates to assign along time dimension

dst = xr.merge([
    dsst['hvof']['CC_K_massFrac_in'],
    dsst['motor']['Motor C Relative'],
    dsst['hvof']['CC_equivalenceRatio']
]).rename({
    'CC_K_massFrac_in': 'kwt',
    'Motor C Relative': 'motor',
    'CC_equivalenceRatio': 'phi'
})

coords_to_assign = gen_coords_to_assign_1(
        dsst,
        kwt_bins =  [-0.01, 0.01, 0.08, 0.15, 0.25 , 0.4, 0.75, 1.5, 2.5, 5],
        motor_bins = np.linspace(-12.5,237.5,11),
        phi_bins = [0.5,0.7,0.9,1.1,1.3,1.5]
    )
dss = [da.to_dataset(name=name) for name, da in coords_to_assign.items()]
ds_coords = xr.merge(dss)

# %%

fp_cuttimes = pjoin(REPO_DIR,'experiment','metadata', 'ct_sequence.csv')
df_ct = mhdlab.fileio.load_df_cuttimes(fp_cuttimes)

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws)
# df_exptw = pd.read_csv(fp_expt_tws)

df_ct['date'] = df_ct['Start Time'].apply(lambda x: x.date())
df_exptw['date'] = df_exptw['Start Time'].apply(lambda x: x.date())
df_exptw = df_exptw.set_index('date')

df_exptw


# %%

def dst_coord_plot(dst, dst_coords, tw):

    fig, axes = plt.subplots(3)

    dst['phi'].sel(time=tw).plot(ax=axes[0])
    dst_coords['phi'].sel(time=tw).plot(color='r', ax=axes[0], alpha=0.5)

    (dst['kwt']*100).sel(time=tw).plot(ax=axes[1])
    dst_coords['kwt'].sel(time=tw).plot(color='r', ax=axes[1], alpha=0.5)

    axes[1].set_ylim(0,1)

    dst['motor'].sel(time=tw).plot(ax=axes[2])

    dst_coords['motor'].sel(time=tw).plot(color='r', ax=axes[2], alpha=0.5)

    axes[0].set_ylim(0.5,0.9)

row = df_exptw.loc[pd.to_datetime('2023-05-24').date()]
tw = slice(row['Start Time'], row['Stop Time'])

dst_coord_plot(dst, ds_coords, tw)




#%%[markdown]

# Downselect to test case sequences

# # still doesn't work becaue when syringe shutds off the torch jumpts to 0.68 equivalence ratio


#%%

from mhdlab.coords.ct import downselect_acq_time
ds_hvof_tcs = downselect_acq_time(dsst['hvof'], df_ct, 'time')
ds_motor = downselect_acq_time(dsst['motor'], df_ct, 'time')

dsst_tcs = {
    'hvof': ds_hvof_tcs,
    'motor': ds_motor,
    'filterwheel': dsst['filterwheel']
}

dst_tcs = xr.merge([
    dsst_tcs['hvof']['CC_K_massFrac_in'],
    dsst_tcs['motor']['Motor C Relative'],
    dsst_tcs['hvof']['CC_equivalenceRatio']
]).rename({
    'CC_K_massFrac_in': 'kwt',
    'Motor C Relative': 'motor',
    'CC_equivalenceRatio': 'phi'
})



#%%

coords_to_assign = gen_coords_to_assign_1(
    dsst_tcs,
    kwt_bins =  [-0.01, 0.01, 0.08, 0.15, 0.25 , 0.4, 0.75, 1.5, 2.5, 5],
    motor_bins = np.linspace(-12.5,237.5,11),
    phi_bins = [0.5,0.7,0.9,1.1,1.3,1.5]
    )

dss = [da.to_dataset(name=name) for name, da in coords_to_assign.items()]
dst_coords_tcs = xr.merge(dss)
# %%
row = df_exptw.loc[pd.to_datetime('2023-05-24').date()]
tw = slice(row['Start Time'], row['Stop Time'])

dst_coord_plot(dst_tcs, dst_coords_tcs, tw)

#%%

set(dst_coords_tcs['phi'].dropna('time').values)
# %%[markdown]

# change bins

#%%

nominal_phi = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6]
bin_width = 0.05

phi_bins = [ [x-bin_width, x+bin_width] for x in nominal_phi]
phi_bins = [x for sublist in phi_bins for x in sublist]
phi_bins

#%%

coords_to_assign = gen_coords_to_assign_1(
    dsst,
    kwt_bins =  [-0.01, 0.01, 0.08, 0.15, 0.25 , 0.4, 0.75, 1.5, 2.5, 5],
    motor_bins = np.linspace(-12.5,237.5,11),
    phi_bins = phi_bins
    # phi_bins = [0.5,0.7,1,1.5]
    )

dss = [da.to_dataset(name=name) for name, da in coords_to_assign.items()]
ds_coords_tcs = xr.merge(dss)
# %%
row = df_exptw.loc[pd.to_datetime('2023-05-24').date()]
tw = slice(row['Start Time'], row['Stop Time'])

dst_coord_plot(dst, ds_coords_tcs, tw)


# %%
set(ds_coords_tcs['phi'].dropna('time').values)