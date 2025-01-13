# Note, kaleido needs to be 1.0.0rc0 for this to work.

#%%
from mhdlab.analysis.standard_import import *
from params import duration_lookup, time_downselect_lookup
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdlab.fileio.TFxr(fp_dsst).as_dsst(convert_to_PT=False)

# var = 'CC_K_massFrac_in'
var = 'CC_total_flow_in'
da = dsst['hvof'][var]

from mhdlab.plot.anim.plotly_dial import gen_movie_da

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdlab.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

for date, tw in df_exptw.ct.iterslices():
    da_plot = da.sel(time=tw)

    fp_out = 'output/dial_{}_{}.mp4'.format(var, date)

    time_downselect = time_downselect_lookup['dials']
    if time_downselect:
        da_plot = da_plot.isel(time=slice(0,-1,time_downselect))
    duration = duration_lookup[date]

    gen_movie_da(da_plot, fp_out, tw, duration)
