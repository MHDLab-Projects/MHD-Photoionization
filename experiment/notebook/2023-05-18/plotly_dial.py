#%%
from mhdpy.analysis.standard_import import *
create_standard_folders()
DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')

#have time axes show up in PST
plt.rcParams['timezone'] = 'US/Pacific'
import matplotlib.dates as mdates

# %%

fp_dsst = pjoin(DIR_EXPT_PROC_DATA, 'dsst.tdms')
dsst = mhdpy.fileio.TFxr(fp_dsst).as_dsst(convert_to_PT=False)

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws)
# df_exptw = pd.read_csv(fp_expt_tws)

df_exptw['date'] = df_exptw['Start Time'].apply(lambda x: str(x.date()))
df_exptw = df_exptw.set_index('date')

# %%

da = dsst['hvof']['CC_K_massFrac_in']

tw = df_exptw.ct.index_slice('2023-05-18')

da.sel(time=tw).plot()  

#%%


from mhdpy.plot.anim.plotly_dial import gen_movie_da

fp_out = 'output/plotly_dial.mp4'

da_sel = da.sel(time=tw).isel(time=slice(0,-1,100))

gen_movie_da(da_sel, fp_out)


# %%
