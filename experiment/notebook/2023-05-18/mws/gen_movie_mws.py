#%%

from mhdpy.analysis.standard_import import *

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_lecroy.cdf')

ds = xr.load_dataset(fp)
ds = ds[['i','q']]

ds = ds.mws.calc_mag_phase()

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

tw = df_exptw.ct.index_slice('2023-05-18')
ds = ds.sel(acq_time=tw)

#%%

import matplotlib.pyplot as plt
import numpy as np
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

da_plot = ds['mag'].isel(acq_time=slice(0,-1,10))
da_plot = da_plot.sel(time=slice(-1,30))

ts = da_plot.coords['acq_time'].values
frames = len(ts)

tw_grid = np.linspace(
    pd.to_datetime(tw.start).timestamp(),
    pd.to_datetime(tw.stop).timestamp(),
    frames
)

tw_grid = pd.to_datetime(tw_grid, unit='s')

duration = 132
fps = frames/duration

fig, ax = plt.subplots(1,1, figsize = (10,5))

def make_frame(t):
    nearest_time = tw_grid[int(t*fps)]
    da_sel = da_plot.sel(acq_time=nearest_time, method='nearest')

    ax.clear()
    da_sel.plot(ax=ax)

    ax.set_ylim(0,0.1)

    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('output/mws_movie.mp4', fps=fps)


