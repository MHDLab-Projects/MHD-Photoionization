#%%

from mhdpy.analysis.standard_import import *

fp = pjoin(REPO_DIR, 'experiment', 'data', 'proc_data', 'ds_absem.cdf')

ds = xr.load_dataset(fp)

ds = ds.absem.calc_alpha()

ds = ds.set_index(acq=['time','mp']).unstack('acq').rename(time='acq_time')

ds

#%%

fp_expt_tws = pjoin(REPO_DIR, 'experiment', 'metadata', 'ct_experiment.csv')
df_exptw = mhdpy.fileio.load_df_cuttimes(fp_expt_tws).set_index('Event')

tw = df_exptw.ct.index_slice('2023-05-18')
ds = ds.sel(acq_time=tw)

#%%

import matplotlib.pyplot as plt
import numpy as np
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

da_plot = ds['alpha']#.isel(acq_time=slice(0,-1,10))
da_plot = da_plot.sel(wavelength=slice(760,780))

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

fig, axes = plt.subplots(1,2, figsize = (10,5))

axes[0].set_title('Barrel')
axes[1].set_title('Seed Port')

def make_frame(t):
    nearest_time = tw_grid[int(t*fps)]
    da_sel = da_plot.sel(acq_time=nearest_time, method='nearest')

    da_sel = da_sel.dropna('mp','all')

    remaining_mp = da_sel.coords['mp'].item()
    if remaining_mp == 'barrel':
        ax = axes[0]
        title= 'Barrel'
    elif remaining_mp == 'mw_horns':
        ax = axes[1]
        title = 'MW Horns'
    else:
        raise ValueError

    ax.clear()
    ax.set_ylim(-1,2)
    da_sel.sel(mp=remaining_mp).plot(ax=ax)
    ax.set_title(title)
    fig.suptitle(str(nearest_time))


    fig.tight_layout()

    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=duration)
animation.write_videofile('output/absem_movie.mp4', fps=fps)



# %%
