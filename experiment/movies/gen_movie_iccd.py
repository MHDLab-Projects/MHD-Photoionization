#%%
from mhdlab.analysis.standard_import import *
import pi_paper_utils as ppu
import matplotlib.ticker as ticker
import numpy as np
import matplotlib.pyplot as plt
from moviepy import VideoClip
from mhdlab.plot.anim import mplfig_to_npimage

create_standard_folders()

# update the dpi to 300 for final figures (see note in mpl.style)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 16

# da_sel_tf.to_netcdf(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))

DIR_DATA_OUT = pjoin(REPO_DIR, 'final', 'dataset', 'topcam', 'output')
da_sel_tf = xr.load_dataarray(pjoin(DIR_DATA_OUT, 'spe_20230524_84_tfm.cdf'))


#%%
# da_sel_tf.mean('gatedelay').plot()


# %%

da_sel_tf

def gen_movie_iccd(da_sel_tf, fp_out, duration):
    """
    Generate a moviepy video of iccd data going through the gate delays.
    """
    da_sel_tf = da_sel_tf.transpose('gatedelay', 'x', 'y')
    assert list(da_sel_tf.dims) == ['gatedelay', 'x', 'y'], "DataArray must have dims ('gatedelay', 'x', 'y')"

    ts = da_sel_tf.coords['gatedelay'].values
    frames = len(ts)

    fps = frames / duration


    def make_frame(t):

        fig, ax = plt.subplots(figsize=(6, 6))
        nearest_time = ts[int(t * fps)]
        da_sel = da_sel_tf.sel(gatedelay=nearest_time, method='nearest')

        da_sel.plot(ax=ax, vmin=0, vmax=500)
        ax.set_title(f'Gate Delay: {nearest_time} ns')

        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')

        fig.tight_layout()

        return mplfig_to_npimage(fig)

    animation = VideoClip(make_frame, duration=duration)
    animation.write_videofile(fp_out, fps=fps)

# Example usage
fp_out = pjoin('output', 'iccd_movie_2025_05_24.mp4')
gen_movie_iccd(da_sel_tf, fp_out, duration=5)
# %%
