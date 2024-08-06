from . import fileio
from .constants import *

import matplotlib.pyplot as plt
import scienceplots

# Apply styles elsewhere after importing package
# plt.style.use(['science', 'ieee', 'pi_paper_utils.mystyle'])

# Apply styles here. need to pass in a Path object to plt.style.use
from pathlib import Path
style_path = (Path(__file__).parent / 'mystyle.mplstyle').as_posix()

# see venv/lib/python3.11/site-packages/matplotlib/mpl-data/stylelib for available styles

plt.style.use(['science', style_path])
# plt.style.use(['science', 'ieee', style_path])
# plt.style.use(['science', 'seaborn-v0_8-deep', style_path])

import numpy as np
import xarray as xr

ndf_dict = {
1: 0,
2: 0.4,
3: 0.7,
4: 1, 
5: 1.5,
6: np.inf
}

ndf_dict = {key: 10**(-ndf_dict[key]) for key in ndf_dict}

def convert_fw_pos_relpower(da_fw_pos):
    """
    Converts dsst['filterwheel']['Filter Position'] to relative power
    Can later be converted to absolute power using the laser calibration
    """
    filter_coord = xr.DataArray(da_fw_pos.to_series().map(ndf_dict))
    filter_coord.attrs = dict(long_name='Rel. Power')

    return filter_coord 