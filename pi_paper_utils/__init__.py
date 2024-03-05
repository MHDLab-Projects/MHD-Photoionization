from . import fileio

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

# constant properties

from pint import Quantity
# Determined in experiment/analysis/various/power_calibration.py
LASER_POWER = Quantity(6.98836, 'mJ')

# Determined in experiment/analysis/mws/laser_profile.py
LASER_AREA = Quantity(40.77, 'mm^2')

MOTOR_OFFSET = Quantity(1 + 3/16, 'inches')