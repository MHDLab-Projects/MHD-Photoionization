# constant properties

from pint import Quantity
# Determined in experiment/analysis/various/power_calibration.py
LASER_POWER = Quantity(6.98836, 'mJ')

# Determined in experiment/analysis/mws/laser_profile.py
LASER_AREA = Quantity(40.77, 'mm^2')

MOTOR_OFFSET = Quantity(1 + 3/16, 'inches')

## final\analysis\mws_block_expt.py
MAG_MWS_BLOCK = Quantity(0.006245581698287483, 'V')
# MAG_MWS_BLOCK = Quantity(0, 'V')

