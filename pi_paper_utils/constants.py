# constant properties

from pint import Quantity
# Determined in experiment/analysis/various/power_calibration.py
LASER_POWER = Quantity(6.98836, 'mJ')

# Determined in experiment/analysis/mwt/laser_profile.py
LASER_AREA = Quantity(40.77, 'mm^2')

MOTOR_OFFSET = Quantity(1 + 3/16, 'inches')

## final\analysis\mwt_block_expt.py
MAG_MWT_BLOCK = Quantity(0.006245581698287483, 'V')
# MAG_MWT_BLOCK = Quantity(0, 'V')

AES_BARREL_OFFSET = Quantity(1/8, 'inches')

# Taken from geometry.toml #TODO: read this file?
CFD_EXIT_OFFSET = Quantity(20.911, 'cm')

# species name used in fitting and plots
CFD_K_SPECIES_NAME = 'K'
CFD_KOH_SPECIES_NAME = 'KOH'
CFD_Kp_SPECIES_NAME = 'Kp'
CFD_allK_SPECIES_NAME = 'all_K'
CFD_OH_SPECIES_NAME = 'OH'

# CFD_K_SPECIES_NAME = 'Yeq_K'
# CFD_KOH_SPECIES_NAME = 'Yeq_KOH'
# CFD_Kp_SPECIES_NAME = 'Yeq_K+'
# CFD_allK_SPECIES_NAME = 'all_K_Yeq'
# CFD_OH_SPECIES_NAME = 'Yeq_OH'
