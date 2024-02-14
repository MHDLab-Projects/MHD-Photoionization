#%%
from mhdpy.analysis.standard_import import *

from mhdpy.pyvista_utils import CFDDatasetAccessor


def test_convert_rho():
    #TODO: create a file like this from a vtk, and then test that it loads correctly. Move tests to mhdpy
    fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

    ds = xr.load_dataset(fp)

    ds.cfd.convert_species_rho()
