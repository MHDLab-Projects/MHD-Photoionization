#%%
from mhdpy.analysis.standard_import import *

from mhdpy.pyvista_utils import CFDDatasetAccessor


def test_convert_rho():
    #TODO: create a file like this from a vtk, and then test that it loads correctly. Move tests to mhdpy
    fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

    ds = xr.load_dataset(fp)

    ds = ds.cfd.quantify_default()

    ds.cfd.convert_all_rho_number()

def test_calc_rho_number():

    #TODO: create a file like this from a vtk, and then test that it loads correctly. Move tests to mhdpy
    fp = pjoin(REPO_DIR, 'modeling', 'cfd','output', 'line_profiles_torchaxis_Yeq.cdf')

    ds = xr.load_dataset(fp)

    ds = ds.cfd.quantify_default()

    ds.cfd.calc_rho_number()
