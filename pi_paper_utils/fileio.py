"""
Standard pipelines for importing processed datasets
"""

from mhdpy.analysis.standard_import import *

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
from mhdpy.pyvista_utils import CFDDatasetAccessor


def load_absem(
        tc, 
        avg_mnum = False,
        wl_slice = slice(750,790)
        ):
    """
    pipeline for loading absem data

    Parameters
    ----------
    tc : str
        test case. determines which file to load
    avg_mnum : bool, optional
        if True, average over mnum before calculation of alpha
    """

    ds_absem = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'absem','{}.cdf'.format(tc)))

    # #TODO: having to do this on office comp?
    # ds_absem.coords['mp'] = ds_absem.coords['mp'].astype(str)
    # ds_absem.coords['date'] = ds_absem.coords['date'].astype(str)

    ds_absem = ds_absem.xr_utils.stack_run()

    ds_absem = ds_absem.sel(wavelength=wl_slice)

    if tc == '53x':
        ds_absem = ds_absem.drop(0, 'kwt')

    if avg_mnum:
        ds_absem = ds_absem.mean('mnum')
    
    ds_absem = ds_absem.absem.calc_alpha()

    return ds_absem

def load_lecroy(tc, avg_mnum = False):
    """
    pipeline for loading lecroy data

    Parameters
    ----------
    tc : str
        test case. determines which file to load
    avg_mnum : bool, optional
        if True, average over mnum before calculation of mag_phase_AS
    """
    ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))


    ds_lecroy = ds_lecroy.xr_utils.stack_run()

    if tc == '53x':
        ds_lecroy = ds_lecroy.drop(0,'kwt')

    if avg_mnum:
        ds_lecroy = ds_lecroy.mean('mnum')

    ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected
    ds_lecroy = ds_lecroy.mws.calc_mag_phase_AS() # calculate before or after mnum mean?

    return ds_lecroy

def load_cfd_centerline(kwt_interp = None):
    """
    pipeline for loading cfd torch axis data

    Parameters
    ----------
    kwt_interp : xr.DataArray, optional
        if not None, interpolate the cfd data to the kwt values in this array
    """

    fp = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_torchaxis_Yeq.cdf')

    ds_cfd = xr.load_dataset(fp)

    if kwt_interp is not None:
        ds_cfd = ds_cfd.interp(kwt=kwt_interp).dropna('kwt', how='all')
        # ds_cfd = ds_cfd.interp(kwt=ds_lecroy.coords['kwt']).dropna('kwt', how='all')

    ds_cfd = ds_cfd.cfd.quantify_default()
    ds_cfd = ds_cfd.cfd.convert_all_rho_number()


    all_K_species = ['Yeq_K','Yeq_K+','Yeq_K2CO3','Yeq_KO','Yeq_KOH']
    ds_cfd['all_K_Yeq'] = sum([ds_cfd[field] for field in all_K_species])

    ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()
    return ds_cfd


