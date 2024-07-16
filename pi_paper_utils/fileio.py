"""
Standard pipelines for importing processed datasets
"""

from mhdpy.analysis.standard_import import *

DIR_EXPT_PROC_DATA = pjoin(REPO_DIR, 'experiment', 'data','proc_data')
from mhdpy.pyvista_utils import CFDDatasetAccessor

from .constants import MAG_MWS_BLOCK


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

from mhdpy.coords.ct import downselect_acq_time

def load_lecroy(tc, 
                avg_mnum = False, 
                AS_calc = None,
                df_ct_downselect = None,
                ):
    """
    pipeline for loading lecroy data

    Parameters
    ----------
    tc : str
        test case. determines which file to load
    avg_mnum : bool, optional
        if True, average over mnum before calculation of mag_phase_AS
    AS_calc : str, optional
        if not None, calculate the AS using this method. Options are 'relative' or 'absolute'
    df_ct_downselect : pd.DataFrame, optional
        if not None, downselect the acquisition times to the times in this dataframe, before averaging over mnum which removes 'acq_time' dimension
    """
    ds_lecroy = xr.load_dataset(pjoin(DIR_EXPT_PROC_DATA, 'lecroy','{}.cdf'.format(tc)))

    if tc == '53x':
        ds_lecroy = ds_lecroy.drop(0,'kwt')

    if df_ct_downselect is not None:
        ds_lecroy = downselect_acq_time(ds_lecroy, df_ct_downselect)

    # Calculate mnum dependent stats, which requires calculating mag.  
    ds_lecroy = ds_lecroy.mws.calc_mag_phase()

    if avg_mnum:
        ds_lecroy = ds_lecroy.mean('mnum', keep_attrs=True)

    ds_lecroy = ds_lecroy.sortby('time') # Needed otherwise pre pulse time cannot be selected

    # calculate before or after mnum mean?
    if AS_calc is not None:
        assert AS_calc in ['relative', 'absolute'], 'AS_calc must be "relative" or "absolute"'
        if AS_calc == 'absolute':
            da_nothing = load_mws_T0()
            ds_lecroy = ds_lecroy.mws.calc_AS_abs(mag_0=da_nothing, mag_block=MAG_MWS_BLOCK)
        else:
            ds_lecroy = ds_lecroy.mws.calc_AS_rel() 

    if 'pd1' in ds_lecroy:
        ds_lecroy = ds_lecroy.mws.calc_dPD()

    ds_lecroy = ds_lecroy.xr_utils.stack_run()

    return ds_lecroy

def load_mws_T0():
    """
    pipeline for loading mws T0 data
    """
    fp_nothing = pjoin(REPO_DIR,'final', 'dataset', 'output', 'mws_T0.csv')

    df_nothing = pd.read_csv(fp_nothing, index_col=0)['0']
    da_nothing = xr.DataArray(df_nothing).pint.quantify('volt')

    return da_nothing

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

    all_K_species = ['K','Kp','K2CO3', 'KO','KOH']
    ds_cfd['all_K'] = sum([ds_cfd[field] for field in all_K_species])

    ds_cfd['rho_number'] = ds_cfd.cfd.calc_rho_number()
    return ds_cfd


from mhdpy.pyvista_utils import CFDDatasetAccessor

def load_cfd_beam(kwt_interp = None, convert_rho_number = True):


    fp_cfd_profiles = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_beam_Yeq.cdf')
    ds_cfd_beam_mobile = xr.load_dataset(fp_cfd_profiles)

    fp_cfd_profiles = pjoin(REPO_DIR, 'final', 'dataset', 'output', 'line_profiles_beam_barrelexit_Yeq.cdf')
    ds_cfd_beam_barrel = xr.load_dataset(fp_cfd_profiles)


    ds_cfd_beam = xr.concat([
                        ds_cfd_beam_mobile.assign_coords(mp='mw_horns'),
                        ds_cfd_beam_barrel.assign_coords(mp='barrel')
                        ], dim='mp')

    if kwt_interp is not None:
        ds_cfd_beam = ds_cfd_beam.interp(kwt=kwt_interp).dropna('kwt', how='all')


    ds_cfd_beam = ds_cfd_beam.cfd.quantify_default()

    if convert_rho_number:
        ds_cfd_beam = ds_cfd_beam.cfd.convert_all_rho_number()

    ds_cfd_beam.coords['dist'] = ds_cfd_beam.coords['dist'].pint.to('cm') # Fitting expects cm

    return ds_cfd_beam


results_dir = pjoin(os.getenv('CFD_RAW_FOLDER'), 'medium_12July24')

basename = 'sweep' #sweepK for coarse run, sweep otherwise
filename = 'frontCyl_chem2.vtk' # chem2 for medium, chem1 for coarse

cfd_fp_dict = {
    '0.8_0.0': pjoin(results_dir, f'{basename}00_mdot0131_phi0801_K0000', filename),
    '0.8_0.06': pjoin(results_dir, f'{basename}01_mdot0132_phi0793_K0006', filename),
    '0.8_0.1': pjoin(results_dir, f'{basename}02_mdot0131_phi0799_K0010', filename),
    '0.8_0.18': pjoin(results_dir, f'{basename}03_mdot0132_phi0789_K0018', filename),
    '0.8_0.31': pjoin(results_dir, f'{basename}04_mdot0131_phi0797_K0031', filename),
    '0.8_0.55': pjoin(results_dir, f'{basename}05_mdot0132_phi0782_K0055', filename),
    '0.8_1.0': pjoin(results_dir, f'{basename}07_mdot0131_phi0787_K0099', filename),
    '0.6_1.0': pjoin(results_dir, f'{basename}06_mdot0132_phi0584_K0099', filename),
    '0.8_1.75': pjoin(results_dir, f'{basename}08_mdot0133_phi0759_K0175', filename),
}

# soi = ['K', 'K_p', 'e_m', 'OH', 'OH_m', 'KOH', 'O2', 'H2O', 'N2', 'CO2'] # works for coarse
# soi = ['K', 'K+', 'e-', 'OH', 'OH-', 'KOH', 'O2', 'H2O', 'N2', 'CO2'] # works for 24Apr24_corr
soi = ['K',  'Kp', 'K2CO3', 'KO', 'KOH', 'em', 'OH', 'OHm', 'O2', 'H2O', 'N2', 'CO2'] # works for july 12
# soi = ['K',  'OH', 'KOH', 'O2', 'H2O', 'N2', 'CO2'] # works for june08
soi_Yeq = ['Yeq_K', 'Yeq_K+', 'Yeq_e-', 'Yeq_OH', 'Yeq_OH-', 'Yeq_KOH', 'Yeq_K2CO3', 'Yeq_KO']
additional = ['T', 'p', 'rho']
cfd_all_fields = [*soi, *soi_Yeq, *additional]

# sp_dir = gen_path('sharepoint')
# results_dir = pjoin(sp_dir, 'Team Member Files', 'DaveH', 'Results', 'axiJP8200_17Jul23')

# fps = {
#     '0.8_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi080_K010', 'frontCyl_chem1.vtk'),
#     '0.8_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi080_K100', 'frontCyl_chem1.vtk'),
#     '0.6_0.1': pjoin(results_dir, 'medium', 'mdot0130_phi060_K010', 'frontCyl_chem1.vtk'),
#     '0.6_1.0': pjoin(results_dir, 'medium', 'mdot0130_phi060_K100', 'frontCyl_chem1.vtk'),
# }