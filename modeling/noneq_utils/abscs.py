# -*- coding: utf-8 -*-
# TODO: fix variable/function cases, improve readability of var/func names
"""
Functions to calculate UV absorption cross sections     
"""
# stdlib imports
import os

# third party imports
import numpy as np
import pandas as pd
import xyzpy

def calc_ds_cs(
        Ts,
        wls=np.linspace(190, 285, 100)
):
    """
    TODO: improve documentation

    Parameters
    ----------
    Ts
    wls : np.ndarray

    Returns
    -------
    ds_cs
    """
    packagepath = os.path.dirname(os.path.realpath(__file__))
    abscs_folder = os.path.join(packagepath, 'Absorption Cross Sections')

    fp_KOH = os.path.join(
        abscs_folder,
        'KOH_Weng2019_1800K.csv'
    )
    fp_K = os.path.join(
        abscs_folder,
        'K cross section.csv'
    )
    fp_CO = os.path.join(
        abscs_folder,
        'CO_Thompson(1963)_298K_175-207nm.txt'
    )

    cs_KOH1800 = pd.read_csv(
        fp_KOH,
        header=None,
        index_col=0
    )[1]
    cs_K = pd.read_csv(
        fp_K,
        header=None,
        index_col=0
    )[1]
    cs_CO = pd.read_csv(
        fp_CO,
        header=None,
        delimiter="\t",
        index_col=0
    )[1].dropna()

    dummy = pd.Series(index=wls)

    cs_K = pd.concat([cs_K, dummy]).sort_index().interpolate(
        'linear',
        limit_direction='both'
    ).reindex(wls)
    cs_K.name = 'K CS'
    cs_K.index.name = 'wavelength'
    cs_KOH1800 = pd.concat([
        cs_KOH1800,
        dummy
    ]).sort_index().interpolate(
        'linear',
        limit_direction='both'
    ).reindex(wls)
    cs_KOH1800.name = 'KOH CS'
    cs_KOH1800.index.name = 'wavelength'
    s = pd.concat([cs_CO, dummy]).sort_index().interpolate()
    cs_CO = s[~s.index.duplicated()].reindex(wls)
    cs_CO.name = 'CO CS'
    cs_CO.index.name = 'wavelength'

    combos = {
        'T': Ts,
        'wl': wls
    }

    r = xyzpy.Runner(
        calc_css,
        var_names=['K', 'CO2', 'H2O', 'KOH', 'O2', 'CO'],
        constants={
            'cs_K': cs_K,
            'cs_KOH1800': cs_KOH1800,
            'cs_CO': cs_CO
        }
    )

    ds_cs = r.run_combos(combos)

    ds_cs.coords['T'].attrs = dict(
        long_name='Temperature',
        units='K'
    )
    ds_cs.coords['wl'].attrs = dict(
        long_name='Wavelength',
        units='nm'
    )

    for var in ds_cs.data_vars:
        ds_cs[var].name = 'cross section'
        ds_cs[var].attrs = dict(
            long_name='Absorption Cross Section',
            units='cm^2/molecule'
        )

    return ds_cs


def calc_cs_co2(
        wl,
        T
):
    """
    TODO: improve documentation
    Parameters
    ----------
    wl
    T

    Returns
    -------

    """
    wl = wl / 100
    T = T / 1000
    a = 0.05449 + 0.13766 * T + (23.529 / T)
    b = 1.991 - 0.17125 * T - (14.694 / T)
    return 1e-19 * np.exp((a + b * wl))


def calc_cs_h2o(
        wl,
        T
):
    """
    TODO: improve documentation

    Parameters
    ----------
    wl
    T

    Returns
    -------

    """
    sig0 = 10 ** (-15.598 - 0.0103 * wl)
    exponent = -(33664 - 202.6 * wl) / T
    return sig0 * np.exp(-exponent)


def calc_cs_o2(
        wl,
        T
):
    """
    TODO: improve documentation

    Parameters
    ----------
    wl
    T

    Returns
    -------

    """
    sig0 = 10 ** (-17.05)
    prefactor = 1 - np.exp(-2240 / T)
    exponent = (wl * 7.74e-2 - 12.31) * (2240 / T)
    return sig0 * prefactor * np.exp(-exponent)


def calc_css(
        wl,
        T,
        cs_K,
        cs_KOH1800,
        cs_CO
):
    """
    TODO: improve documentation

    Parameters
    ----------
    wl
    T
    cs_K
    cs_KOH1800
    cs_CO

    Returns
    -------

    """
    H2O = calc_cs_h2o(wl, T)
    CO2 = calc_cs_co2(wl, T)
    O2 = calc_cs_o2(wl, T)
    K = cs_K[wl]
    KOH = cs_KOH1800[wl]
    CO = cs_CO[wl]

    return K, CO2, H2O, KOH, O2, CO


if __name__ == '__main__':
    ds = calc_ds_cs(np.linspace(1000,3000))
    print(ds)