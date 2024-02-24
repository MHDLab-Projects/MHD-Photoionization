import os
import pandas as pd
import numpy as np
import xarray as xr

e = 1.6e-19
E_IP = 4.34
K_MHD = 0.5

def calc_G_NE(P_in, eta):
    G_NE = eta*P_in/(e*E_IP)
    return G_NE

def calc_ne(G_th, G_NE, krb):
    ne = np.sqrt((G_NE + G_th)/krb)
    return ne

def calc_alpha(ne, mue, krb, u, B, eta):
    """
    calculates alpha=dP_net/dP_in
    """

    DPmhd_Dne = e*K_MHD*(1-K_MHD)*(u**2)*(B**2)*mue

    DPR_Dne = 2*e*E_IP*krb*ne

    alpha = eta*(DPmhd_Dne/DPR_Dne)

    return alpha

def calc_alpha_const_nx(krm, mue, u, B, eta):
    """
    calculates alpha=dP_net/dP_in for a constant nx, where x is recombination partner
    """

    DPmhd_Dne = e*K_MHD*(1-K_MHD)*(u**2)*(B**2)*mue

    DPR_Dne = e*E_IP*krm

    alpha = eta*(DPmhd_Dne/DPR_Dne)

    return alpha

def calc_ne_const_nx(ne0, krm, G_NE):
    """
    calculates ne for a constant nx, where x is recombination partner 
    """

    #TODO: replace ne0 with calc_ne G_NE=0 and compare

    ne = ne0 + G_NE/krm

    return ne


def calc_NE_all_const_nx(P_in, eta, krm, mue_cant, B, u, ne0):

    G_NE = calc_G_NE(P_in, eta)

    ne = calc_ne_const_nx(ne0, krm, G_NE)
    ne.name = 'ne'
    ne.attrs = dict(units = '$\#/cm^3$', long_name = 'Calculated electron density')

    sig_NE = calc_sig(ne, mue_cant)
    sig_NE.name = 'sigma'

    alpha = calc_alpha_const_nx(krm, mue_cant, u, B, eta)
    alpha.name = 'alpha'

    ds = xr.merge([ne, sig_NE, alpha])

    return ds


def calc_NE_all(P_in, eta, G_th, krb, mue_cant, B, u):
    """
    # If wl in combos calculate eta based on photoionization, vs if not use eta = 1
    # Further, If L and R in combos calculate fractional absorption (FA) based on optical cavity approach, else FA = 1
    # If P_in not in combos P_in = 0
    # If l_bk (fractional length of bulk vs contact) not in combos l_bk = 0 (Tbulk = 3000K)
    """
    
    G_NE = calc_G_NE(P_in, eta)

    # G_NE.name = 'G_NE'
    # G_NE.attrs = dict(units = '$\#/cm^3 s$', long_name = 'Maximum Generation Rate')
    # G_NE.coords['P_in'].attrs = dict(units = '$W/cm^3$', long_name = 'Input Power Density')

    ne = calc_ne(G_th, G_NE, krb)
    ne.name = 'ne'
    ne.attrs = dict(units = '$\#/cm^3$', long_name = 'Calculated electron density')

    sig_NE = calc_sig(ne, mue_cant)
    sig_NE.name = 'sigma'

    alpha = calc_alpha(ne, mue_cant, krb, u, B, eta)
    alpha.name = 'alpha'

    #TODO: Can't figure out to get xyzpy to work when passing these as separate returns, had to use var_names= None 
    ds = xr.merge([ne, sig_NE, alpha])

    return ds

## Legacy

def calc_P_mhd(sig, u, B):
    """
    calculates mhd power density from simple 0d formula
    u is fluid velocity in cm/s
    """
    fac = K_MHD*(1-K_MHD)*(u**2)*(B**2)
    P_mhd = fac*sig
    P_mhd = P_mhd.assign_attrs(units = '$W/cm^3$', long_name = 'MHD power')
    return P_mhd

def calc_ureq(P_mhd_set, sig, B):
    """
    calculates velocity required for specific power output
    P_mhd_set: set power output in W/cm^3
    """
    fac = K_MHD*(1-K_MHD)*(B**2)*sig
    ureq = np.sqrt(P_mhd_set/fac)
    ureq = ureq.assign_attrs(long_name = 'required velocity for ' + str(P_mhd_set) + " W/cm^3", units='cm/s')
    return ureq


def calc_P_net(P_mhd):
    P_in = P_mhd.coords['P_in']
    P_net = P_mhd - P_in
    P_net.name = 'P_net'
    P_net = P_net.assign_attrs(units = '$W/cm^3$', long_name = 'Net MHD Power (Total - Laser in)')
    return P_net


def calc_eff_NE(P_mhd):
    delta_P_mhd = P_mhd- P_mhd.sel(P_in = 0)
    P_in = P_mhd.coords['P_in']
    # eff_NE = (delta_P_mhd - P_in)/P_in
    eff_NE =   delta_P_mhd/P_in
    eff_NE.name = 'eff_NE'
    eff_NE.attrs = dict(long_name = 'NonEq Generation Efficiency')
    return eff_NE

def calc_ne_finite(G_th, G_NE, krb, K):
    """Calculates ne but with a finite number of K that can be ionized"""
    k_NE = G_NE/K
    k_th = G_th/K

    frac = (k_NE + k_th)/(krb)
    ne_finite = (frac/2)*(np.sqrt(1+(4/frac)*K) - 1)
    ne_finite.name = 'ne_finite'
    ne_finite.attrs = dict(units = '$\#/cm^3$', long_name = 'Calculated finite electron density')

    return ne_finite

def calc_sig(ne,mue):
    #Need to compare this conductivity to cantera? 
    sig = e*mue*ne
    sig = sig.assign_attrs(units = '$S/cm$', long_name = 'Conductivity')
    return sig


def calc_sigma_tot(sigma_bl, sigma_bk, l_bk):
    num = sigma_bl*sigma_bk
    dem = sigma_bk + l_bk*(sigma_bl - sigma_bk)
    return num/dem


def calc_dsigma_tot(sigma_bl, sigma_bk, l_bk):
    num = (1-l_bk)*sigma_bk**2
    dem = ( sigma_bk + l_bk*(sigma_bl - sigma_bk) )**2
    return num/dem

### Photoionization ###

def calc_atten_lengths(ds_cs, ds_species_rho):
    """Caclulate optical attenuation length for gas mixture defined by ds_species_rho (#/cm^3)"""

    abstot = ds_species_rho['K'].where(False).fillna(0)
    for species in ds_cs.data_vars:
    #     if species == 'K':
    #         abstot = abstot + ds_cs[species]*ds_species_rho[species]
        abstot = abstot + ds_cs[species]*ds_species_rho[species]
        
    lamtot = 1/abstot
        
    lamtot.name = 'tot'
    lamtot.attrs = dict(long_name = 'Attenuation Length', units = 'cm')

    gas_lam = lamtot.to_dataset()

    for species in ds_cs.data_vars:
        lam = 1/(ds_species_rho[species]*ds_cs[species])
        lam.attrs = dict(long_name = 'Attenuation Length from ' + species, units = 'cm')
        gas_lam = gas_lam.assign(temp=lam).rename(temp=species)

    return gas_lam


def calc_FA(lamtot, L, R):
    FAL = (1-np.exp(-L/lamtot))
    FA = FAL/(1-(1-FAL)*R)
    FA.attrs = dict(long_name = 'Total Fraction Absorbed')
    FA.name = 'FA'

    return FA

def calc_eta_PI(lamtot, lamK, FA):
    E_ph = 1240/lamtot.coords['wl']*e
    # Ratio of attenuation length of total mixture to potassum is the  the fractional absorbance of potassium. Have a onenote page that proves this is valid
    eta_PI = (lamtot/lamK)*FA*(e*E_IP/E_ph)
    return eta_PI


