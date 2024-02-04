"""
Summary
-------
Modification of cross-section code in absem.fitting
- split function into 2 seperate function
    1. cross-section
    2. absorption
- also calculate intermediate absorption coefficient kappa

TODO
- make modifications to mhdpy, remove code and update link to mhdpy
    - code copied over to avoid issues with git-modules

"""
from pint import Quantity
import numpy as np
from lmfit.lineshapes import voigt
from scipy.ndimage import gaussian_filter1d


def Q_2peak(wave, **kwargs):
    """Model of the absorption of potassium peaks with two Voigt profiles.

    Parameter
    ---------
    wave : float array or pint float array
        wavelength [nm]

    Returns
    -------
    Q : float array
        absorption cross section [1/nm^2]
    """

    if isinstance(wave, Quantity):
        x = wave.to('nm').magnitude
    else:
        x = wave

    p1_stat_weight = Quantity((4/6)*0.7 , 'dimensionless')
    p1_center = Quantity(766.504333 , 'nm')
    p1_delta_f_C = Quantity( 3.19e9 , 'Hz')
    p2_stat_weight = Quantity(0.35 , 'dimensionless')
    p2_center = Quantity( 769.913219 , 'nm')
    p2_delta_f_C = Quantity(3.04e9 , 'Hz')
    const_absorb = Quantity(0, '1/nm')

    # Note there is a c**2 in this prefactor due to conversion of Sv (voigt) from frequency to wavelength space. See Alemade 1982 II.209
    e = Quantity(1.6e-19, 'C')
    me = Quantity(9.1e-31, 'kg')
    epsilon0 = Quantity(8.85e-12, 'F/m')
    c = Quantity(3e8, 'm/s')
    f0 = (e**2)/(4*epsilon0*me*(c**2))

    T = Quantity(3000, 'K')
    kb = Quantity(1.38e-23, 'J/K')
    K_MW        = Quantity(39.1, 'g/mol')
    NA = Quantity(6.022E23, '1/mol')
    mK = K_MW/NA

    p1_delta_wl_D = 2*np.sqrt(2*np.log(2)*kb*T/mK)*(p1_center/c)
    p1_delta_wl_C = ((p1_center**2)/c)*p1_delta_f_C
    p1_prefactor = p1_stat_weight*f0*p1_center**2

    p2_delta_wl_D = 2*np.sqrt(2*np.log(2)*kb*T/mK)*(p2_center/c)
    p2_delta_wl_C = ((p2_center**2)/c)*p2_delta_f_C
    p2_prefactor = p2_stat_weight*f0*p2_center**2

    # voigt profile over wavelength integrates to 1. int[S_{lambda}d_{lambda}] = 1 Can check with np.trapz(voigt, x)
    # Spectral distribution function, Called S_lambda in Alkemade 1982
    S_p1 = voigt(x, 1, p1_center.to('nm').magnitude, p1_delta_wl_D.to('nm').magnitude, p1_delta_wl_C.to('nm').magnitude)
    S_p1 = Quantity(S_p1, '1/nm')

    S_p2 = voigt(x, 1, p2_center.to('nm').magnitude, p2_delta_wl_D.to('nm').magnitude, p2_delta_wl_C.to('nm').magnitude)
    S_p2 = Quantity(S_p2, '1/nm')

    Q = p1_prefactor*S_p1 + p2_prefactor*S_p2

    Q.ito("nm^2")
    return Q

    #Atomic absorption cross section (~nm^2)
    #kappa = p1_prefactor*S_p1 + p2_prefactor*S_p2

def alpha_2peak(x,
    nK =1e-7,
    L=1e7,
    assert_xs_equal_spacing=True,
    **kwargs):
    """
    Parameter
    ---------
    x : float array
        wavelength [nm]
    nK : float
        number density [1/nm^3]
    L : float
        length [nm]

    Returns
    -------
    alpha : float array
        absorption []
    tau : float array
        transmission []
    kappa : float array
        absorption coefficient [1/nm]
    Q : float array
        cross section
    """
    Q = Q_2peak(x)

    kappa0 = Quantity(0, '1/nm')


    nK = Quantity(nK, '1/nm^3')
    L = Quantity(L, 'nm')

    kappa = Q*nK + kappa0

    tau = Quantity(np.exp(-kappa*L), 'dimensionless')

    alpha = -tau + Quantity(1, 'dimensionless')

    # Final gaussian blurring of the spectrum.
    # The final spectrum theoretically determined by alpha_2peak needs to be convoluted with the spectrometer response i.e. blurring.
    # See interp_alpha for more information.
    blur_pixels_amount = 1 # Was a parameter...
    if blur_pixels_amount:

        print("before assert")

        if assert_xs_equal_spacing:
            x_deltas = np.diff(x)
            assert np.close(x_deltas, x_deltas[0])
            assert np.isclose(x_deltas[0], 0.026*blur_pixels_amount)
        alpha[:] = gaussian_filter1d(alpha.magnitude, blur_pixels_amount)

    return {"alpha":alpha, "tau":tau, "kappa":kappa, "Q":Q}


