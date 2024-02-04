# -*- coding: utf-8 -*-
"""
Summary
-------
Absorption Cross Sections


"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants

unit_nm = 1e-9

def mesh_grid_ndim(a, b):
    new_shape = a.shape + b.shape
    a1, b1 = np.meshgrid(a, b)
    a1 = a1.T.reshape(new_shape)
    b1 = b1.T.reshape(new_shape)
    return a1, b1

def Q_absorption_constant(T, lam_nm):
    """absorption cross-section


    """
    qe = ct.electron_charge
    me = ct.electron_mass
    c0 = ct.light_speed
    eps0 = ct.permeability_0
    C = qe**2/(4.0*me*c0)
    return (1.0 + T*0.0 + lam_nm*0.0)

class GaussAbsorption:
    """

    Calculates the relative cross-section assuming a series of
    Gaussian peaks defined by
        A_i   - amplitude
        w_i   - width
        lam_i - wavelength center
        b_i   - temperature exponent

        Q0 = 1e-18 m^2 [default]
        T0 = 1000 K

        $$
        beta_i = beta0_i (T/T0)^b
        Q/Q0  = sum_i A_i exp ( -(beta_i (lam - lam_i)^ 2)
        $$


    """

    def __init__(self, A=[2.0/3.0,1.0], lam=[765.0,770.0], w=[8.0,8.0], b=[0.0,0.0]):
        """

        Parameters
        ----------
        A : float array, optional
            peak amplitude []. The default is [2.0/3.0,1.0].
        lam : float array, optional
            peak wavelength[nm]. The default is [765.0,770.0].
        w : float array, optional
            peak width @ 1000K. The default is [8.0,8.0].
        b : TYPE, optional
            temperature exponent. The default is 0.0.

        Returns
        -------
        None.

        """
        self.A = np.array(A)
        self.lam = np.array(lam)
        self.w = np.array(w)
        self.beta = np.sqrt(-np.log(0.1))/self.w*4.0
        self.b = np.array(b)
        self.T0 = 1000.0
        self.Q0 = 1e-18
        self.A_min = 1e-9



    def calc(self, T, lam_nm):
        """calculate cross section"""

        T0, b, beta, A, lam = self.T0, self.b, self.beta, self.A, self.lam

        t = T/T0
        t1, w1 = mesh_grid_ndim(t, lam_nm)
        beta_list = [beta_i*t1**b[i] for i, beta_i in enumerate(beta)]
        #beta_list = beta*t**b
        print(t.shape)
        print(beta_list[0].shape)
        #

        Q_list = [ A[i]*np.exp( -(beta_i*(w1-lam[i]))**2) for i,beta_i in enumerate(beta_list)]
        #Q_list =
        Q = np.array(Q_list)
        return Q.sum(axis=0)*self.Q0


    def __call__(self, T, lam_nm):
        """return nornmalized cross section"""
        beta = self.beta[:,np.newaxis]*(T/self.T0)**self.b[:,np.newaxis]
        Q = self.A[:,np.newaxis]*np.exp(-(beta*(lam_nm-self.lam[:,np.newaxis]))**2)
        Q += self.A_min
        return np.sum(Q, axis=0)

    def cross_section(self, T, lam_nm):
        """return dimensional cross section"""
        return self(T, lam_nm)*self.Q0

    def plot(self, T=[300.0,1000.0,2000.0], ax=None):
        """

        Parameters
        ----------
        T : float array
            temperature [K]
        ax : matplotlib figure axes, optional none
            if ax == None, figure and axes are generated
        Returns
        -------
        ax

        """
        lam, w = self.lam, self.w
        p = []
        n_peaks = len(lam)
        for i in range(n_peaks):
            p.append( np.linspace(lam[i]-w[i],lam[i]+w[i]))
        p = np.sort(np.unique(np.array(p)))
        if ax is None:
            fig, ax = plt.subplots()
        for t in T:
            Q = self(t,p)
            ax.plot(p, Q, label="{}".format(t))
            ax.set_xlabel("wavelength [nm]")
            ax.set_ylabel("cross section []")
        ax.legend(title="T [K]")
        return ax


class VoigtAbsorption:

    p_atm = 101325.0

    def __init__(self):
        self.center = np.array([766.504333, 769.913219])
        self.stat_weight = np.array([4.0/6.0*0.7, 0.35])
        self.Deltaf = np.array([3.19e9,3.04e9])
        self.MW = 39.1
        self.update()

    def update(self, p=1e5, T_ref=1000.0):
        unit_nm = 1e-9
        k_B = constants.Boltzmann/unit_nm**2
        c = constants.c/unit_nm
        N_A = constants.Avogadro

        self.mass = self.MW*N_A
        self.T_ref = T_ref
        self.sigma = 2.0*np.sqrt(2.0*np.log(2.0)*k_B*T_ref/self.mass)*(self.center/c)
        self.gamma = (self.center**2/c)*self.Deltaf

    def calc(self, T=300.0, p=1e5, lam_nm=760.0, verbose=True, do_plot=False):
        """calculate absorption cross-section

        Parameters
        ----------
        T : float
            temperture [K]
        p : float
            pressure [Pa]
        lam_nm : float
            wavelength [nm] - not used
        verbose : int or bool
            output to screen
        do_plot : bool
            plot coefficients

        Ref
        ---
        Hanson, Lecture Notes
        https://cefrc.princeton.edu/sites/g/files/toruqf1071/files/Files/2013%20Lecture%20Notes/Hanson/pLecture6.pdf

        """
        k_B = constants.Boltzmann
        c = constants.c

        lam_0 = self.center*unit_nm
        nu_0 = c/lam_0

        tau_u = 1e-8  # electronic
        tau_u = 1e-2  # vib-rot
        delta_nu_N = 1.0/(2.0*np.pi)/tau_u
        # collisional broadening, Lorentzian
        # Z_B = p*np.pi*Q_AB*np.sqrt(8.0/np.pi*mu_AB*k_B*T)
        # delta_nu_C = Z_B/np.pi
        unit_cm = 1e-2
        #     [1/(atm-m)] = [1/(atm-cm)] / [cm/m]
        two_gamma = 0.1/unit_cm
        #     [1/m] = [1/(atm-m)] [ ] [Pa/(Pa/atm)]
        delta_k_C = two_gamma*(300.0/T)**0.5*(p/self.p_atm)
        #     [1/s] = [1/m]*[m/s]
        delta_nu_C = delta_k_C*c

        delta_lam_C = (lam_0/nu_0)*delta_nu_C
        if verbose:
            print("Broadening")
            print("collisional")
            print("freq [1/s]", delta_nu_C)
            print("wavelength [nm]", delta_lam_C/unit_nm)
            print("wave# [1/cm]", delta_k_C*unit_cm)

        # dopler broadening
        mass = constants.m_p*self.MW
        nu_0 = c/(self.center*unit_nm)
        delta_nu_D = 2.0*np.sqrt(2.0*k_B*T*np.log(2)/(mass*c**2))*nu_0


        lam_0 = self.center*unit_nm
        delta_lam_D = delta_nu_D*lam_0/nu_0
        delta_k_D = delta_nu_D/c
        if verbose:
            print("doppler")
            print("freq [1/s]", delta_nu_D)
            print("wavelength [nm]", delta_lam_D/unit_nm)
            print("wave# [1/cm]", delta_k_D*unit_cm)

        if do_plot:
            fig, ax = plt.subplots(2,2)

        for i in range(2):
            wave = np.linspace(self.center[0] - 5*delta_lam_D[0]/unit_nm, self.center[1] + delta_lam_D[1]/unit_nm*5.0, 100)
            wave = np.linspace(761,775,1000)
            nu = c/(wave*unit_nm)
            phi = [ delta_nu_C/( (nu - nu_0[i])**2 + (delta_nu_C/2)**2) for i in range(2)]

            phi_D = [ np.exp( -( 2.0*np.sqrt(np.log(2))*(nu - nu_0[i])/delta_nu_D[i] )**2 ) for i in range(2)]
            #print(np.max(phi_D))

            if do_plot:
                ax[0,0].plot(wave, phi[i]/phi[i].max(), label="coll.".format(i))
                ax[1,0].plot(wave, phi_D[i]/phi_D[i].max(), label="dopler{}".format(i))

                if i == 0:
                    for j, f in enumerate([phi[i], phi_D[i]]):
                        kappa = f/f.max()
                        for A in [1e2,1e3,1e4]:
                            ax[j,1].plot(wave, np.exp(-kappa*A), label="{}".format(A))
                        ax[j,1].legend()

        if do_plot:
            ax[0,1].legend(title="A")
            ax[1,1].legend(title="A")

            for i, name in enumerate(["Collisional","Doppler"]):
                ax[i,0].set_ylabel(r"$\kappa$")
                ax[i,1].set_ylabel(r"$\exp(-\kappa A)$")
                ax[i,0].set_title(name)
                ax[i,1].set_title(name)

            fig.tight_layout()

    def __call__(self, T, lam):
        """voigt(x, amplitude=1.0, center=0.0, sigma=1.0, gamma=None)
        kappa = p1_prefactor*voigt(x, 1, p1_center, p1_delta_wl_D, p1_delta_wl_C) + p2_prefactor*voigt(x, 1, p2_center, p2_delta_wl_D, p2_delta_wl_C)
        """
        val = calc(lam, A[0], center[0], sigma[0], gamma[0])
        return val


def blurred_alpha(
    x, nK =1e-7, L=1e7, p1_stat_weight = (4/6)*0.7,
    p1_center = 766.504333, p1_delta_f_C = 3.19e9,p2_stat_weight = 0.35,
    p2_center = 769.913219, p2_delta_f_C = 3.04e9, const_absorb = 0
    ):
    """

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    nK : TYPE, optional
        DESCRIPTION. The default is 1e-7.
    L : TYPE, optional
        DESCRIPTION. The default is 1e7.
    p1_stat_weight : TYPE, optional
        DESCRIPTION. The default is (4/6)*0.7.
    p1_center : TYPE, optional
        DESCRIPTION. The default is 766.504333.
    p1_delta_f_C : TYPE, optional
        DESCRIPTION. The default is 3.19e9.
    p2_stat_weight : TYPE, optional
        DESCRIPTION. The default is 0.35.
    p2_center : TYPE, optional
        DESCRIPTION. The default is 769.913219.
    p2_delta_f_C : TYPE, optional
        DESCRIPTION. The default is 3.04e9.
    const_absorb : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    """
    # units converted from m to nm
    e = 1.6e-19
    me = 9.1e-31
    epsilon0 = 8.85e-39
    #c = 3e17
    c = constants.c/unit_nm
    f0 = (e**2)/(4*epsilon0*me*(c**2))

    T = 1000.0
    kb_old = 1.38e-5
    kb = constants.Boltzmann/unit_nm**2
    print("boltzmann",kb/kb_old)
    K_MW        = 39.1
    NA = 6.022E23
    mK_old         = (K_MW/1000)/NA
    mK = K_MW*constants.m_p
    print("mass K",mK/mK_old)
    # mK = 39*1.67e-27

    #TODO: convert to array math

    p1_delta_wl_D = 2.0*np.sqrt(2.0*np.log(2.0)*kb*T/mK)*(p1_center/c)
    p1_delta_wl_C = ((p1_center**2)/c)*p1_delta_f_C
    p1_prefactor = p1_stat_weight*f0*p1_center**2

    p2_delta_wl_D = 2*np.sqrt(2*np.log(2)*kb*T/mK)*(p2_center/c)
    p2_delta_wl_C = ((p2_center**2)/c)*p2_delta_f_C
    p2_prefactor = p2_stat_weight*f0*p2_center**2

    print("lam/c",p1_center/c)
    print("sigma {:8.2e}".format(p1_delta_wl_D))
    print("gamma {:8.2e}".format(p1_delta_wl_C))

    # rint("lam/c",p1_center/c)
    print("sigma {:8.2e}".format(p2_delta_wl_D))
    print("gamma {:8.2e}".format(p2_delta_wl_C))

    #kappa = p1_prefactor*voigt(x, 1, p1_center, p1_delta_wl_D, p1_delta_wl_C) + p2_prefactor*voigt(x, 1, p2_center, p2_delta_wl_D, p2_delta_wl_C)