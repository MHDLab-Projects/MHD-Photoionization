"""
Summary
-------
Use scipy tools for beam integration with absorption.

Testing of different methods to perform the ray integration.

Routines
--------


Notes
-----
Simpsons rule is much faster than the quadrature but it does not have
the strict error control

TODO
----
- add timer
- output beam meshes
- 2D plots for I(x_CL,lambda)
- should ionized potassium be included
- finite diameter beam
- vertical perturbation of source
- optimization testing
    - assume gaussian T and n profiles
    $$
        (T,n)(x,r) = (T,n)_CL(x) exp( -\beta_{T,n}(x) r^2 )
        (T,n)_CL(x) = x^a_{T,n}
        \beta_{T,n}(x) = x^b_{T,n}
    $$


https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.PolyDataFilters.ray_trace.html#pyvista.PolyDataFilters.ray_trace

"""
import os
from os.path import join as pjoin
import pyvista as pv
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

from mhdpy.pyvista_utils.axi import pv_interpolator, AxiInterpolator, line_interpolator
import cross_section

try:
    from mhdpy.fileio import gen_path
except Exception as e:
    print(e)

case = "mdot0130_phi080_K010"
try:
    sp_dir = gen_path('sharepoint')
    results_rel = r"Team Member Files\DaveH\Results\axiJP8200_17Jul23"
    results_dir = pjoin(sp_dir, results_rel)
    fname = os.path.join(results_dir, "medium",case,"frontCyl.vtk")

except Exception as e:
    print(e)
    results_dir = r"C:\Users\Huckaby\Desktop\MHD\Simulations\viz_example"
    sim_run = "."
    study = "."
    fname = os.path.join(results_dir,study,sim_run,"frontCyl.vtk")

n_dim = 3

class Source:

    def __init__(self, start=[1e-2,0,-1e-2], finish=[1e-2,0,1e-2], diameter=1e-3):

        self.start = np.array(start)
        self.finish = np.array(finish)
        self.diameter = diameter


if __name__ == "__main__":
    mesh = pv.read(fname)

    p = mesh.point_data
    R_u = ct.gas_constant
    N_A = ct.avogadro
    p['C_mix'] = p['p']/(p['T']*R_u)
    # TODO calculate mole fraction
    p['C_K'] = p['C_mix']*p['K']
    f_mesh = pv_interpolator(mesh, ["C_K","T"])

    f_axi = AxiInterpolator(mesh, ["C_K", "T"])


if __name__ == "__main__":
    Q_absorption = cross_section.GaussAbsorption(b=[-0.3,-0.3])
    Q_absorption.plot()


def f_quad(s, f_val, f_Qabs, lam_nm, kappa_CL):
    """integrand for integrate.quad

    Parameters
    ----------
    s : float
        position coordinate, *s*
    f_val(s) : function
        calculates an array containing field variables as a function of
        the ray position coordinate, *s*
        f_val(s)[0] : molar concentration or number density
        f_val(s)[1] : temperature
    f_Qabs(T,lam_nm) : function
        calculates the absorption cross-section as a function of
        temperature, *T* and wavelength, *lam*
    lam_nm : float
        ray wavelength [nm]
    kapaa_CL : float
        centerline absorption coefficient


    Returns
    -------
    kappa_star : float
        normalized absorption coefficient []

    """
    #pos = x0 + s*(x1 - x0)
    #r = (pos[1]**2 + pos[2]**2)**0.5
    #x = pos[0]
    #out = f_mesh((x, r))
    out = f_val(s)
    C_K, T = out[0], out[1]
    return (f_Qabs(T, lam_nm)*C_K)/kappa_CL


def axi_array(pos):
    r = (pos[:,1]**2 + pos[:,2]**2)**0.5
    pos_axi = pos*1.0
    pos_axi[:,1] = r
    pos_axi[:,2] = np.atan2(pos[:,1],pos[:2])
    return pos_axi

def test_quad_single(wavelength_nm=760.00):
    """integration of a single ray using quad function

    Parameters
    ----------
    wavelength_nm : float
        ray wavelength [nm]

    """
    test_quad_range([wavelength_nm])

_wave_range = np.linspace(760,775,31)
def test_range(wave_range=_wave_range):
    """integration of ray of different wavelengths using quad function

    Parameters
    ----------
    wave_range : float array
        ray wavelengths [nm]

    Returns
    -------
    kappa_star : float array
        normalized effective absorption coeffient
    kappaL : float array
        effective absorption coeffient beam length
    intensity_ratio : float array
        target intensity / source intensity

    Notes
    -----

    $ Q_K $ - absorption cross-section non-dimensionalized with $Q_ref$
    $ Q_ref $ - nominal absorption cross-section [m^2]


    $$
        \kappa_{CK} = Q_K( T_CL, \lambda) C_{K,CL}
    $$

    $$
        \kappa^* = \int_{a,b} ds ( Q_K( T(s), \lambda) C_K(s) )/ \kappa_{CL}
    $$

    $$
        (\kappa L) = \kappa^* (\kappa_{CL} N_A Q_{ref}) L
    $$

    $$
        \dfrac{I_t}{I_s} = \log \left(-\kappa L \right)
    $$



    """
    x_CL = np.array([0.3,0.0,0.0])
    x_offset = np.array([-0.05,0.0,-0.1])
    x_source = x_CL + x_offset
    x_target = x_CL - x_offset

    L = np.sum((x_target-x_source)**2)**0.5
    s_CL = np.sum((x_CL-x_source)**2)**0.5/L


    # barrel exit diameter
    d = 1.2e-3
    s = np.concatenate( [np.linspace(0,1,100),s_CL + np.linspace(-0.5,0.5,100)*(d/L)] )
    s = np.unique(s)
    s = np.sort(s)

    def f_pos(s):
        return x_source + s[:np.newaxis]*(x_target - x_source)

    pos = f_pos(s)
    pos_axi = axi_array(pos)
    out = f_mesh((pos_axi[:,:2]))
    fig, ax = plt.subplots(2,1,sharex=True)
    ax[0].semilogy(pos[:,2],out[:,0]+1e-9,'-')
    ax[1].plot(pos[:,2],out[:,1],'-')

    out = f_axi(x_CL[:])
    T_CL = out[1]
    C_CL = out[0]

    def f_val(s):
        pos = f_pos(s)
        r = (pos[1]**2 + pos[2]**2)**0.5
        x = pos[0]
        return f_mesh((x, r))

    A = wave_range*0.0
    err = A*0.0
    Q_CL = Q_absorption(T_CL, wave_range)
    kappa_CL = Q_CL*C_CL
    #
    # If the temperature exponent is constant then all the integrals are the same
    #
    n_wave = len(wave_range)
    # additional points to help integration
    points = [s_CL-0.5*(d/L), s_CL, s_CL+0.5*(d/L)]
    for i, wave in enumerate(wave_range):
        out = integrate.quad(f_quad, 0, 1, args=(wave, f_val, Q_absorption, kappa_CL[i]),
                                 full_output=True, points=points, limit=200)
        A[i] = out[0]
        err[i] = out[1]

        val, err = out[0], out[1]
        quad_out = out[2]

        if len(out) > 3:
            # print error message
            print(out[3])

        if len(out) > 2:
            last = quad_out["last"]
            n_p = last + 1
            p_left = quad_out["alist"][:last]
            p_right = quad_out["blist"][:last]
            val_cell = quad_out["rlist"][:last]

            n_dim = 3
            s = np.unique( np.concatenate([p_left,p_right]) )
            pos = f_pos(s)
            pos_axi = axi_array(pos)
            out = f_mesh(pos_axi[:,:2])
            color = cmap(float(i)/float(n_wave))
            ax[0].semilogy(pos[:,2],out[:,0]+1e-9, 'o', label='quad end-points', color=color)
            ax[1].plot(pos[:,2],out[:,1],'o')

    ax[1].set_xlabel("x [m]")
    ax[0].set_ylabel("C_K [kmol/m^3]")
    ax[1].set_ylabel("T [K]")
    ax[0].legend()

    # L is needed because the integral is calculated using a non-dimensional distance, s
    # N_A is needed because the concentration is used instead of the number density
    A1 = (A*kappa_CL)*N_A*Q_absorption.Q0*L
    I_ratio = np.exp(-A1)
    if (n_wave > 1):
        fig, ax = plt.subplots(3,1,sharex=True)
        ax[0].plot(wave_range, A)
        ax[0].set_ylabel("raw integral")
        ax[1].plot(wave_range, A1)
        ax[1].set_ylabel(r"$\kappa L_B$")
        ax[2].plot(wave_range, I_ratio)
        ax[2].set_ylabel("intensity ratio")
        ax[2].set_xlabel("wavelength [nm]")
        tw = ax[0].twinx()
        tw.plot(wave_range, np.log10(err), 'g--')
        tw.set_ylabel("log10(error)")
    return {"kappaL_star":A, "kappaL":A1, "intensity_ratio":I }

_x_CL = np.zeros((9,n_dim))
_x_CL[:,0] = np.linspace(0.2,0.6,9)
_lam_nm = np.arange(760,775+0.1,0.5)
_method = {"integrate":"simpson", "interpolate":"beam"}

#line_interpolator = interpolate.LinearNDInterpolator
#import beam_utils
#line_interpolator = beam_utils.adapt_array 
#(s, v, pos, f_v, max_iter=10, rtol=1e-2, do_plot=False):
    
def integrate_beams(x_CL=_x_CL, lam_nm=_lam_nm, method=_method):
    """calculate relative beam intensity for wavelengths, lam[] and
    and centerline positions x_CL[]

    Parameters
    ----------
    x_CL : float array
        centerline positions [m]
    lam_nm : float arrray
        wave lengths [nm]

    Returns
    -------
    dictionary of 2D float arrays [i_pos,i_wave]
        I_ratio : intensity ratio
        kappaL  : effective absorption coefficeint x beam length
        kappa_star : float array

    TODO:
    1. split-out code from inner loop into 1 or more functions
    2. vectorize integration
    """

    n_pos, n_dim1 = x_CL.shape
    n_wave = len(lam_nm)

    # source_offset
    x_offset = np.array([-0.05,0.0,-0.1])
    x_source = x_CL + x_offset
    x_target = x_CL - x_offset
    L = np.sum((x_source - x_target)**2,axis=1)
    wave_range = lam_nm

    # numerical integrand
    A = np.zeros((n_pos, n_wave))
    err = np.zeros((n_pos, n_wave))
    out = f_axi(x_CL)
    C_CL = out[:,0]
    T_CL = out[:,1]
    A1 = A*0.0
    Q0 = Q_absorption.Q0

    m1 = method.get("integrate","simpsons")
    m2 = method.get("interpolate","beam")

    for j in range(n_pos):

        s_beam = np.linspace(0,1.0,21)

        def f_pos(s):
            return x_source[j] + (x_target - x_source)[j]*s[:,np.newaxis]

        # create a direction 1D interpolating function
        f_beam, fig = line_interpolator(f_pos, f_axi, max_iter=15, do_plot=1)

        fig.axes[0].set_title("x_CL = {}".format(x_CL[j]))
        s_beam = f_beam.x
        #ds = s_beam[1:] - s_beam[:-1]
        s_CL = np.sum((x_CL[j]-x_source[j])**2)**0.5/L[j]
        points = [s_CL]
        Q_CL = Q_absorption(T_CL[j], wave_range)

        print("Integrating @ x_CL = {} for {} wavelengths".format(x_CL[j], n_wave))

        def f_val(s):
            pos = f_pos(s)
            r = ((pos.T[1])**2 + (pos.T[2])**2).T
            return f_mesh((pos[:,0], r))

        kappa_CL = C_CL[j]*Q_CL
        for i, wave in enumerate(wave_range):

            if m1 == "quad":
                if m2 == "axi":
                    out = integrate.quad(f_quad, 0, 1, epsabs=1e-3, epsrel=1e-4,
                                     args=(f_val, Q_absorption, wave, kappa_CL[i]),
                                     full_output=True, points=points, limit=200)
                else:
                    out = integrate.quad(f_quad, 0, 1,
                                    args=(f_beam, Q_absorption, wave, kappa_CL[i]),
                                    full_output=True, points=points, limit=200)
                A[j,i] = out[0]
                err[j,i] = out[1]
            else:
                if m2 == "axi":
                    vals = f_val(s_beam)
                else:
                    vals = f_beam(s_beam)
                C_K, T = vals[0], vals[1]
                Q = Q_absorption(T, wave)/Q_CL[i]
                a = Q*(C_K/C_CL[j])
                out = integrate.simpson(a, s_beam)
                A[j,i] = out

        print("Done x = ",x_CL[j])
        A1[j,:] = A[j,:]*kappa_CL*L[j]

    A1 *= Q0*N_A
    I_ratio = np.exp(-A1)
    cmap = plt.colormaps["viridis"]
    fig, ax = plt.subplots(1,2, sharey=True)
    for j in range(n_pos):
        color = cmap( float(j)/(n_pos-1) )
        ax[0].plot(wave_range, I_ratio[j,:], color=color, label="{:3.0f}".format(100*x_CL[j,0]))

    for i in range(0,n_wave,4):
        color = cmap( float(i)/(n_wave-1))
        ax[1].plot(x_CL[:,0]*100, I_ratio[:,i], color=color, label="{}".format(wave_range[i]))

    ax[0].set_xlabel("wavelength [nm]")
    ax[1].legend(title="wavelength [nm]")
    ax[0].legend(title="centerline\nposition [cm]")
    ax[1].set_xlabel("centerline position [cm]")
    ax[0].set_ylabel("intensity ratio []")

    return {"kappa_star":A, "kappaL":A1, "intensity_ratio":I_ratio }


if __name__ == "__main__":
    integrate_beams(method={"integrate":"simpsons"})
